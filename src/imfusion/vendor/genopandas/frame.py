# -*- coding: utf-8 -*-
"""Dataframe-related functions/classes."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from collections import OrderedDict
import math

import numpy as np
import pandas as pd

from .tree import GenomicIntervalTree
from .util import reorder_columns


class GenomicDataFrame(pd.DataFrame):
    """DataFrame with fast indexing by genomic position.

    Requires columns 'chromosome', 'start' and 'end' to be present in the
    DataFrame, as these columns are used for indexing.

    Examples
    --------

    Constructing from scratch:

    >>> df = pd.DataFrame.from_records(
    ...    [('1', 10, 20), ('2', 10, 20), ('2', 30, 40)],
    ...    columns=['chromosome', 'start', 'end'])
    >>> GenomicDataFrame(df)

    Constructing with non-default columns:

    >>> df = pd.DataFrame.from_records(
    ...    [('1', 10, 20), ('2', 10, 20), ('2', 30, 40)],
    ...    columns=['chrom', 'chromStart', 'chromEnd'])
    >>> GenomicDataFrame(
    ...    df,
    ...    chromosome_col='chrom',
    ...    start_col='start',
    ...    end_col='end')

    Reading from a GTF file:

    >>> GenomicDataFrame.from_gtf('/path/to/reference.gtf.gz')

    Querying by genomic position:

    >>> genomic_df.gi.search('2', 30, 50)

    """

    _internal_names = pd.DataFrame._internal_names + ['_gi']
    _internal_names_set = set(_internal_names)

    _metadata = ['_gi_metadata']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        chrom_lengths = kwargs.pop('chrom_lengths', None)

        self._gloc = None
        self._gloc_metadata = {'chrom_lengths': chrom_lengths}

    @property
    def gloc(self):
        """Genomic indexer for querying the dataframe."""
        if self._gloc is None:
            self._gloc = GenomicIndexer.from_df(self, **self._gloc_metadata)
        return self._gloc

    @property
    def _constructor(self):
        return GenomicDataFrame

    @classmethod
    def from_csv(cls,
                 file_path,
                 index_col=None,
                 drop_index_col=True,
                 chrom_lengths=None,
                 **kwargs):
        """Creates a GenomicDataFrame from a csv file using ``pandas.read_csv``.

        Parameters
        ----------
        file_path : str
            Path to file.
        index_col : List[str] or List[int]
            Columns to use for index. Columns can be indicated by their name
            (str) or index (int). Should contain two entries for positioned
            data, three entries for ranged data. If not given, the first three
            columns of the DataFrame are used by default.
        drop_index_col : bool
            Whether to drop the index columns in the DataFrame (True, default)
            or to drop them from the dataframe (False).
        chrom_lengths : Dict[str, int]
            Chromosomes sizes to use in the genomic index.
        **kwargs
            Any extra kwargs are passed to ``pandas.read_csv``.

        Returns
        -------
        GenomicDataFrame
            DataFrame containing the file contents.

        """

        data = pd.read_csv(file_path, index_col=None, **kwargs)

        if index_col is None:
            index_col = list(data.columns[:3])
        elif isinstance(index_col[0], int):
            index_col = [data.columns[i] for i in index_col]

        # Convert chromosome to str.
        data[index_col[0]] = data[index_col[0]].astype(str)

        data = data.set_index(index_col, drop=drop_index_col)

        return cls(data, chrom_lengths=chrom_lengths)

    @classmethod
    def from_gtf(cls, gtf_path, filter_=None):
        """Creates a GenomicDataFrame from a GTF file."""

        try:
            import pysam
        except ImportError:
            raise ImportError('Pysam needs to be installed for '
                              'reading GTF files')

        # Parse records into rows.
        gtf_file = pysam.TabixFile(str(gtf_path), parser=pysam.asGTF())
        records = (rec for rec in gtf_file.fetch())

        # Filter records if needed.
        if filter_ is not None:
            records = (rec for rec in records if filter_(rec))

        # Build dataframe.
        def _record_to_row(record):
            row = {
                'contig': record.contig,
                'source': record.source,
                'feature': record.feature,
                'start': record.start,
                'end': record.end,
                'score': record.score,
                'strand': record.strand,
                'frame': record.frame
            }
            row.update(dict(record))
            return row

        gdf = cls.from_records(
            (_record_to_row(rec) for rec in records),
            index_col=['contig', 'start', 'end'],
            drop_index_col=False)

        # Reorder columns to correspond with GTF format.
        columns = ('contig', 'source', 'feature', 'start', 'end', 'score',
                   'strand', 'frame')
        gdf = reorder_columns(gdf, order=columns)

        return gdf

    @classmethod
    def from_records(cls,
                     data,
                     index_col,
                     columns=None,
                     drop_index_col=True,
                     **kwargs):
        """Creates a GenomicDataFrame from a structured or record ndarray."""

        if not 2 <= len(index_col) <= 3:
            raise ValueError('index_col should contain 2 entries'
                             ' (for positioned data or 3 entries'
                             ' (for ranged data)')

        df = super().from_records(data, columns=columns, **kwargs)
        df = df.set_index(index_col, drop=drop_index_col)

        return df

    @property
    def is_positioned(self):
        """Returns true if dataframe is positioned."""
        return self.index.nlevels == 2

    @property
    def is_ranged(self):
        """Returns true if dataframe is ranged."""
        return self.index.nlevels == 3

    def as_positioned(self):
        """Converts a ranged frame (with starts/ends in the index) to a
           positioned df (with single positions in the index)."""

        if self.is_positioned:
            return self.copy()

        chromosomes = self.index.get_level_values(0)
        starts = self.index.get_level_values(1)
        ends = self.index.get_level_values(2)

        positions = (starts + ends) // 2

        names = [self.index.names[0], 'position']
        new_index = pd.MultiIndex.from_arrays(
            [chromosomes, positions], names=names)

        new_df = self.copy()
        new_df.index = new_index

        return new_df

    def as_ranged(self, width=1):
        """Converts a positioned frame (with positions in the index)
           to a ranged df (with starts/ends in the index)."""

        if not self.is_positioned:
            raise ValueError('Only positioned frames can be '
                             'converted to a ranged frame')

        chromosomes = self.index.get_level_values(0)
        positions = self.index.get_level_values(1)

        starts = positions - (width // 2)
        ends = positions + math.ceil(width / 2)

        names = [self.index.names[0], 'start', 'end']
        new_index = pd.MultiIndex.from_arrays(
            [chromosomes, starts, ends], names=names)

        new_df = self.copy()
        new_df.index = new_index

        return new_df


class GenomicIndexer(object):
    """Base GenomicIndexer class used to index GenomicDataFrames."""

    def __init__(self, df, chrom_lengths=None):
        self._df = df
        self._lengths = chrom_lengths

        self._trees = None

    def __getitem__(self, item):
        """Accessor used to query the dataframe by position.

        If a list of chromosomes is given, the dataframe is subset to the
        given chromosomes. Note that chromosomes are also re-ordered to
        adhere to the given order. If a single chromosome is given, a
        GenomicSlice is returned. This slice object can be sliced to query
        a specific genomic range.
        """

        if isinstance(item, list):
            return self._df.reindex(index=[item], level=0)

        return GenomicSlice(self, chromosome=item)

    @classmethod
    def from_df(cls, df, chrom_lengths=None):
        """Constructs the appropriate indexer sub-class for the given frame."""

        if df.index.nlevels == 3:
            return GenomicRangeIndexer(df, chrom_lengths=chrom_lengths)
        elif df.index.nlevels == 2:
            return GenomicPositionIndexer(df, chrom_lengths=chrom_lengths)
        else:
            raise ValueError('DataFrame index should contain two levels '
                             '(for positional data) or three levels '
                             '(for ranged data)')

    @property
    def df(self):
        """The indexed DataFrame."""
        return self._df

    @property
    def chromosome(self):
        """Chromosome values."""
        return self._df.index.get_level_values(0)

    @property
    def chromosomes(self):
        """Available chromosomes."""
        return list(self.chromosome.unique())

    @property
    def chromosome_lengths(self):
        """Chromosome lengths."""

        if self._lengths is None:
            self._lengths = self._calculate_lengths()

        return self._lengths

    def _calculate_lengths(self):
        raise NotImplementedError()

    @property
    def chromosome_offsets(self):
        """Chromosome offsets (used when plotting chromosomes linearly)."""

        # Sort lengths by chromosome.
        chromosomes = self.chromosomes
        lengths = self.chromosome_lengths

        # Record offsets in ordered dict.
        sorted_lengths = [lengths[chrom] for chrom in chromosomes]

        cumsums = np.concatenate([[0], np.cumsum(sorted_lengths)])
        offsets = OrderedDict(zip(chromosomes, cumsums[:-1]))

        # Add special marker for end.
        offsets['_END_'] = cumsums[-1]

        return offsets

    def _offset_positions(self, positions):
        offsets = pd.Series(self.chromosome_offsets)
        return positions + offsets.loc[self.chromosome].values

    @property
    def trees(self):
        """Trees used for indexing the DataFrame."""

        if self._trees is None:
            self._trees = self._build_trees()

        return self._trees

    def rebuild(self):
        """Rebuilds the genomic interval trees."""
        self._trees = self._build_trees()

    def _build_trees(self):
        raise NotImplementedError()

    def search(self,
               chromosome,
               start,
               end,
               strict_left=False,
               strict_right=False):
        """Searches the DataFrame for rows within given range."""

        overlap = self.trees.search(
            chromosome,
            start,
            end,
            strict_left=strict_left,
            strict_right=strict_right)

        indices = [interval[2] for interval in overlap]

        return self._df.iloc[indices].sort_index()


class GenomicSlice(object):
    """Supporting class used by the GenomicIndexer for slicing chromosomes."""

    def __init__(self, indexer, chromosome):
        self._indexer = indexer
        self._chromosome = chromosome

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self._indexer.search(
                self._chromosome, start=item.start, end=item.stop)
        return self._indexer.search(self._chromosome, start=item)


class GenomicRangeIndexer(GenomicIndexer):
    """GenomicIndexer class for querying ranged (start/end) data."""

    def __init__(self, df, chrom_lengths=None):

        if not df.index.nlevels == 3:
            raise ValueError('Dataframe must have three levels')

        super().__init__(df, chrom_lengths=chrom_lengths)

    @property
    def start(self):
        """Start positions."""
        return self._df.index.get_level_values(1)

    @property
    def start_offset(self):
        """Start positions, offset by chromosome lengths."""
        return self._offset_positions(self.start)

    @property
    def end(self):
        """End positions."""
        return self._df.index.get_level_values(2)

    @property
    def end_offset(self):
        """End positions, offset by chromosome lengths."""
        return self._offset_positions(self.end)

    def _calculate_lengths(self):
        lengths = pd.Series(self.end).groupby(self.chromosome).max()
        return dict(zip(lengths.index, lengths.values))

    def _build_trees(self):
        tuples = zip(self.chromosome, self.start, self.end,
                     range(self._df.shape[0]))
        return GenomicIntervalTree.from_tuples(tuples)


class GenomicPositionIndexer(GenomicIndexer):
    """GenomicIndexer class for querying positioned (single position) data."""

    def __init__(self, df, chrom_lengths=None):

        if not df.index.nlevels == 2:
            raise ValueError('Dataframe must have two levels')

        super().__init__(df, chrom_lengths=chrom_lengths)

    @property
    def position(self):
        """Positions."""
        return self._df.index.get_level_values(1)

    @property
    def position_offset(self):
        """Positions, offset by chromosome lengths."""
        return self._offset_positions(self.position)

    def _calculate_lengths(self):
        lengths = pd.Series(self.position).groupby(self.chromosome).max()
        return dict(zip(lengths.index, lengths.values))

    def _build_trees(self):
        positions = self.position
        tuples = zip(self.chromosome, positions, positions + 1,
                     range(self._df.shape[0]))
        return GenomicIntervalTree.from_tuples(tuples)
