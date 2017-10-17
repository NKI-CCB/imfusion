# -*- coding: utf-8 -*-
"""Module containing model classes for fusions and insertions."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import collections
import itertools
import operator
import warnings

from future.utils import native_str
import numpy as np
import pandas as pd
import toolz

from imfusion.vendor.frozendict import frozendict


class RecordSet(object):
    """Base class that provides functionality for serializing and
       deserializing namedtuple records into a DataFrame format.

    Subclasses should override the ``_tuple_class`` method to
    return the namedtuple class that should be used as a record.
    """

    def __init__(self, values):
        self._values = self._check_frame(values)

    @classmethod
    def _check_frame(cls, values):
        if len(values) > 0:
            fields = cls._tuple_fields()

            for field in fields:
                if field not in values.columns:
                    raise ValueError(
                        'Missing required column {}'.format(field))

        return values.reindex(columns=fields)

    @classmethod
    def _tuple_class(cls):
        """Returns namedtuple class used to instantiate records."""
        raise NotImplementedError()

    @classmethod
    def _tuple_fields(cls):
        """Returns the fields in the named tuple class."""
        return cls._tuple_class()._fields

    @property
    def values(self):
        """Internal DataFrame representation of records."""
        return self._values

    def __getitem__(self, item):
        return self._values[item]

    def __setitem__(self, idx, value):
        self._values[idx] = value

    def __len__(self):
        return len(self._values)

    @property
    def loc(self):
        """Label-based indexer (similar to pandas .loc)."""
        return LocWrapper(self._values.loc, constructor=self._loc_constructor)

    @property
    def iloc(self):
        """Label-based indexer (similar to pandas .loc)."""
        return LocWrapper(self._values.iloc, constructor=self._loc_constructor)

    def _loc_constructor(self, values):
        if len(values.shape) != 2:
            return values
        return self.__class__(values)

    @classmethod
    def from_tuples(cls, tuples):
        """Builds a record set instance from the given tuples."""
        records = (tup._asdict() for tup in tuples)
        return cls(pd.DataFrame.from_records(records))

    def to_tuples(self):
        """Converts the record set into an iterable of tuples."""

        tuple_class = self._tuple_class()

        for row in self._values.itertuples():
            row_dict = row._asdict()
            row_dict.pop('Index', None)

            yield tuple_class(**row_dict)

    @classmethod
    def from_csv(cls, file_path, **kwargs):
        """Reads a record set from a csv file using pandas.read_csv."""
        values = pd.read_csv(native_str(file_path), **kwargs)
        return cls(values)

    def to_csv(self, file_path, **kwargs):
        """Writes the record set to a csv file using pandas' to_csv."""
        self._values.to_csv(file_path, **kwargs)

    def groupby(self, by, **kwargs):
        """Groups the set by values of the specified columns."""
        for key, group in self._values.groupby(by, **kwargs):
            yield key, self.__class__(group)

    def query(self, expr, **kwargs):
        """Queries the columns of the set with a boolean expression."""
        return self.__class__(self._values.query(expr, **kwargs))

    @classmethod
    def concat(cls, record_sets):
        """Concatenates multiple records sets into a single set."""
        return cls(pd.concat((rs.values for rs in record_sets), axis=0))


class LocWrapper(object):
    """Wrapper class that wraps an objects loc/iloc accessor."""

    def __init__(self, loc, constructor=None):
        if constructor is None:
            constructor = lambda x: x

        self._loc = loc
        self._constructor = constructor

    def __getitem__(self, item):
        result = self._loc[item]
        return self._constructor(result)


class MetadataRecordSet(RecordSet):
    """Base RecordSet that supports record metadata.

    Extension of the RecordSet class, which assumes that records contain
    a dict 'metadata' field which contains variable metadata. The
    MetadataRecordSet class ensures that this data is expanded from the
    original record when converted to the set's DataFrame format, and
    converted back again when transforming back to tuples.
    """

    METADATA_FIELD = 'metadata'

    @property
    def metadata_columns(self):
        """Available metadata columns."""
        return set(self._values.columns) - set(self._tuple_fields())

    @classmethod
    def _check_frame(cls, values):
        fields = [
            field for field in cls._tuple_fields()
            if field != cls.METADATA_FIELD
        ]

        if len(values) > 0:
            for field in fields:
                if field not in values.columns:
                    raise ValueError(
                        'Missing required column {}'.format(field))

        extra_cols = set(values.columns) - set(fields)
        col_order = list(fields) + sorted(extra_cols)

        return values.reindex(columns=col_order)

    @classmethod
    def from_tuples(cls, tuples):
        """Builds a record set instance from the given tuples."""

        metadata_field = cls.METADATA_FIELD

        def _to_record(tup):
            record = tup._asdict()
            record.update(record.pop(metadata_field))
            return record

        records = (_to_record(tup) for tup in tuples)
        return cls(pd.DataFrame.from_records(records))

    def to_tuples(self):
        """Converts the record set into an iterable of tuples."""

        tuple_class = self._tuple_class()

        for row in self._values.itertuples():
            row_dict = row._asdict()
            row_dict.pop('Index', None)

            yield tuple_class(**row_dict)


_Fusion = collections.namedtuple('Fusion', [
    'chromosome_a', 'position_a', 'strand_a', 'chromosome_b', 'position_b',
    'strand_b', 'flank_a', 'flank_b', 'support_junction', 'support_spanning',
    'sample'
])


class Fusion(_Fusion):
    """Model Fusion class."""

    __slots__ = ()

    @property
    def support(self):
        """Support score for the fusion."""
        return self.support_junction + self.support_spanning

    def normalize(self, chromosome=None):
        """Normalizes fusions so that the side whose chromosome has the
           lowest lexical ordering is used as donor.
        """

        if chromosome is not None:
            if self.chromosome_a == chromosome:
                is_norm = True
            else:
                if self.chromosome_b != chromosome:
                    raise ValueError(
                        'Fusion does not include given chromosome')
                is_norm = False
        elif self.chromosome_a != self.chromosome_b:
            is_norm = self.chromosome_a < self.chromosome_b
        else:
            is_norm = self.position_a < self.position_b

        if is_norm:
            return self
        else:
            return self._replace(
                chromosome_a=self.chromosome_b,
                position_a=self.position_b,
                strand_a=self.strand_b * -1,
                chromosome_b=self.chromosome_a,
                position_b=self.position_a,
                strand_b=self.strand_a * -1,
                flank_a=self.flank_b,
                flank_b=self.flank_a)

    def distance(self, other):
        """Determine distance to other fusion."""

        if (self.chromosome_a != other.chromosome_a
                or self.chromosome_b != other.chromosome_b):
            raise ValueError('Fusions are on different reference sequences')

        return (abs(self.position_a - other.position_a) +
                abs(self.position_b - other.position_b))

    @classmethod
    def merge(cls, junctions, max_dist):
        """Merges fusion within genomic distance into a single fusion."""
        for grp in cls._groupby_genomic_position(junctions, max_dist):
            if len(grp) == 1:
                yield grp[0]
            else:
                get_support = operator.attrgetter('support_junction')
                sorted_grp = sorted(grp, key=get_support)

                yield sorted_grp[-1]._replace(
                    support_junction=sum(f.support_junction
                                         for f in sorted_grp),
                    support_spanning=sum(f.support_spanning
                                         for f in sorted_grp),
                    flank_a=max(f.flank_a for f in sorted_grp),
                    flank_b=max(f.flank_b for f in sorted_grp))

    @classmethod
    def _groupby_genomic_position(cls, junctions, max_dist):

        # Group junctions by strand and sequence.
        def _keyfunc(fusion):
            return (fusion.chromosome_a, fusion.strand_a, fusion.chromosome_b,
                    fusion.strand_b)

        sorted_juncs = sorted(junctions, key=_keyfunc)
        grouped_juncs = itertools.groupby(sorted_juncs, key=_keyfunc)

        for _, grp in grouped_juncs:
            for grp_a in cls._groupby_position(grp, 'a', max_dist):
                for grp_b in cls._groupby_position(grp_a, 'b', max_dist):
                    yield grp_b

    @staticmethod
    def _groupby_position(junctions, side, max_dist):
        get_loc = operator.attrgetter('position_' + side)
        sorted_juncs = sorted(junctions, key=get_loc)

        grp, prev_pos = [], np.nan
        for junc in sorted_juncs:
            if (get_loc(junc) - prev_pos) > max_dist:
                yield grp
                grp = [junc]
            else:
                grp.append(junc)
            prev_pos = get_loc(junc)

        if len(grp) > 0:
            yield grp


_TransposonFusion = collections.namedtuple('TransposonFusion', [
    'chromosome', 'anchor_genome', 'anchor_transposon', 'strand_genome',
    'strand_transposon', 'flank_genome', 'flank_transposon',
    'support_junction', 'support_spanning', 'sample', 'metadata'
])


class TransposonFusion(_TransposonFusion):
    """Model class representing a gene-transposon fusion."""

    __slots__ = ()

    def __new__(cls, chromosome, anchor_genome, anchor_transposon,
                strand_genome, strand_transposon, flank_genome,
                flank_transposon, support_junction, support_spanning, sample,
                **kwargs):

        return super().__new__(
            cls,
            chromosome=chromosome,
            anchor_genome=int(anchor_genome),
            anchor_transposon=int(anchor_transposon),
            strand_genome=int(strand_genome),
            strand_transposon=int(strand_transposon),
            flank_genome=int(flank_genome),
            flank_transposon=int(flank_transposon),
            support_junction=int(support_junction),
            support_spanning=int(support_spanning),
            sample=sample,
            metadata=frozendict(kwargs))

    def __getattr__(self, name):
        return self.metadata[name]

    @property
    def support(self):
        """Support score for the fusion."""
        return self.support_junction + self.support_spanning

    @property
    def genome_region(self):
        """Returns the spanned region on the genome."""
        if self.flank_genome > 0:
            start = self.anchor_genome
            end = self.anchor_genome + self.flank_genome
        else:
            start = self.anchor_genome + self.flank_genome
            end = self.anchor_genome
        return self.chromosome, start, end

    @property
    def transposon_region(self):
        """Returns the spanned region on the transposon."""
        if self.flank_transposon > 0:
            start = self.anchor_transposon
            end = self.anchor_transposon + self.flank_transposon
        else:
            start = self.anchor_transposon + self.flank_transposon
            end = self.anchor_transposon
        return start, end

    @classmethod
    def from_fusion(cls, fusion, transposon_name, metadata=None):
        """Converts the fusion to a transposon fusion object."""

        if (fusion.chromosome_a == transposon_name
                and fusion.chromosome_b == transposon_name):
            raise ValueError('Fusion does not involve genomic sequence')
        elif (fusion.chromosome_a != transposon_name
              and fusion.chromosome_b != transposon_name):
            raise ValueError('Fusion does not involve transposon')
        elif fusion.chromosome_a == transposon_name:
            tr_key, gen_key = 'a', 'b'
            tr_flank = fusion.flank_a * -fusion.strand_a
            gen_flank = fusion.flank_b * fusion.strand_b
        else:
            tr_key, gen_key = 'b', 'a'
            tr_flank = fusion.flank_b * fusion.strand_b
            gen_flank = fusion.flank_a * -fusion.strand_a

        return cls(
            chromosome=getattr(fusion, 'chromosome_' + gen_key),
            anchor_genome=getattr(fusion, 'position_' + gen_key),
            anchor_transposon=getattr(fusion, 'position_' + tr_key),
            strand_genome=getattr(fusion, 'strand_' + gen_key),
            strand_transposon=getattr(fusion, 'strand_' + tr_key),
            flank_genome=gen_flank,
            flank_transposon=tr_flank,
            support_junction=fusion.support_junction,
            support_spanning=fusion.support_spanning,
            sample=fusion.sample,
            **(metadata or {}))


_Insertion = collections.namedtuple('Insertion', [
    'id', 'chromosome', 'position', 'strand', 'support', 'sample', 'metadata'
])


class Insertion(_Insertion):
    """Model class representing an insertion."""

    __slots__ = ()

    def __new__(cls, id, chromosome, position, strand, support, sample,
                **kwargs):

        return super().__new__(
            cls,
            id=id,
            chromosome=chromosome,
            position=int(position),
            strand=int(strand),
            support=int(support),
            sample=sample,
            metadata=frozendict(kwargs))

    def __getattr__(self, name):
        return self.metadata[name]

    @classmethod
    def from_transposon_fusion(cls, fusion, id=None, drop_metadata=None):
        """Converts (annotated) transposon fusion to an insertion.

        Parameters
        ----------
        id : str
            ID to use for the insertion.

        Returns
        -------
        Insertion
            Insertion object derived from the transposon fusion.

        """

        if drop_metadata is None:
            drop_metadata = {'strand_genome', 'strand_transposon'}

        # Fusion should have annotated transposon feature.
        if 'feature_name' not in fusion.metadata:
            raise ValueError(
                'Fusion does not have an annotated transposon feature')

        strand = fusion.strand_genome * fusion.strand_transposon

        if 'gene_strand' in fusion.metadata:
            gene_strand = fusion.metadata['gene_strand']
            orientation = 'sense' if strand == gene_strand else 'antisense'
        else:
            orientation = None

        ins_metadata = toolz.keyfilter(lambda k: k not in drop_metadata,
                                       fusion.metadata)

        ins_metadata['transposon_anchor'] = fusion.anchor_transposon
        ins_metadata['support_junction'] = fusion.support_junction
        ins_metadata['support_spanning'] = fusion.support_spanning

        if orientation is not None:
            ins_metadata['orientation'] = orientation

        return Insertion(
            id=id,
            chromosome=fusion.chromosome,
            position=fusion.anchor_genome,
            strand=strand,
            support=fusion.support,
            sample=fusion.sample,
            **ins_metadata)

    @classmethod
    def from_transposon_fusions(cls,
                                fusions,
                                id_fmt_str=None,
                                drop_metadata=None):
        """Converts annotated transposon fusions to insertions."""

        if id_fmt_str is not None:
            for i, fusion in enumerate(fusions):
                yield cls.from_transposon_fusion(
                    fusion,
                    id=id_fmt_str.format(i + 1, sample=fusion.sample),
                    drop_metadata=drop_metadata)
        else:
            for fusion in fusions:
                yield cls.from_transposon_fusion(
                    fusion, drop_metadata=drop_metadata)


class InsertionSet(MetadataRecordSet):
    """Class that represents an insertion dataset."""

    @property
    def samples(self):
        return set(self._values['sample'])

    @classmethod
    def _tuple_class(cls):
        return Insertion

    @classmethod
    def from_csv(cls, file_path, **kwargs):
        """Reads a InsertionSet from a csv file using pandas.read_csv."""
        values = pd.read_csv(native_str(file_path), **kwargs)

        if 'chromosome' not in values and 'seqname' in values:
            warnings.warn('Using \'seqname\' to describe insertion positions '
                          'is deprecated, use \'chromosome\' instead',
                          DeprecationWarning)
            values = values.rename(columns={'seqname': 'chromosome'})

        return cls(values)

    @classmethod
    def concat(cls, insertion_sets, check_overlap=True):
        """Concatenates multiple records sets into a single set."""

        seen = set()
        if check_overlap:
            for ins_set in insertion_sets:
                if len(ins_set) > 0 and (seen & ins_set.samples):
                    raise ValueError('Duplicate samples in insertion sets')
                seen |= ins_set.samples

        return super().concat(insertion_sets)
