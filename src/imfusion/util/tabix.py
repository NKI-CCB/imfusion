# -*- coding: utf-8 -*-
"""Provides functionality for creating and dealing with tabix-indexed files."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import contextlib
import csv
from itertools import chain, groupby, tee
import subprocess
from typing import Callable, Iterable, Any

from future.utils import native_str
import pandas as pd
import pathlib2 as pathlib
import pysam
from sortedcontainers import SortedList

GTF_COLUMNS = [
    'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
    'attribute'
]

GTF_DTYPES = {
    'seqname': str,
    'source': 'category',
    'feature': 'category',
    'start': int,
    'end': int,
    'strand': 'category',
    'frame': 'category',
    'attribute': 'str'
}


class GtfIterator(object):
    """Iterator that iterates over records in a GTF file using pysam."""

    def __init__(self, gtf_path):
        # type: (pathlib.Path) -> None
        self._gtf_path = gtf_path

    @contextlib.contextmanager
    def _open_file(self):  # type: (...) -> pysam.TabixFile
        # Open gtf file.
        gtf_file = pysam.TabixFile(
            native_str(self._gtf_path), parser=pysam.asGTF())

        # Yield file object and ensure it is closed.
        try:
            yield gtf_file
        finally:
            gtf_file.close()

    @property
    def contigs(self):  # type: (...) -> List[str]
        """Contigs present in the gtf file."""

        with self._open_file() as gtf_file:
            return gtf_file.contigs

    def fetch(
            self,
            reference=None,  # type: str
            start=None,  # type: int
            end=None,  # type: int
            filters=None  # type: Iterable[Callable]
    ):  # type (...) -> Iterable[pysam.ctabixproxies.GTFProxy]
        """Fetches gtf records (gtf rows) from the gtf file."""

        with self._open_file() as gtf_file:
            # For some reason pysam does not fetch all records if reference
            # is None under Python 2.7. To fix this, here we simply chain all
            # the contig records into one iterable.
            if reference is None:
                contigs = gtf_file.contigs
                records = chain.from_iterable((gtf_file.fetch(
                    reference=ref, start=start, end=end) for ref in contigs))
            else:
                records = gtf_file.fetch(
                    reference=reference, start=start, end=end)

            # Filter records on additional filters.
            if filters is not None:
                for filter_func in filters:
                    records = (rec for rec in records if filter_func(rec))

            for record in records:
                yield record

    def fetch_genes(self, gene_ids=None, **kwargs):
        # type: (List[str], **Any) -> Iterable[pysam.ctabixproxies.GTFProxy]
        """Fetches gene records from the GTF file."""

        # Fetch gene records.
        gene_filter = lambda rec: rec.feature == 'gene'
        records = self.fetch(filters=[gene_filter], **kwargs)

        if gene_ids is not None:
            # Run through records.
            sought_ids = set(gene_ids)

            hits = {}
            for rec in records:
                if rec['gene_id'] in sought_ids:
                    # Record hits.
                    hits[rec['gene_id']] = rec

                    # Stop if we found everything.
                    if len(hits) == len(gene_ids):
                        break

            # Yield results, returning None for genes that weren't found.
            for gene_id in gene_ids:
                yield hits.get(gene_id, None)

        else:
            for record in records:
                yield record


def index_gtf(file_path, output_path, sort=True, force=True):
    # type: (pathlib.Path, pathlib.Path) -> None
    """Compresses and indexes a gtf file using bgzip and tabix."""

    # Sort file before compressing and indexing.
    if sort:
        sorted_path = _append_suffix(output_path, '.srt')
        sort_gtf(file_path, output_path=sorted_path)
    else:
        sorted_path = file_path

    # Gzip and index file.
    pysam.tabix_compress(
        native_str(sorted_path),
        filename_out=native_str(output_path),
        force=force)
    pysam.tabix_index(native_str(output_path), preset='gff', force=force)

    # Clean up sort temp file.
    if sort:
        sorted_path.unlink()


def sort_gtf(file_path, output_path):
    # type: (pathlib.Path, pathlib.Path) -> None
    """Sorts a gtf file by position, required for indexing by tabix."""

    try:
        # Try sorting using shell command.
        _sort_gtf_shell(file_path, output_path)
    except subprocess.CalledProcessError:
        # Fall back into in-memory sort if shell command fails.
        # TODO: Log warning.
        _sort_gtf_pandas(file_path, output_path)


def _sort_gtf_shell(file_path, output_path):
    """Sorts gtf file externally using grep/sort shell command."""

    with open(str(output_path), 'w') as out_file:
        cmd = '(grep ^"#" {0}; grep -v ^"#" {0} | sort -k1,1 -k4,4n)'
        subprocess.check_call(
            cmd.format(file_path), stdout=out_file, shell=True)


def _sort_gtf_pandas(file_path, output_path):
    """Sorts gtf file in-memory using pandas."""

    gtf_frame = read_gtf_frame(file_path)
    gtf_frame = sort_gtf_frame(gtf_frame)
    write_gtf_frame(gtf_frame, output_path)


def _append_suffix(path_obj, suffix):
    # type: (pathlib.Path, str) -> pathlib.Path
    if len(path_obj.suffixes) == 0:
        return path_obj.with_suffix(suffix)
    else:
        return path_obj.with_suffix(path_obj.suffixes[-1] + suffix)


def read_gtf_frame(file_path):
    """Reads an uncompressed GTF file into a pandas DataFrame."""

    return pd.read_csv(
        str(file_path),
        sep='\t',
        names=GTF_COLUMNS,
        dtype=GTF_DTYPES,
        comment='#')


def write_gtf_frame(gtf_frame, file_path):
    """Writes gtf frame to a GTF file."""

    gtf_frame.to_csv(
        str(file_path),
        sep='\t',
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE)


def sort_gtf_frame(gtf_frame):
    """Sorts gtf data by position."""
    return gtf_frame.sort_values(['seqname', 'start'], ascending=True)


def flatten_gtf_frame(gtf_frame):
    """Flattens exons in given GTF in the same style as DEXSeq."""

    # Extract exons and add gene_id attribute to group by.
    gtf_frame = gtf_frame.loc[gtf_frame['feature'] == 'exon'].copy()
    gtf_frame['gene_id'] = gtf_frame['attribute'].str.extract(
        r'gene_id "(\w+)"', expand=False)

    # Flatten per gene and merge into single frame.
    flattened = pd.concat(
        (_flatten_gene_group(grp) for _, grp in gtf_frame.groupby('gene_id')),
        axis=0,
        ignore_index=True)

    return flattened


def _flatten_gene_group(grp):
    """Flattens exons for a given gene."""

    first = grp.iloc[0]

    # Flatten intervals.
    grp = grp.sort_values('start', ascending=True)
    intervals = list(zip(grp.start, grp.end))
    flattened = list(_flatten_intervals(intervals))

    # Reverse if on negative strand.
    if first.strand == '-':
        flattened = flattened[::-1]

    # Extract starts and ends.
    flat_starts, flat_ends = zip(*flattened)

    # Build attribute strings.
    attr_fmt = 'gene_id "' + first.gene_id + '"; exonic_part_number "{:03d}"'
    attributes = [attr_fmt.format(i) for i in range(1, len(flattened) + 1)]

    # Build frame.
    return pd.DataFrame(
        {
            'seqname': first.seqname,
            'source': 'flatten_gtf',
            'feature': 'exonic_part',
            'start': flat_starts,
            'end': flat_ends,
            'score': '.',
            'strand': first.strand,
            'frame': '.',
            'attribute': attributes
        },
        columns=GTF_COLUMNS)


def _flatten_intervals(intervals):
    """Flattens a list of (start, end) intervals."""

    # Remove duplicate intervals and sort.
    intervals = sorted(set(intervals))

    if len(intervals) > 0:
        # Start first run.
        start, end = intervals[0]
        positions = SortedList((start - 1, end))

        for start, end in intervals[1:]:
            if positions[-1] < start:
                # We have left the current run, yield everything in positions.
                for interval in _intervals_from_positions(positions):
                    yield interval

                # Start new run.
                positions = SortedList((start - 1, end))
            else:
                # Record current end as part of the current run.
                positions.update((start - 1, end))

        # Yield any remaining intervals.
        for interval in _intervals_from_positions(positions):
            yield interval


def _intervals_from_positions(positions):
    """Returns intervals from positions."""

    # Remove duplicates and correct first position.
    pos = list(_uniq(positions))

    # Yield intervals.
    for start, end in _pairwise(_uniq(pos)):
        yield start + 1, end


def _uniq(iterable):
    """Returns uniq elements for a sorted iterable.

    Example: [1, 1, 2, 2, 2, 3, 3, 1] -> [1, 2, 3, 1

    """

    for k, _ in groupby(iterable):
        yield k


def _pairwise(iterable):
    """Returns pairs of current, next values from list.

    Example: s -> (s0,s1), (s1,s2), (s2, s3), ...

    """

    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)
