# -*- coding: utf-8 -*-
"""Provides functionality for creating and dealing with tabix-indexed files."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import contextlib
import csv
import itertools
import subprocess
from typing import Callable, Iterable, Any

try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib

from future.utils import native_str
import pandas as pd
import pysam


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
                records = itertools.chain.from_iterable((gtf_file.fetch(
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


def index_gtf(file_path, output_path, force=True):
    # type: (pathlib.Path, pathlib.Path) -> None
    """Compresses and indexes a gtf file using bgzip and tabix."""

    # Sort file before compressing and indexing.
    sorted_path = _append_suffix(output_path, '.srt')
    sort_gtf(file_path, output_path=sorted_path)

    # Gzip and index file.
    pysam.tabix_compress(
        str(sorted_path), filename_out=str(output_path), force=force)
    pysam.tabix_index(str(output_path), preset='gff', force=force)

    # Clean up temp file.
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

    gtf = read_gtf(file_path)
    gtf = gtf.sort_values(['seqname', 'start'], ascending=True)
    gtf.to_csv(
        str(output_path),
        sep='\t',
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE)


def _append_suffix(path_obj, suffix):
    # type: (pathlib.Path, str) -> pathlib.Path
    if len(path_obj.suffixes) == 0:
        return path_obj.with_suffix(suffix)
    else:
        return path_obj.with_suffix(path_obj.suffixes[-1] + suffix)


def read_gtf(file_path):
    """Reads an uncompressed GTF file into a pandas DataFrame."""

    names = [
        'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
        'frame', 'attribute'
    ]

    dtypes = {
        'seqname': str,
        'source': 'category',
        'feature': 'category',
        'start': int,
        'end': int,
        'strand': 'category',
        'frame': 'category',
        'attribute': 'str'
    }

    return pd.read_csv(
        str(file_path), sep='\t', names=names, dtype=dtypes, comment='#')
