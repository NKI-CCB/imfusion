# -*- coding: utf-8 -*-
"""Provides functionality for creating and dealing with tabix-indexed files."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import contextlib
import itertools
import subprocess
from typing import Callable, Iterable, Any

try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib

from future.utils import native_str
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


def index_gtf(file_path, output_path):
    # type: (pathlib.Path, pathlib.Path) -> None
    """Compresses and indexes a gtf file using bgzip and tabix."""

    # Sort file before compressing and indexing.
    sorted_path = _append_suffix(output_path, '.srt')
    sort_gtf(file_path, output_path=sorted_path)

    # Gzip and index file.
    bgzip(sorted_path, output_path=output_path)
    tabix(output_path, preset='gff')

    # Clean up temp file.
    sorted_path.unlink()


def sort_gtf(file_path, output_path):
    # type: (pathlib.Path, pathlib.Path) -> None
    """Sorts a gtf file by position, required for indexing by tabix."""
    with open(str(output_path), 'w') as out_file:
        cmd = '(grep ^"#" {0}; grep -v ^"#" {0} | sort -k1,1 -k4,4n)'
        subprocess.check_call(
            cmd.format(file_path), stdout=out_file, shell=True)


def _append_suffix(path_obj, suffix):
    # type: (pathlib.Path, str) -> pathlib.Path
    if len(path_obj.suffixes) == 0:
        return path_obj.with_suffix(suffix)
    else:
        return path_obj.with_suffix(path_obj.suffixes[-1] + suffix)


def bgzip(file_path, output_path=None):
    # type: (pathlib.Path, pathlib.Path) -> None
    """Uses bgzip to compress given file."""

    output_path = str(output_path)

    if output_path is None:
        output_path = file_path + '.gz'

    with open(output_path, 'w') as out_file:
        subprocess.check_call(['bgzip', '-c', str(file_path)], stdout=out_file)


def tabix(file_path, preset):
    # type: (pathlib.Path, str) -> None
    """Uses tabix to index given file."""
    subprocess.check_call(['tabix', '-p', preset, str(file_path)])
