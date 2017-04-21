# -*- coding: utf-8 -*-
"""Implements functions for testing for Commonly Targeted Genes (CTGs)."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from typing import Tuple, Optional, Iterable, Callable, Pattern
import itertools
import logging
import operator
import re

import pathlib2 as pathlib

import numpy as np
import pandas as pd
import pyfaidx
import pysam

from intervaltree import IntervalTree
from scipy.stats import poisson

from imfusion.build import Reference
from imfusion.model import Insertion
from imfusion.util.genomic import GenomicIntervalTree
from imfusion.util.tabix import GtfIterator


def test_ctgs(
        insertions,  # type: List[Insertion]
        reference,  # type: Reference
        gene_ids=None,  # type: Set[str]
        chromosomes=None,  # type: Set[str]
        pattern=None,  # type: str
        per_sample=True,  # type: bool
        window=None  #type: Tuple[int, int]
):
    """Identifies genes that are significantly enriched for insertions (CTGs).

    This function takes a DataFrame of insertions, coming from multiple samples,
    and identifies if any genes are more frequently affected by an
    insertion than would be expected by chance. These genes are called
    Commonly Targeted Genes (CTGs). CTGs are selected by comparing the
    number of insertions within the gene to the number of insertions
    that would be expected from the background insertion rate, which is
    modeled using a Poisson distribution.

    Parameters
    ----------
    insertions : List[Insertion]
        Insertions to test.
    reference : Reference
        Reference index used by the aligner to identify insertions.
    genes : List[str]
        List of genes to test (defaults to all genes with an insertion).
    chromosomes : List[str]
        List of chromosomes to include, defaults to all chromosomes
        in the given reference GTF.
    pattern : str
        Specificity pattern of the used transposon.
    per_sample : bool
        Whether to perform the per sample test (recommended), which
        effectively collapes insertions per sample/gene combination.
        This avoids issues in which insertions that are detected
        multiple times or that may have hopped inside the gene locus
        are counted multiple times.
    window : Tuple[int, int]
        Window to include around gene (in bp). Specified as (upstream_dist,
        downstream_dist). For example: (2000, 2000) specifies in a 2KB
        window around each gene.

    Returns
    -------
    pandas.DataFrame
        Results of CTG test for tested genes. Contains two columns:
        p_value and q_value. The last column (q_value)
        represents the p-value of the gene after correcting for
        multiple testing using bonferroni correction.

    """

    # Determine gene windows using GTF.
    logging.info('Generating gene windows')
    gene_windows = _build_gene_windows(
        reference.indexed_gtf_path, window=window, chromosomes=chromosomes)

    # Subset insertions to gene intervals.
    insertions = _subset_to_windows(insertions, gene_windows)

    if gene_ids is None:
        gene_ids = set(ins.metadata['gene_id'] for ins in insertions)

    # Collapse insertions per gene/sample (recommended).
    # Corrects for hopping/multiple detection issues.
    if per_sample:
        logging.info('Collapsing insertions')
        insertions = list(_collapse_per_sample(insertions))

    # Calculate total number of pattern occurrences within intervals.
    logging.info('Counting pattern occurrences')
    reference_seq = pyfaidx.Fasta(str(reference.fasta_path))

    total = count_total(
        reference_seq, pattern=pattern, intervals=gene_windows.values())

    # Calculate p-values for each gene.
    logging.info('Calculating significance for genes')
    insertion_trees = GenomicIntervalTree.from_objects_position(
        insertions, chrom_attr='seqname')

    p_values = {
        gene_id: test_region(
            insertions=insertions,
            reference_seq=reference_seq,
            region=gene_windows[gene_id],
            total=total,
            pattern=pattern,
            filters=[lambda ins, gid=gene_id: ins.metadata['gene_id'] == gid],
            insertion_trees=insertion_trees)
        for gene_id in gene_ids
    }

    # Build result frame.
    result = pd.DataFrame.from_records(
        iter(p_values.items()), columns=['gene_id', 'p_value'])

    # Calculate corrected p-value using bonferroni correction.
    result['q_value'] = (result['p_value'] * len(result)).clip_upper(1.0)

    # Sort by q-value and p-value.
    result.sort_values(by=['q_value', 'p_value'], inplace=True)

    # Annotate with gene_name if possible.
    if 'gene_name' in insertions[0].metadata:
        name_map = {
            ins.metadata['gene_id']: ins.metadata['gene_name']
            for ins in insertions
        }
        result.insert(1, 'gene_name', result['gene_id'].map(name_map))

    # Annotate with frequency.
    frequency = (Insertion.to_frame(insertions)
                 .groupby('gene_id')['sample'].nunique()
                 .reset_index(name='n_samples'))
    result = pd.merge(result, frequency, on='gene_id', how='left')

    return result


def _build_gene_windows(
        gtf_path,  # type: pathlib.Path
        window=None,  # type: Optional[Tuple[int, int]]
        chromosomes=None  # type: Set[str]
):
    gtf_iter = GtfIterator(gtf_path)

    if chromosomes is None:
        chromosomes = set(gtf_iter.contigs)

    records = itertools.chain.from_iterable(
        gtf_iter.fetch_genes(reference=chrom) for chrom in chromosomes)

    return {rec['gene_id']: _apply_gene_window(rec, window) for rec in records}


def _apply_gene_window(
        gene,  # type: pysam.ctabixproxies.GTFProxy
        window=None  # type: Tuple[int, int]
):  # type: (...) -> Tuple[str, int, int]

    if window is None:
        return gene.contig, gene.start, gene.end,
    else:
        upstream_offset, downstream_offset = window

        if gene.strand == '-':
            start = gene.start - downstream_offset
            end = gene.end + upstream_offset
        elif gene.strand == '+':
            start = gene.start - upstream_offset
            end = gene.end + downstream_offset
        else:
            raise ValueError('Unknown value for strand')

        return gene.contig, start, end


def _subset_to_windows(
        insertions,  # type: List[Insertion]
        gene_windows  # type: Dict[str, Tuple[str, int, int]]
):  # type: (...) -> List[Insertion]
    """Subsets insertions for given gene windows."""

    # Create lookup trees.
    trees = {
        chrom: IntervalTree.from_tuples((i[1:]) for i in chrom_int)
        for chrom, chrom_int in itertools.groupby(
            sorted(gene_windows.values()), operator.itemgetter(0))
    }

    # Determine which insertions overlap tree intervals and
    # correspond to genes with known gene window.
    def _in_windows(ins, trees):
        try:
            return trees[ins.seqname].overlaps(ins.position)
        except KeyError:
            return False

    return [
        ins for ins in insertions
        if ins.metadata['gene_id'] in gene_windows and _in_windows(ins, trees)
    ]


def _collapse_per_sample(insertions):
    # Type: (List[Insertion]) -> Generator
    def _keyfunc(insertion):
        return (insertion.metadata['sample'],
                str(insertion.metadata['gene_id']))

    grouped = itertools.groupby(sorted(insertions, key=_keyfunc), key=_keyfunc)

    for _, grp in grouped:
        grp = list(grp)

        if len(grp) > 1:
            mean_pos = int(np.mean([ins.position for ins in grp]))
            yield grp[0]._replace(position=mean_pos)
        else:
            yield grp[0]


def test_region(
        insertions,  # type: List[Insertion]
        reference_seq,  # type: pyfaidx.Fasta
        region,  # type: Tuple[str, int, int]
        pattern=None,  # type: Optional[str]
        intervals=None,  # type: Optional[Iterable[Tuple[str, int, int]]]
        total=None,  # type: Optional[int]
        filters=None,  # type: Optional[List[Callable]]
        insertion_trees=None  # type: GenomicIntervalTree
):  # type: (...) -> float
    """Tests a given genomic region for enrichment in insertions."""

    if total is None:
        total = count_total(
            reference_seq, pattern=pattern, intervals=intervals)

    # Count pattern in region.
    region_count = count_region(reference_seq, region=region, pattern=pattern)

    # Sub-select insertions for region.
    if insertion_trees is None:
        insertion_trees = GenomicIntervalTree.from_objects_position(
            insertions, chrom_attr='seqname')

    region_ins = set(interval[2]
                     for interval in insertion_trees.search(*region))

    # Apply additional filter functions to insertions if given
    # (such as filtering on gene name/id for example).
    if filters is not None:
        for filter_func in filters:
            region_ins = set(ins for ins in region_ins if filter_func(ins))

    # Calculate p-value.
    x = len(list(region_ins))
    mu = len(insertions) * (region_count / total)

    # Note here we use loc=1, because we are interested in
    # calculating P(X >= x), not P(X > x) (the default
    # surivival function).
    p_val = poisson.sf(x, mu=mu, loc=1)  # type: float

    return p_val


def count_region(
        reference_seq,  # type: pyfaidx.Fasta
        region,  # type: Tuple[str, int, int]
        pattern=None  # type: Optional[str]
):  # type: (...) -> int
    """Counts occurrences of pattern within given genomic region.

    Parameters
    ----------
    reference : pyfaidx.Fasta
        Reference to count occurrences for.
    region: Tuple[str, int, int]
        Genomic region to search in.
    pattern : str
        Nucleotide sequence to count occurences for. If None, the
        length of the region is used.

    Returns
    -------
    int
        Number of occurrences of the pattern within the given region.


    """

    chrom, start, end = region
    seq = reference_seq[chrom][int(start):int(end)]

    return _count_sequence(seq, regex=_build_regex(pattern))


def _build_regex(pattern):
    # type: (str) -> Pattern[str]
    if pattern is not None:
        return re.compile(pattern + '|' + pattern[::-1])
    return None


def _count_sequence(sequence, regex=None):
    # type: (pyfaidx.Sequence, Pattern[str]) -> int
    """Counts occurrences of pattern in sequence.

    Parameters
    ----------
    sequence : pyfaidx.Sequence
        Sequence to search.
    regex : Pattern[str]
        Pattern to count.

    Returns
    -------
        int: Number of occurrences of pattern.

    """

    if regex is None:
        count = len(sequence)
    else:
        count = sum((1 for _ in regex.finditer(str(sequence))))

    return count


def count_total(
        reference_seq,  # type: pyfaidx.Sequence
        pattern=None,  # type: str
        intervals=None  # type: Iterable[Tuple[str, int, int]]
):  # type: (...) -> int
    """Counts total occurrences of pattern in reference.

    Parameters
    ----------
    reference : pyfaidx.Fasta
        Reference to count occurrences for.
    pattern : str
        Nucleotide sequence to search for. If None, the length of
        sequences is counted instead of pattern of occurrences.
    intervals : List[tuple(str, int, int)]
        List of genomic intervals to which search should be restricted.
        If None, the entire reference is used.

    Returns
    -------
    int
        Number of occurrences of the pattern within the given reference,
        or within the given intervals (if applicable).


    """

    regex = _build_regex(pattern)

    if intervals is None:
        # Simply count for the entire sequence.
        count = sum(_count_sequence(reference_seq[seq], regex=regex)
                    for seq in reference_seq.keys()) # yapf: disable
    else:
        # Flatten intervals, and then only count for sequences
        # within the flattened intervals.
        merged_intervals = list(merge_genomic_intervals(intervals))

        seqs = [
            reference_seq[chrom][start:end]
            for chrom, start, end in merged_intervals
        ]

        count = sum(_count_sequence(seq, regex=regex) for seq in seqs)

    return count


def merge_genomic_intervals(intervals):
    # type: (Iterable[Tuple[str, int, int]]) -> Iterable[Tuple[str, int, int]]
    """Merges overlapping genomic intervals.

    Parameters
    ----------
    intervals : List[tuple(str, int, int)]
        List of intervals to merge. Intervals are specified as tuples
        of (chomosome, start_position, end_position) values.

    Yields:
    tuple(str, int, int)
        Next merged interval.

    """

    # Group intervals by chromosome.
    grouped_intervals = itertools.groupby(
        sorted(intervals), operator.itemgetter(0))

    # Now yield merged intervals per chromosome.
    for chrom, grp in grouped_intervals:
        chrom_intervals = [interval[1:] for interval in grp]
        for low, high in merge_intervals(chrom_intervals, is_sorted=True):
            yield chrom, low, high


def merge_intervals(
        intervals,  # type: Iterable[Tuple[int, int]]
        is_sorted=False  # type: Optional[bool]
):  # type: (...) -> Iterable[Tuple[int, int]]
    """Merges overlapping intervals.

    Parameters
    ----------
    intervals : List[tuple(int, int)]
        List of intervals to merge. Intervals are specified as tuples
        of (start_position, end_position) values.
    is_sorted : bool
        Indicates if the intervals have already been sorted. Avoids
        re-sorting an already sorted list.

    Yields:
    tuple(int, int)
        Next merged interval.
    """

    if not is_sorted:
        intervals = sorted(intervals, key=operator.itemgetter(0))
    else:
        # Ensure intervals is a list.
        intervals = list(intervals)

    if not intervals:
        # No intervals to merge.
        return

    # Low and high represent the bounds of
    # the current run of merges.
    low, high = intervals[0]

    for iv in intervals[1:]:
        if iv[0] <= high:
            # New interval overlaps current run,
            # merge with the current run.
            high = max(high, iv[1])
        else:
            # Current run is over, yield accumulated
            # interval and start new run.
            yield low, high
            low, high = iv

    yield low, high
