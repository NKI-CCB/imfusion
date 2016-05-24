# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import itertools
import logging
import operator
import re

import toolz
import pandas as pd

from intervaltree import IntervalTree
from scipy.stats import poisson

from imfusion.expression.de_test import de_exon


def test_ctgs(insertion_frame, reference_seq, reference_gtf,
              chromosomes=None, pattern=None, gene_ids=None,
              per_sample=True, window=None, threshold=0.05):
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
    insertion_frame : pd.DataFrame
        Insertions to test (in DataFrame format).
    reference_sequence : pyfaidx.Fasta
        Fasta sequence of the reference genome.
    reference_gtf : GtfFile
        GtfFile containing reference genes.
    chromosomes : list[str]
        List of chromosomes to include, defaults to all chromosomes
        in reference_gtf.
    pattern : str
        Specificity pattern of the used transposon.
    genes : list[str]
        List of genes to test (defaults to all genes with an insertion).
    per_sample : bool
        Whether to perform the per sample test (recommended), which
        effectively collapes insertions per sample/gene combination.
        This avoids issues in which insertions that are detected
        multiple times or that may have hopped inside the gene locus
        are counted multiple times.
    window : tuple(int, int)
        Window to include around gene (in bp). Specified as (upstream_dist,
        downstream_dist). For example: (-2000, 2000) specifies in a 2KB
        window around each gene.
    threshold : float
        Maximum p-value for selected CTGs.

    Returns
    -------
    pandas.DataFrame
        Results of CTG test for tested genes. Contains three columns:
        gene_id, p_val and p_val_corr. The last column, p_val_corr,
        represents the p-value of the gene after correcting for
        multiple testing using bonferroni correction.

    """

    # Read reference genes from GTF.
    logging.info('-- Reading reference genes')

    ref_genes = reference_gtf.get_region(
        filters={'feature': 'gene'})
    ref_genes.set_index('gene_id', drop=False, inplace=True)

    # Subset reference genes to given chromosomes.
    if chromosomes is not None:

        ref_genes = ref_genes[ref_genes['contig'].isin(chromosomes)]

    # Determine gene intervals to include.
    intervals = (_apply_gene_window(gene, window)
                 for _, gene in ref_genes.iterrows())
    intervals = [int_ for int_ in intervals if int_[1] < int_[2]]

    # Subset insertions to gene intervals.
    insertion_frame = _subset_to_intervals(insertion_frame, intervals)

    # Collapse insertions per gene/sample (recommended).
    # Corrects for hopping/multiple detection issues.
    if per_sample:
        insertion_frame = _collapse_per_sample(insertion_frame)

    # Calculate total number of pattern occurrences within intervals.
    logging.info('-- Counting pattern occurrences')
    total = count_total(reference_seq, pattern=pattern, intervals=intervals)

    # Determine which genes to test, defaulting to all genes with an insertion.
    if gene_ids is None:
        gene_ids = set(insertion_frame['gene_id'].dropna())

    # Calculate p-values for each gene.
    logging.info('-- Calculating significance for genes')

    result = pd.DataFrame.from_records(
        ((gene_id, ref_genes.ix[gene_id]['gene_name'],
          test_region(*_apply_gene_window(ref_genes.ix[gene_id], window),
                     insertion_frame=insertion_frame,
                     reference=reference_seq,
                     total=total,
                     pattern=pattern,
                     filters={'gene_id': gene_id}))
         for gene_id in gene_ids),
        columns=['gene_id', 'gene_name', 'pval'])

    # Calculate corrected p-value using bonferroni correction.
    result['pval_corr'] = (result['pval'] * len(result)).clip_upper(1.0)

    # Threshold results and sort by p-value.
    result = result.ix[result['pval_corr'] <= threshold]
    result.sort_values(by='pval_corr', inplace=True)

    return result


def test_de(ctgs, insertions, dexseq_gtf, exon_counts_path, threshold=0.05):
    """Tests identified CTGs for differential expression.

    This function takes CTG frame produced by `test_ctgs` and tests each
    of the identified CTGs for differential expression using the groupwise
    exon-level differential expression test (`de_exon`). The resulting
    DE p-values are added to the DataFrame and CTGs that are not
    differentially expressed are dropped.

    Parameters
    ----------
    ctgs : pandas.DataFrame
        DataFrame containing the identified CTGs (as generated by `test_ctgs`).
    insertions : List[insertions]
        List of insertions to use in the test. Should be the same insertions
        as used to identify CTGs.
    dexseq_gtf : imfusion.util.tabix.GtfFile
        GtfFile instance containing the flattened exon representation of
        the original reference_gtf. The corresponding gtf file is typically
        generated using DEXSeqs script for preparing exon annotations.
    exon_counts_path : pathlib:Path
        Path to the file containing exon counts for all samples.
    threshold : float
        Maximum p-value for differential expression.

    Returns
    -------
    pandas.DataFrame
        CTG dataFrame containing the differential expression test results.

    """

    # Test each CTG genes for differential expression.
    de_results = {}
    for gene_id in ctgs['gene_id']:
        try:
            de_results[gene_id] = de_exon(insertions, gene_id,
                                          dexseq_gtf, exon_counts_path)
        except ValueError:
            pass

    # Summarize result in a DataFrame.
    de_results = pd.DataFrame(((gene_id, res.p_value, res.direction)
                               for gene_id, res in de_results.items()),
                              columns=['gene_id', 'de_pval', 'de_direction'])

    # Merge with CTG frame.
    merged = pd.merge(ctgs, de_results, on='gene_id')

    # Select significant diff. expr. CTGs.
    merged = merged.query('de_pval < {}'.format(threshold))

    return merged


def _subset_to_intervals(insertion_frame, intervals):
    """Subsets insertions for a list of genomic intervals."""

    # Create lookup trees.
    trees = {chrom: IntervalTree.from_tuples((i[1:]) for i in chrom_int)
             for chrom, chrom_int in
             itertools.groupby(sorted(intervals), operator.itemgetter(0))}

    # Determine which insertions overlap trees.
    def _in_intervals(ins):
        if ins.seqname not in trees:
            return False
        else:
            return trees[ins.seqname].overlaps(ins.position)

    return insertion_frame.ix[insertion_frame.apply(_in_intervals, axis=1)]


def _collapse_per_sample(insertion_frame):
    return (insertion_frame
            .groupby(['sample_id', 'gene_id'])
            .agg({'seqname': 'first',
                  'position': 'mean'})
            .reset_index())


def _apply_gene_window(gene, window=None):
    if window is None:
        return gene.contig, gene.start, gene.end,
    else:
        upstream_offset, downstream_offset = window

        if gene.strand == '-':
            start = gene.start - downstream_offset
            end = gene.end - upstream_offset
        elif gene.strand == '+':
            start = gene.start + upstream_offset
            end = gene.end + downstream_offset
        else:
            raise ValueError('Unknown value for strand')

        return gene.contig, start, end


def test_region(chrom, start, end, insertion_frame, reference,
                pattern=None, intervals=None, total=None, filters=None):
    """Tests a given genomic region for enrichment in insertions."""

    if total is None:
        total = count_total(reference, pattern=pattern, intervals=intervals)

    # Count pattern in region.
    region_count = count_region(chrom, start, end, reference, pattern)

    # Sub-select insertions for region.
    region_ins = insertion_frame.query(
        'seqname == {!r} and position > {} and position < {}'
        .format(chrom, start, end))

    # Apply additional filters to insertions if given
    # (such as filtering on gene name/id for example).
    if filters is not None:
        filter_str = ' and '.join('{} == {!r}'.format(field, value)
                                  for field, value in filters.items())
        region_ins = region_ins.query(filter_str)

    # Calculate p-value.
    x = len(region_ins)
    mu = len(insertion_frame) * (region_count / total)

    # Note here we use loc=1, because we are interested in
    # calculating P(X >= x), not P(X > x) (the default
    # surivival function).
    p_val = poisson.sf(x, mu=mu, loc=1)

    return p_val


def count_region(chrom, start, end, reference, pattern=None):
    """Counts occurrences of pattern within given genomic region.

    Parameters
    ----------
    chrom : str
        Region chromosome name.
    start : int
        Region start position.
    end : int
        Region end position.
    reference : pyfaidx.Fasta
        Reference to count occurrences for.
    pattern : str
        Regex pattern of string to search for. If None, the length
        of the region is returned.

    Returns
    -------
    int
        Number of occurrences of the pattern within the given region.


    """

    seq = reference[chrom][start:end]

    if pattern is not None:
        count = count_pattern(seq, pattern)
    else:
        count = len(seq)

    return count


def count_pattern(sequence, pattern):
    """Counts occurrences of pattern in sequence.

    Parameters
    ----------
    sequence : pyfaidx.Sequence or str
        Sequence to search.
    pattern : str
        Pattern to count.

    Return:
        int: Number of occurrences of pattern.

    """
    regex = re.compile(pattern)
    return count_matches(sequence, regex)


def count_matches(sequence, regex):
    """Counts matches of regex in string.

    Parameters:
    sequence : pyfaidx.Sequence or str
        Sequence to search.
    regex : SRE_Pattern
        Regex to apply.

    Returns
    -------
    int
        Number of occurrences of regex.

    """
    return sum((1 for _ in regex.finditer(str(sequence))))


def count_total(reference, pattern=None, intervals=None):
    """Counts total occurrences of pattern in reference.

    Parameters
    ----------
    reference : pyfaidx.Fasta
        Reference to count occurrences for.
    pattern : str
        Regex pattern of string to search for. If None, the length of
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

    if pattern is not None:
        count_func = toolz.curry(count_matches, regex=re.compile(pattern))
    else:
        # If no pattern is specified, use sequence length.
        count_func = len

    if intervals is None:
        # Simply count for the entire sequence.
        count = sum(count_func(reference[seq])
                    for seq in reference.keys())
    else:
        # Flatten intervals, and then only count for sequences
        # within the flattened intervals.
        merged_intervals = list(merge_genomic_intervals(intervals))

        seqs = [reference[chrom][start:end]
                for chrom, start, end in merged_intervals]

        count = sum(count_func(seq) for seq in seqs)

    return count


def merge_genomic_intervals(intervals):
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


def merge_intervals(intervals, is_sorted=False):
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
