"""Module containing functions for testing for differential expression."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import bisect
import itertools
from typing import Iterable, Tuple

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns
import toolz

from imfusion.model import Insertion
from .counts import estimate_size_factors


def de_exon(
        insertions,  # type: Iterable[Insertion]
        exon_counts,  # type: pd.DataFrame
        gene_id,  # type: str
        pos_samples=None,  # type: Set[str]
        neg_samples=None  # type: Set[str]
):  # type: (...) -> DeExonResult
    """Performs the groupwise exon-level differential expression test.

    Tests if the expression of exons after the insertion site(s) in a
    gene is significantly increased or decreased in samples with an
    insertion (pos_samples) compared to samples without an insertion
    (neg_samples). The test is performed by comparing normalized counts
    after the insertion sites between samples with and without an insertion
    in the gene, using the non-parametric Mann-Whitney-U test.

    Note that the before/after split for the groupwise test is taken as the
    common set of before/after exons over all samples with an insertion. In
    cases where either set is empty, for example due to insertions before the
    first exon of the gene, we attempt to drop samples that prevent a proper
    split and perform the test without these samples.

    Parameters
    ----------
    insertions : List[Insertion]
        List of insertions.
    exon_counts : pandas.DataFrame
        Matrix containing exon counts, with samples along the columns and
        exons along the rows. The DataFrame should have a multi-index
        containing the chromosome, start, end and strand of the exon. The
        samples should correspond with samples in the insertions frame.
    gene_id : str
        ID of the gene of interest. Should correspond with a
        gene in the count matrix.
    pos_samples : Optional[Set[str]]
        Set of positive samples (with insertion) to use in the test. Defaults
        to all samples with an insertion in the gene of interest.
    neg_samples : Optional[Set[str]]
        Set of negative samples (without insertion) to use in the test.
        Defaults to all samples not in the positive set.

    Returns
    -------
    DeExonResult
        Result of the differential expression test.

    """

    # Subset counts and exons for gene.
    counts = exon_counts.ix[gene_id]

    insertions = [ins for ins in insertions
                  if ins.metadata['gene_id'] == gene_id] # yapf: disable

    # Split counts by insertions.
    before, after, dropped_samples = split_counts(counts, insertions)

    # Define postive/negative sample groups (positive = with insertion).
    if pos_samples is None:
        pos_samples = set(ins.metadata['sample'] for ins in insertions)

    if neg_samples is None:
        neg_samples = set(counts.columns) - pos_samples

    pos_samples -= dropped_samples
    neg_samples -= dropped_samples

    # Normalize counts using before counts.
    size_factors = estimate_size_factors(before + 1)
    norm_after = after / size_factors

    # Check for missing samples.
    missing_samples = (pos_samples | neg_samples) - set(exon_counts.columns)
    if len(missing_samples) > 0:
        raise ValueError('Missing samples in counts ({})'
                         .format(', '.join(missing_samples)))

    # Calculate p-value between groups.
    pos_sums = norm_after[list(pos_samples)].sum()
    neg_sums = norm_after[list(neg_samples)].sum()

    p_value = mannwhitneyu(pos_sums, neg_sums)[1]

    direction = 1 if pos_sums.mean() > neg_sums.mean() else -1

    # Return result.
    return DeExonResult(before, after, pos_samples, neg_samples,
                        dropped_samples, direction, p_value)


def split_counts(
        counts,  # type: pd.DataFrame,
        insertions,  # type: List[Insertion]
        min_before=1,  # type: int
        min_after=1  # type: int
):  # type: (...) -> Tuple[pd.DataFrame, pd.DataFrame, Set[str]]
    """Splits count frame for exons before and after insertion sites.

    Parameters
    ----------
    counts : pd.DataFrame
        Matrix of exon expression counts, with positions along the rows
        and samples along the columns. The index of the DataFrame should
        be a multi-level index containing the following levels: gene_id,
        chromosome, start, end and strand.
    insertions : List[Insertion]
        List of identified insertions.
    min_before : int
        Minimum number of exons to retain before the split. Samples with less
        exons before their insertion sites will be dropped.
    min_after : int
        Minimum number of exons to retain after the split. Samples with less
        exons after their insertion sites will be dropped.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, Set[str]]

        Returns a tuple of two DataFrames, representing the counts before
        and after the split respectively, and a set of sample names that
        were dropped from the split for violating the min_before/min_after
        restrictions.

    """

    # TODO: tests for min_before/min_after.

    # Extract exon information from counts.
    exons = _get_exons(counts)
    strand = exons.iloc[0].strand

    # Switch limits if gene on - strand.
    if strand == -1:
        min_after, min_before = min_before, min_after

    max_after = len(exons) - min_after

    # Search for common split.
    dropped = set()  # type: Set[str]

    curr_min, curr_max = len(exons), 0
    for insertion in insertions:
        idx = bisect.bisect_right(exons.end, insertion.position)

        if idx < min_before or idx > max_after:
            # Add to invalid samples.
            dropped |= {insertion.metadata['sample']}
        else:
            if idx > curr_max:
                curr_max = idx

            if idx < curr_min:
                curr_min = idx

    if curr_min > curr_max:
        raise ValueError('No valid split found')

    # Apply split.
    before = counts.iloc[:curr_min]
    after = counts.iloc[curr_max:]

    # Switch if gene on - strand.
    if strand == -1:
        before, after = after, before

    return before, after, dropped


def _get_exons(counts):
    # type: (pd.DataFrame) -> pd.DataFrame
    """Extracts exon position information from given count frame."""

    exons = pd.DataFrame.from_records(
        list(counts.index.get_values()),
        columns=['chromosome', 'start', 'end', 'strand'])
    exons['strand'] = exons['strand'].map({'-': -1, '+': 1})

    return exons


class DeExonResult(object):
    """Class embodying the results of the groupwise exon-level DE test.

    Attributes
    ----------
    before : pandas.DataFrame
        Expression counts of exons before the split.
    after : pandas.DataFrame
        Expression counts of exons after the split.
    positive_samples : Set[str]
        Samples with an insertion.
    negative_samples : Set[str]
        Samples without an insertion.
    dropped_samples : Set[str]
        Samples that were dropped from the analysis, because of violating
        the constraints on the minimum number of exons before/after the split.
    direction : int
        Direction of the differential expression (1 = positive, -1 = negative).
    p_value : float
        P-value of the differential expression test.

    """

    def __init__(self, before, after, pos_samples, neg_samples,
                 dropped_samples, direction, p_value):
        self.before = before
        self.after = after
        self.positive_samples = pos_samples
        self.negative_samples = neg_samples
        self.dropped_samples = dropped_samples
        self.direction = direction
        self.p_value = p_value

    def plot_boxplot(self,
                     log=False,
                     ax=None,
                     show_points=True,
                     box_kws=None,
                     strip_kws=None):
        """Plots boxplot of 'after' expression for samples with/without
        insertions in the gene."""

        ax = ax or plt.subplots()[1]

        # Collect plot data.
        plot_data = pd.DataFrame({'after': self.after.sum()})
        plot_data['insertion'] = [
            s in self.positive_samples for s in plot_data.index
        ]

        # Log transform data if needed.
        if log:
            plot_data['after'] = np.log2(plot_data['after'] + 1)

        # Plot using Seaborn.
        sns.boxplot(
            data=plot_data,
            x='insertion',
            y='after',
            showfliers=not show_points,
            order=[False, True],
            ax=ax,
            **(box_kws or {}))

        if show_points:
            default_kws = dict(color='black', size=3.5, jitter=True)
            strip_kws = toolz.merge(default_kws, strip_kws or {})

            sns.stripplot(
                data=plot_data,
                x='insertion',
                y='after',
                order=[False, True],
                ax=ax,
                **strip_kws)

        ax.set_xlabel('')
        ax.set_xticklabels(['With\ninsertion', 'Without\ninsertion'])

        ylabel = 'Normalized expression after insertion sites'
        ax.set_ylabel(ylabel if not log else ylabel + ' (log2)')

        return ax

    def plot_sums(self, log=False, ax=None, line_kws=None):
        """Plots the distribution of before/after counts for the samples."""

        return _plot_sums(
            before=self.before,
            after=self.after,
            pos_samples=self.positive_samples,
            neg_samples=self.negative_samples,
            log=log,
            ax=ax,
            line_kws=line_kws)


def _plot_sums(before,
               after,
               pos_samples,
               neg_samples,
               ax=None,
               log=False,
               width=0.5,
               line_kws=None):
    """Helper function for plotting expression sums in line graph."""

    if ax is None:
        _, ax = plt.subplots()

    line_kws = line_kws or {}

    plot_data = pd.DataFrame({'before': before.sum(), 'after': after.sum()})

    if log:
        plot_data = np.log2(plot_data + 1)

    # Plot background lines for negative samples.
    for tup in plot_data.ix[neg_samples].itertuples():
        _plot_sums_sample(
            tup.before, tup.after, width, ax=ax, color='lightgrey', **line_kws)

    # Plot colored lines for positive samples.
    pos_samples = plot_data.ix[pos_samples].itertuples()
    colors = itertools.cycle(sns.color_palette())

    for tup, color in zip(pos_samples, colors):
        _plot_sums_sample(
            tup.before, tup.after, width, ax=ax, color=color, **line_kws)

    # Set ticks.
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['With\ninsertion', 'Without\ninsertion'])

    # Set axis limits/labels.
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel('Expression' if not log else 'Expression (log2)')

    return ax


def _plot_sums_sample(before, after, width, ax, **kwargs):
    width = width / 2
    ax.plot([-width, width], [before] * 2, **kwargs)
    ax.plot([-width + 1, width + 1], [after] * 2, **kwargs)
    ax.plot([width, -width + 1], [before, after], linestyle='dashed', **kwargs)


# def de_exon_single(insertions, gene_id, insertion_id, dexseq_gtf, exon_counts):
#     """Performs the single-sample exon-level differential expression test.

#     Tests if the expression of exons after the insertion site of the
#     given sample is significantly increased or decreased compared to samples
#     without an insertion. This test is performed by comparing the (normalized)
#     after count of the given sample to a background distribution of normalized
#     counts of samples without an insertion, which is modeled using a
#     negative binomial distribution.

#     Note: this function requires Rpy2 to be installed, as R functions
#     are used to fit the negative binomial distribution.

#     Parameters
#     ----------
#     insertions : pandas.DataFrame
#         DataFrame containing all insertions.
#     gene_id : str
#         ID of the gene of interest. Should correspond with a
#         gene in the DEXSeq gtf file.
#     insertion_id : str
#         ID of the insertion of interest. Should correspond with
#         an insertion in the list of insertions.
#     dexseq_gtf : GtfFile
#         Gtf file containing exon features generated using DEXSeq. Can either
#         be given as a GtfFile object or as a string specifying the path to
#         the gtf file.
#     exon_counts : pandas.DataFrame
#         DataFrame containing exon counts. The DataFrame is expected
#         to contain samples as columns, and have a multi-index containing
#         the chromosome, start, end and strand of the exon. This index
#         should correspond with the annotation in the DEXSeq gtf. The
#         samples should correspond with samples in the insertions frame.

#     Returns
#     -------
#     DeExonResult
#         Result of the differential expression test.

#     """

#     raise NotImplementedError()

#     # Calculate before/after expression for gene.
#     sums, exon_split = _exon_sums(gene,
#                                   pd.DataFrame([insertion]), exons,
#                                   exon_counts)

#     # Split samples into positive/negative samples.
#     pos_sample = insertion['sample']

#     insertions_gene = insertions.ix[insertions['gene_id'] == gene.gene_id]
#     neg_samples = set(exon_counts.columns) - set(insertions_gene['sample'])

#     # Fit model and calculate p-value.
#     nb_fit = NegativeBinomial(data=sums['after'][neg_samples])

#     cdf = nb_fit.cdf(sums['after'][pos_sample])
#     sf = nb_fit.sf(sums['after'][pos_sample])

#     if cdf > sf:
#         p_value, direction = sf, 1
#     else:
#         p_value, direction = cdf, -1

#     # Return result.
#     return DeExonResult(sums, ([pos_sample], neg_samples), exon_split, nb_fit,
#                         direction, p_value)

# class DeExonSingleResult(object):
#     """Class containing the results of the single-sample exon-level DE test.

#     Attributes
#     ----------
#     sums : pandas.DataFrame
#         DataFrame of before/after expression counts for all samples.
#     sample_split : tuple(List[str], List[str])
#         Split of samples into positive/negative samples.
#     exon_split : tuple
#         Split of exons into before/after groups.
#     nb_fit : imfusion.expression.de_test.stats.NegativeBinomial
#         Fit negative-binomial background distribution.
#     direction : int
#         Direction of the differential expression (1 = positive, -1 = negative).
#     p_value : float
#         P-value of the differential expression test.

#     """

#     def __init__(self, sums, sample_split, exon_split, nb_fit, direction,
#                  p_value):
#         self.sums = sums
#         self.sample_split = sample_split
#         self.exon_split = exon_split
#         self.nb_fit = nb_fit
#         self.direction = direction
#         self.p_value = p_value

#     def plot_fit(self, ax=None):
#         """Plots the sample expression on the background distribution."""
#         ax = ax or plt.subplots()[1]

#         ax.hist(
#             self.sums['after'][self.sample_split[1]],
#             bins=25,
#             normed=True,
#             alpha=0.6)

#         ax.plot(
#             [self.sums['after'][self.sample_split[0]]] * 2,
#             ax.get_ylim(),
#             linestyle='dashed',
#             color='red')

#         x = np.arange(*ax.get_xlim())
#         ax.plot(x, self.nb_fit.pdf(x), 'k', linewidth=2)

#         return ax

#     def plot_sums(self, log=False, **kwargs):
#         """Plots the distribution of before/after counts for the samples."""
#         sums = np.log2(self.sums + 1) if log else self.sums
#         return _plot_sums(sums, [self.sample_split[0]], self.sample_split[1],
#                           **kwargs)
