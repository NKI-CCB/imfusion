# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

from future.utils import native_str, native

import itertools
from functools import reduce

try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt

from imfusion.util.tabix import GtfFile

from .util import read_exon_counts
from .norm import estimate_size_factors
from .stats import NegativeBinomial


# Exon tests #

def de_exon(insertions, gene_id, dexseq_gtf, exon_counts,
            pos_samples=None, neg_samples=None):
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
    insertions : pandas.DataFrame
        DataFrame containing all insertions.
    gene_id : str
        ID of the gene of interest. Should correspond with a
        gene in the DEXSeq gtf file.
    dexseq_gtf : GtfFile
        Gtf file containing exon features generated using DEXSeq. Can either
        be given as a GtfFile object or as a string specifying the path to
        the gtf file.
    exon_counts : pandas.DataFrame or pathlib.Path
        DataFrame containing exon counts. The DataFrame is expected
        to contain samples as columns, and have a multi-index containing
        the chromosome, start, end and strand of the exon. This index
        should correspond with the annotation in the DEXSeq gtf. The
        samples should correspond with samples in the insertions frame.
        If a Path is given, it should point to a TSV file containing the counts.
    pos_samples : set[str]
        Set of positive samples (with insertion) to use in the test. Defaults
        to all samples with an insertion in the gene of interest.
    neg_samples : set[str]
        Set of negative samples (without insertion) to use in the test.
        Defaults to all samples not in the positive set.

    Returns
    -------
    DeExonResult
        Result of the differential expression test.

    """

    # Lookup gene.
    # TODO: speed-up using chromosome from insertions?
    gene = dexseq_gtf.get_gene(gene_id, feature_type='aggregate_gene')

    if isinstance(exon_counts, (str, native_str, pathlib.Path)):
        exon_counts = read_exon_counts(exon_counts, gene_id=gene.gene_id)

    # Subset insertions for gene.
    ins_gene = insertions.ix[insertions['gene_id'] == gene.gene_id]

    # Calculate before/after expression for gene.
    sums, exon_split, dropped_samples = \
        _exon_sums(gene, ins_gene, dexseq_gtf, exon_counts)

    # Define pos/neg sample groups.
    if pos_samples is None:
        pos_samples = set(ins_gene['sample_id'])

    if neg_samples is None:
        neg_samples = set(exon_counts.columns) - pos_samples

    pos_samples -= dropped_samples
    neg_samples -= dropped_samples

    # Check for missing samples.
    missing_samples = (pos_samples | neg_samples) - set(exon_counts.columns)
    if len(missing_samples) > 0:
        raise ValueError('Missing samples in counts ({})'
                         .format(', '.join(missing_samples)))

    # Calculate p-value between groups.
    pos_sums = sums.ix[pos_samples]['after']
    neg_sums = sums.ix[neg_samples]['after']

    p_value = mannwhitneyu(pos_sums, neg_sums)[1]

    direction = 1 if pos_sums.mean() > neg_sums.mean() else -1

    # Return result.
    return DeExonResult(sums, (pos_samples, neg_samples),
                        exon_split, direction, p_value)


def _exon_sums(gene, insertions, exons, exon_counts):
    # Split exons into exons before/after insertions.
    before, after, dropped_samples = \
        split_exons_common(exons, gene, insertions)
    before, after = map(_to_count_index, (before, after))

    # Convert to native types.
    before = native(before)
    after = native(after)

    # Normalize counts using before exons.
    counts = exon_counts.ix[gene.gene_id]

    size_factors = estimate_size_factors(counts.ix[before] + 1)
    norm_counts = counts / size_factors

    # Calculate before/after sums.
    sums = pd.DataFrame({
            'before': norm_counts.ix[before].sum(),
            'after': norm_counts.ix[after].sum()
        }, columns=['before', 'after']).round().astype(np.int)

    return sums, (before, after), dropped_samples


def de_exon_single(insertions, gene_id, insertion_id, dexseq_gtf, exon_counts):
    """Performs the single-sample exon-level differential expression test.

    Tests if the expression of exons after the insertion site of the
    given sample is significantly increased or decreased compared to samples
    without an insertion. This test is performed by comparing the (normalized)
    after count of the given sample to a background distribution of normalized
    counts of samples without an insertion, which is modeled using a
    negative binomial distribution.

    Note: this function requires Rpy2 to be installed, as R functions
    are used to fit the negative binomial distribution.

    Parameters
    ----------
    insertions : pandas.DataFrame
        DataFrame containing all insertions.
    gene_id : str
        ID of the gene of interest. Should correspond with a
        gene in the DEXSeq gtf file.
    insertion_id : str
        ID of the insertion of interest. Should correspond with
        an insertion in the list of insertions.
    dexseq_gtf : GtfFile
        Gtf file containing exon features generated using DEXSeq. Can either
        be given as a GtfFile object or as a string specifying the path to
        the gtf file.
    exon_counts : pandas.DataFrame
        DataFrame containing exon counts. The DataFrame is expected
        to contain samples as columns, and have a multi-index containing
        the chromosome, start, end and strand of the exon. This index
        should correspond with the annotation in the DEXSeq gtf. The
        samples should correspond with samples in the insertions frame.

    Returns
    -------
    DeExonResult
        Result of the differential expression test.

    """

    raise NotImplementedError()

    # Calculate before/after expression for gene.
    sums, exon_split = _exon_sums(gene, pd.DataFrame([insertion]),
                                  exons, exon_counts)

    # Split samples into positive/negative samples.
    pos_sample = insertion['sample_id']

    insertions_gene = insertions.ix[insertions['gene_id'] == gene.gene_id]
    neg_samples = set(exon_counts.columns) - set(insertions_gene['sample_id'])

    # Fit model and calculate p-value.
    nb_fit = NegativeBinomial(data=sums['after'][neg_samples])

    cdf = nb_fit.cdf(sums['after'][pos_sample])
    sf = nb_fit.sf(sums['after'][pos_sample])

    if cdf > sf:
        p_value, direction = sf, 1
    else:
        p_value, direction = cdf, -1

    # Return result.
    return DeExonResult(sums, ([pos_sample], neg_samples),
                        exon_split, nb_fit, direction, p_value)


# Result classes #

class DeExonResult(object):
    """Class embodying the results of the groupwise exon-level DE test.

    Attributes
    ----------
    sums : pandas.DataFrame
        DataFrame of before/after expression counts for all samples.
    sample_split : tuple(List[str], List[str])
        Split of samples into positive/negative samples.
    exon_split : tuple
        Split of exons into before/after groups.
    direction : int
        Direction of the differential expression (1 = positive, -1 = negative).
    p_value : float
        P-value of the differential expression test.

    """

    def __init__(self, sums, sample_split, exon_split,
                 direction, p_value):
        self.sums = sums
        self.sample_split = sample_split
        self.exon_split = exon_split
        self.direction = direction
        self.p_value = p_value

    def plot_boxplot(self, log=False, ax=None, show_points=True, **kwargs):
        """Plots boxplot of 'after' expression for samples with/without
        insertions in the gene."""

        ax = ax or plt.subplots()[1]

        # Log transform data if needed.
        if log:
            data = np.log2(self.sums + 1)
        else:
            data = self.sums.copy()

        # Plot using Seaborn.
        pos_samples = set(self.sample_split[0])
        data['insertion'] = [s in pos_samples for s in data.index]

        sns.boxplot(data=data, x='insertion', y='after',
                    showfliers=not show_points, order=[False, True],
                    ax=ax, **kwargs)

        if show_points:
            sns.stripplot(data=data, x='insertion', y='after',
                          order=[False, True], ax=ax, color='black',
                          size=3.5, jitter=True)

        ax.set_xlabel('Has insertion')
        ax.set_ylabel('Normalized expression after')

        return ax

    def plot_sums(self, log=False, **kwargs):
        """Plots the distribution of before/after counts for the samples."""
        sums = np.log2(self.sums + 1) if log else self.sums
        return plot_sums(sums, self.sample_split[0],
                         self.sample_split[1], **kwargs)


class DeExonSingleResult(object):
    """Class containing the results of the single-sample exon-level DE test.

    Attributes
    ----------
    sums : pandas.DataFrame
        DataFrame of before/after expression counts for all samples.
    sample_split : tuple(List[str], List[str])
        Split of samples into positive/negative samples.
    exon_split : tuple
        Split of exons into before/after groups.
    nb_fit : imfusion.expression.de_test.stats.NegativeBinomial
        Fit negative-binomial background distribution.
    direction : int
        Direction of the differential expression (1 = positive, -1 = negative).
    p_value : float
        P-value of the differential expression test.

    """

    def __init__(self, sums, sample_split, exon_split,
                 nb_fit, direction, p_value):
        self.sums = sums
        self.sample_split = sample_split
        self.exon_split = exon_split
        self.nb_fit = nb_fit
        self.direction = direction
        self.p_value = p_value

    def plot_fit(self, ax=None):
        """Plots the sample expression on the background distribution."""
        ax = ax or plt.subplots()[1]

        ax.hist(self.sums['after'][self.sample_split[1]],
                bins=25, normed=True, alpha=0.6)

        ax.plot([self.sums['after'][self.sample_split[0]]] * 2,
                ax.get_ylim(), linestyle='dashed', color='red')

        x = np.arange(*ax.get_xlim())
        ax.plot(x, self.nb_fit.pdf(x), 'k', linewidth=2)

        return ax

    def plot_sums(self, log=False, **kwargs):
        """Plots the distribution of before/after counts for the samples."""
        sums = np.log2(self.sums + 1) if log else self.sums
        return plot_sums(sums, [self.sample_split[0]],
                         self.sample_split[1], **kwargs)


# Splitting code #

def split_exons(exons, gene, insertion):
    assert insertion['gene_id'] == gene['gene_id']

    # Retrieve exons for gene.
    gene_exons = exons.get_region(
        gene.contig, gene.start, gene.end,
        filters={'feature': 'exonic_part',
                 'gene_id': gene.gene_id})

    # Split into exons before/after.
    before_exons = gene_exons.get_region(
        gene.contig, gene.start, insertion.position,
        incl_left=True, incl_right=False)

    after_exons = gene_exons.get_region(
        gene.contig, insertion.position, gene.end,
        incl_left=False, incl_right=True)

    # Mirror for genes on reverse strand.
    if gene.strand == '-':
        before_exons, after_exons = after_exons, before_exons

    return before_exons, after_exons


def split_exons_common(exons, gene, insertions):
    # Select insertions for gene.
    insertions = insertions.ix[insertions['gene_id'] == gene.gene_id]

    # Generate each split, tagged with sample name.
    sample_splits = [(ins['sample_id'], split_exons(exons, gene, ins))
                     for _, ins in insertions.iterrows()]

    # Check for samples with improper splits.
    improper_samples = set()
    for sample, split in sample_splits:
        if len(split[0]) == 0 or len(split[1]) == 0:
            improper_samples |= {sample}

    # Drop any splits belonging to improper samples.
    splits = [split for sample, split in sample_splits
              if sample not in improper_samples]

    if len(splits) > 0:
        # Reduce remaining splits into common features.
        before, after = zip(*splits)
        common_before = reduce(lambda x, y: pd.merge(x, y), before)
        common_after = reduce(lambda x, y: pd.merge(x, y), after)
    else:
        # Raise error if we have no samples to split by.
        raise ValueError('No proper samples to split by')

    return common_before, common_after, improper_samples


# Count lookup #

def _to_count_index(exons):
    return list(zip(exons['contig'], exons['start'],
                    exons['end'], exons['strand']))


# Sum plot #

def plot_sums(sums, pos_samples, neg_samples,
              width=0.5, **kwargs):
    fig, ax = plt.subplots(**kwargs)

    # Plot lines for positive and negative samples.
    # We keep track of plotted lines for mpld3.
    lines = dict()
    for sample in neg_samples:
        lines[sample] = _plot_sums_sample(
            ax, sums.ix[sample], width, color='lightgrey')

    palette = sns.color_palette()
    for sample, color in zip(pos_samples, palette):
        lines[sample] = _plot_sums_sample(
            ax, sums.ix[sample], width, color=color)

    # Set ticks.
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Before', 'After'])

    # Set axis limits/labels.
    ax.set_xlim(-0.5, 1.5)

    ax.set_xlabel('Location (relative to integration)')
    ax.set_ylabel('Expression')

    return fig, ax


def _plot_sums_sample(ax, sample, width, **kwargs):
    width /= 2

    l1 = ax.plot([-width, width], [sample.before] * 2, **kwargs)
    l2 = ax.plot([-width + 1, width + 1], [sample.after] * 2, **kwargs)

    l3 = ax.plot([width, -width+1], [sample.before, sample.after],
                 linestyle='dashed', **kwargs)

    return list(itertools.chain(l1, l2, l3))
