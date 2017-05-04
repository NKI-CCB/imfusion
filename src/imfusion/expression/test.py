# -*- coding: utf-8 -*-
"""Contains functions that test for differential expression."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import bisect
import itertools
import logging
from typing import Iterable, Tuple, Union

from future.utils import native
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import toolz

from imfusion.model import Insertion
from .counts import estimate_size_factors, normalize_counts
from .stats import NegativeBinomial

MATPLOTLIB_IMPORT_ERR_MSG = (
    'Unable to import matplotlib/seaborn. Please make sure that both packages '
    'are installed and can be imported without any errors. (Often issues '
    'are due to misconfigured matplotlib backends.)')


def test_de(
        insertions,  # type: List[Insertion]
        exon_counts,  # type: pd.DataFrame
        gene_ids,  # type: List[str]
        fallback_to_gene=False,  # type: bool
        gene_counts=None  # type: pd.DataFrame
):  # type: (...) -> pd.DataFrame
    """Performs DE test for multiple genes and summarizes in table.

    Parameters
    ----------
    insertions : List[Insertion]
        Insertions to use for test.
    counts : pd.DataFrame
        Exon expression counts to use for test. Expected to conform to
        the format returned by the `read_exon_counts` function.
    gene_ids : List[str]
        IDs of genes to test for differential expression. Expected to conform
        with gene_ids used for the insertions and the expression counts.

    Returns
    -------
    pd.DataFrame
        DataFrame summarizing the results of the DE tests. Contains three
        columns: 'gene_id', 'p_value' and 'direction'. The p-value indicates
        the significance of the differential expression for the given gene,
        whilst the direction indicates whether the expression goes up (=1)
        in samples with an insertion or goes down (=-1).
    """

    rows = []

    for gene_id in gene_ids:
        try:
            # Determine p-value and direction.
            result = test_de_exon(insertions, exon_counts, gene_id=gene_id)
            p_value, dir_, type_ = result.p_value, result.direction, 'exon'
        except ValueError:
            # Failed to find split.
            if fallback_to_gene:
                logging.warning('Using gene test for %s', gene_id)

                # Fallback to gene test.
                if gene_counts is None:
                    logging.warning('Using summed exon counts to approximate '
                                    'gene counts for gene DE test')
                    gene_counts = exon_counts.groupby(level=0).sum()

                # Use helper function to avoid nesteed try/except.
                p_value, dir_, type_ = _test_gene(insertions, gene_counts,
                                                  gene_id)
            else:
                # Return NaNs for failed test.
                p_value, dir_, type_ = np.nan, np.nan, 'exon'

        rows.append((gene_id, p_value, dir_, type_))

    return pd.DataFrame.from_records(
        rows, columns=['gene_id', 'p_value', 'direction', 'test_type'])


def _test_gene(insertions, gene_counts, gene_id):
    try:
        result = test_de_gene(insertions, gene_counts, gene_id=gene_id)
        p_value, dir_, type_ = (result.p_value, result.direction, 'gene')
    except ValueError:
        p_value, dir_, type_ = np.nan, np.nan, 'exon'
    return p_value, dir_, type_


def test_de_exon(
        insertions,  # type: Union[List[Insertion], pd.DataFrame]
        exon_counts,  # type: pd.DataFrame
        gene_id,  # type: str
        pos_samples=None,  # type: Set[str]
        neg_samples=None  # type: Set[str]
):  # type: (...) -> DeResult
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
    insertions : List[Insertion] or pd.DataFrame
        List of insertions.
    exon_counts : pandas.DataFrame
        Matrix containing exon counts, with samples along the columns and
        exons along the rows. The DataFrame should have a multi-index
        containing the chromosome, start, end and strand of the exon. The
        samples should correspond with samples in the insertions frame.
    gene_id : str
        ID of the gene of interest. Should correspond with a
        gene in the count matrix.
    pos_samples : Set[str]
        Set of positive samples (with insertion) to use in the test. Defaults
        to all samples with an insertion in the gene of interest.
    neg_samples : Set[str]
        Set of negative samples (without insertion) to use in the test.
        Defaults to all samples not in the positive set.

    Returns
    -------
    DeResult
        Result of the differential expression test.

    """

    # Convert insertions to objects (if needed) and subset for gene.
    insertion_objs = _preprocess_insertions(insertions, gene_id)

    # Split counts by insertions.
    before, after, dropped_samples = split_counts(
        exon_counts, insertion_objs, gene_id=gene_id)

    # Define postive/negative sample groups (positive = with insertion).
    if pos_samples is None:
        pos_samples = set(ins.metadata['sample'] for ins in insertion_objs)

    if neg_samples is None:
        neg_samples = set(exon_counts.columns) - pos_samples

    pos_samples -= dropped_samples
    neg_samples -= dropped_samples

    if len(pos_samples) == 0:
        raise ValueError('No samples in positive set')

    if len(neg_samples) == 0:
        raise ValueError('No samples in negative set')

    # Normalize counts using before counts.
    size_factors = estimate_size_factors(before + 1)
    norm_before = before / size_factors
    norm_after = after / size_factors

    # Check for missing samples.
    missing_samples = (pos_samples | neg_samples) - set(exon_counts.columns)
    if len(missing_samples) > 0:
        raise ValueError('Missing samples in counts ({})'
                         .format(', '.join(missing_samples)))

    # Calculate p-value between groups.
    pos_sums = norm_after[list(pos_samples)].sum()
    neg_sums = norm_after[list(neg_samples)].sum()

    p_value = mannwhitneyu(pos_sums, neg_sums, alternative='two-sided')[1]

    direction = 1 if pos_sums.mean() > neg_sums.mean() else -1

    # Return result.
    return DeResult(norm_before, norm_after, pos_samples, neg_samples,
                    dropped_samples, direction, p_value)


def _preprocess_insertions(insertions, gene_id):
    """Converts insertions into common object format and subsets for gene."""

    if isinstance(insertions, pd.DataFrame):
        insertions = insertions.loc[insertions['gene_id'] == gene_id]
        insertions = list(Insertion.from_frame(insertions))
    else:
        insertions = [ins for ins in insertions
                      if ins.metadata['gene_id'] == gene_id]  # yapf: disable

    return insertions


def split_counts(
        counts,  # type: pd.DataFrame,
        insertions,  # type: Union[List[Insertion], pd.DataFrame]
        gene_id,  # type: str
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
    insertions : List[Insertion] or pandas.DataFrame
        List of identified insertions.
    gene_id : str
        Gene identifier.
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

    # Convert insertions to objects (if needed) and subset for gene.
    insertion_objs = _preprocess_insertions(insertions, gene_id)

    # Extract exon information from counts.
    exons = _get_exons(counts.loc[gene_id])
    strand = exons.iloc[0].strand

    # Switch limits if gene on - strand.
    if strand == -1:
        min_after, min_before = min_before, min_after

    max_after = len(exons) - min_after

    # Search for common split.
    dropped = set()  # type: Set[str]

    curr_min, curr_max = len(exons), 0
    for insertion in insertion_objs:
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

    # Apply split to full frame (includes gene_id).
    before = counts.loc[[gene_id]].iloc[:curr_min]
    after = counts.loc[[gene_id]].iloc[curr_max:]

    # Switch if gene on - strand.
    if strand == -1:
        before, after = after, before

    return before, after, dropped


def _get_exons(counts):
    # type: (pd.DataFrame) -> pd.DataFrame
    """Extracts exon position information from given count frame."""

    exons = pd.DataFrame.from_records(
        native(list(counts.index.get_values())),
        columns=['chromosome', 'start', 'end', 'strand'])
    exons['strand'] = exons['strand'].map({'-': -1, '+': 1})

    return exons


class DeResult(object):
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

        # Lazy load matplotlib/seaborn.
        try:
            from matplotlib import pyplot as plt
            import seaborn as sns
        except ImportError:
            raise ImportError(MATPLOTLIB_IMPORT_ERR_MSG)

        ax = ax or plt.subplots()[1]

        # Calculate sums.
        plot_data = pd.DataFrame({'after': self.after.sum(axis=0)})

        # Subset to samples and assign insertion status.
        samples = self.positive_samples | self.negative_samples
        plot_data = plot_data.loc[list(samples)]
        plot_data['insertion'] = [s in self.positive_samples
                                  for s in plot_data.index] # yapf: disable

        # Sanity check.
        stat = mannwhitneyu(
            plot_data.query('insertion == True').after,
            plot_data.query('insertion == False').after,
            alternative='two-sided')
        assert abs(stat[1] - self.p_value) < 1e-6

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
        ax.set_xticklabels(['Without\ninsertion', 'With\ninsertion'])

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

    # Lazy load matplotlib/seaborn.
    try:
        from matplotlib import pyplot as plt
        import seaborn as sns
    except ImportError:
        raise ImportError(MATPLOTLIB_IMPORT_ERR_MSG)

    if ax is None:
        _, ax = plt.subplots()

    line_kws = line_kws or {}

    plot_data = pd.DataFrame({
        'before': before.sum(axis=0),
        'after': after.sum(axis=0)
    })

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
    ax.set_xticklabels(['Before insertion', 'After insertion'])

    # Set axis limits/labels.
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel('Expression' if not log else 'Expression (log2)')

    return ax


def _plot_sums_sample(before, after, width, ax, **kwargs):
    width = width / 2
    ax.plot([-width, width], [before] * 2, **kwargs)
    ax.plot([-width + 1, width + 1], [after] * 2, **kwargs)
    ax.plot([width, -width + 1], [before, after], linestyle='dashed', **kwargs)


def test_de_exon_single(
        insertions,  # type: Union[List[Insertion], pd.DataFrame]
        exon_counts,  # type: pd.DataFrame
        insertion_id,  # type: str
        gene_id,  # type: str
        neg_samples=None  # type: Set[str]
):  # type: (...) -> DeSingleResult
    """Performs the single sample exon-level differential expression test.

    Tests if the expression of exons after the insertion site(s) in a
    gene is significantly increased or decreased in the sample with the given
    insertion compared to samples without an insertion (neg_samples).
    The test is performed by comparing normalized counts between the given
    sample and the set of background samples using the negative binomial
    distribution.

    Note that the before/after split for the groupwise test is taken as the
    common set of before/after exons over all samples with an insertion. In
    cases where either set is empty, for example due to insertions before the
    first exon of the gene, we attempt to drop samples that prevent a proper
    split and perform the test without these samples.

    Parameters
    ----------
    insertions : Union[List[Insertion], pd.DataFrame]
        List of insertions.
    exon_counts : pandas.DataFrame
        Matrix containing exon counts, with samples along the columns and
        exons along the rows. The DataFrame should have a multi-index
        containing the chromosome, start, end and strand of the exon. The
        samples should correspond with samples in the insertions frame.
    gene_id : str
        ID of the gene of interest. Should correspond with a
        gene in the count matrix.
    neg_samples : Set[str]
        Set of negative samples (without insertion) to use in the test.
        Defaults to all samples without an insertion in the gene.

    Returns
    -------
    DeSingleResult
        Result of the differential expression test.

    """

    # Convert insertions to objects (if needed) and subset for gene.
    insertion_objs = _preprocess_insertions(insertions, gene_id)

    # Extract selected insertion.
    insertions_by_id = {ins.id: ins for ins in insertion_objs}
    selected_ins = insertions_by_id[insertion_id]

    # Split counts by selected insertion.
    before, after, _ = split_counts(
        exon_counts, [selected_ins], gene_id=gene_id)

    # Define postive/negative sample groups (positive = with the insertion,
    # negative = all samples without an insertion).
    pos_sample = selected_ins.metadata['sample']

    if neg_samples is None:
        samples_with_ins = set(ins.metadata['sample']
                               for ins in insertion_objs)
        neg_samples = set(exon_counts.columns) - samples_with_ins

    if len(neg_samples) == 0:
        raise ValueError('No samples in negative set')

    # Normalize counts using before counts.
    size_factors = estimate_size_factors(before + 1)
    norm_before = before / size_factors
    norm_after = after / size_factors

    # Check for missing samples.
    missing_samples = ({pos_sample} | neg_samples) - set(exon_counts.columns)
    if len(missing_samples) > 0:
        raise ValueError('Missing samples in counts ({})'
                         .format(', '.join(missing_samples)))

    # Calculate p-value between groups.
    pos_sum = int(norm_after[pos_sample].sum())
    neg_sums = norm_after[list(neg_samples)].sum().astype(int)

    nb_distr = NegativeBinomial.fit(neg_sums)

    cdf_value = nb_distr.cdf(pos_sum)
    sf_value = nb_distr.sf(pos_sum)

    if cdf_value > sf_value:
        p_value, direction = sf_value, 1
    else:
        p_value, direction = cdf_value, -1

    # Return result.
    return DeSingleResult(norm_before, norm_after, pos_sample, neg_samples,
                          direction, p_value)


class DeSingleResult(object):
    """Class embodying the results of the single-sample exon-level DE test.

    Attributes
    ----------
    before : pandas.DataFrame
        Expression counts of exons before the split.
    after : pandas.DataFrame
        Expression counts of exons after the split.
    positive_sample : str
        Samples with an insertion.
    negative_samples : Set[str]
        Samples without an insertion.
    direction : int
        Direction of the differential expression (1 = positive, -1 = negative).
    p_value : float
        P-value of the differential expression test.

    """

    def __init__(self, before, after, pos_sample, neg_samples, direction,
                 p_value):
        self.before = before
        self.after = after
        self.positive_sample = pos_sample
        self.negative_samples = neg_samples
        self.direction = direction
        self.p_value = p_value

    def plot_sums(self, log=False, ax=None, line_kws=None):
        """Plots the distribution of before/after counts for the samples."""

        return _plot_sums(
            before=self.before,
            after=self.after,
            pos_samples=[self.positive_sample],
            neg_samples=self.negative_samples,
            log=log,
            ax=ax,
            line_kws=line_kws)


def test_de_gene(
        insertions,  # type: Union[List[Insertion], pd.DataFrame]
        gene_counts,  # type: pd.DataFrame
        gene_id,  # type: str
        pos_samples=None,  # type: Set[str]
        neg_samples=None  # type: Set[str]
):  # type: (...) -> DeGeneResult
    """Performs the groupwise gene-level differential expression test.

    Tests if the expression of the given gene is signficantly increased
    or decreased in samples with an insertion in the gene. Significance
    is calculated using a mannwhitneyu test between the two groups, after
    normalizing for sequencing depth using median-of-ratios normalization
    (as implemented in DESeq2).

    Parameters
    ----------
    insertions : Union[List[Insertion], pd.DataFrame]
        List of insertions.
    gene_counts : pandas.DataFrame
        Matrix containing gene counts, with samples along the columns and
        genes along the rows.
    gene_id : str
        ID of the gene of interest. Should correspond with a
        gene in the count matrix.
    pos_samples : Set[str]
        Set of positive samples (with insertion) to use in the test. Defaults
        to all samples with an insertion in the gene of interest.
    neg_samples : Set[str]
        Set of negative samples (without insertion) to use in the test.
        Defaults to all samples not in the positive set.

    Returns
    -------
    DeResult
        Result of the differential expression test.

    """

    # Normalize gene expression counts.
    norm_counts = normalize_counts(gene_counts)

    # Split into positive/negative samples.
    if pos_samples is None:
        if isinstance(insertions, pd.DataFrame):
            mask = insertions['gene_id'] == gene_id
            pos_samples = set(insertions.loc[mask]['sample'])
        else:
            pos_samples = set([
                ins.metadata['sample'] for ins in insertions
                if ins.metadata['gene_id'] == gene_id
            ])

    if neg_samples is None:
        neg_samples = set(gene_counts.columns) - pos_samples

    if len(pos_samples) == 0:
        raise ValueError('No samples in positive set')

    if len(neg_samples) == 0:
        raise ValueError('No samples in negative set')

    # Perform test.
    pos_counts = norm_counts.loc[gene_id, pos_samples]
    neg_counts = norm_counts.loc[gene_id, neg_samples]

    p_value = mannwhitneyu(pos_counts, neg_counts, alternative='two-sided')[1]
    direction = 1 if pos_counts.mean() > neg_counts.mean() else -1

    return DeGeneResult(
        counts=norm_counts.loc[gene_id],
        pos_samples=pos_samples,
        neg_samples=neg_samples,
        direction=direction,
        p_value=p_value)


class DeGeneResult(object):
    """Class embodying the results of the groupwise gene-level DE test.

    Attributes
    ----------
    counts : pandas.Series
        Expression counts of gene for the different samples.
    positive_samples : Set[str]
        Samples with an insertion.
    negative_samples : Set[str]
        Samples without an insertion.
    direction : int
        Direction of the differential expression (1 = positive, -1 = negative).
    p_value : float
        P-value of the differential expression test.

    """

    def __init__(self, counts, pos_samples, neg_samples, direction, p_value):
        self.counts = counts
        self.positive_samples = pos_samples
        self.negative_samples = neg_samples
        self.direction = direction
        self.p_value = p_value

    def plot_boxplot(self, ax=None, log=False, box_kws=None):
        """Plots boxplot comparing expression between the two groups."""

        # Lazy load matplotlib/seaborn.
        try:
            from matplotlib import pyplot as plt
            import seaborn as sns
        except ImportError:
            raise ImportError(MATPLOTLIB_IMPORT_ERR_MSG)

        if ax is None:
            _, ax = plt.subplots()

        # Assemble plot data.
        samples = list(self.positive_samples | self.negative_samples)

        plot_data = pd.DataFrame({'expression': self.counts.loc[samples]})
        plot_data['insertion'] = [sample in self.positive_samples
                                  for sample in plot_data.index] # yapf: disable

        # Log transform if needed.
        if log:
            plot_data['expression'] = np.log2(plot_data['expression'] + 1)

        # Draw boxplot with seaborn.
        sns.boxplot(
            data=plot_data,
            x='insertion',
            order=[False, True],
            y='expression',
            ax=ax,
            **(box_kws or {}))

        ax.set_xlabel('')
        ax.set_xticklabels(['Without\ninsertion', 'With\ninsertion'])

        ax.set_ylabel('Expression' + (' (log2)' if log else ''))
