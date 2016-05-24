# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import numpy as np
import pandas as pd
import pytest

from imfusion.expression.de_test import test as de_exon
from imfusion.util.tabix import GtfFile

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


def _build_insertions(records):
    return pd.DataFrame.from_records(
        records, columns=['seqname', 'position', 'strand',
                          'gene_id', 'sample_id'])


@pytest.fixture
def dexseq_gtf_path():
    return pytest.helpers.data_path('dexseq.gtf.gz', relative_to=__file__)


@pytest.fixture()
def dexseq_gtf(dexseq_gtf_path):
    return GtfFile(dexseq_gtf_path)


@pytest.fixture()
def gene(dexseq_gtf):
    return dexseq_gtf.get_gene('ENSMUSG00000026510',
                               feature_type='aggregate_gene')


@pytest.fixture()
def insertions():
    return _build_insertions([
        ('1', 182429173, 1, 'ENSMUSG00000026510', 'S1'),
        ('1', 182457042, -1, 'ENSMUSG00000026510', 'S2')])


@pytest.fixture()
def gene_rev(dexseq_gtf):
    return dexseq_gtf.get_gene('ENSMUSG00000030849',
                               feature_type='aggregate_gene')


@pytest.fixture()
def insertions_rev():
    return _build_insertions([
        ('7', 130204552, 1, 'ENSMUSG00000030849', 'S1'),
        ('7', 130168882, -1, 'ENSMUSG00000030849', 'S2')])


@pytest.fixture()
def insertions_both(insertions, insertions_rev):
    return pd.concat([insertions, insertions_rev], ignore_index=True)


class TestSplitExonsCommon(object):

    def test_proper(self, dexseq_gtf, gene, insertions):
        before, after, dropped_samples = \
            de_exon.split_exons_common(dexseq_gtf, gene, insertions)

        assert all(before['end'] < insertions['position'].min())
        assert all(after['start'] > insertions['position'].max())
        assert dropped_samples == set()

    def test_proper_rev(self, dexseq_gtf, gene_rev, insertions_rev):
        before, after, dropped_samples = \
            de_exon.split_exons_common(dexseq_gtf, gene_rev, insertions_rev)

        assert all(before['start'] > insertions_rev['position'].max())
        assert all(after['end'] < insertions_rev['position'].min())
        assert dropped_samples == set()

    def test_proper_comb(self, dexseq_gtf, gene, insertions, insertions_rev):
        insertions_comb = pd.concat([insertions, insertions_rev],
                                    ignore_index=True)

        before, after, dropped_samples = \
            de_exon.split_exons_common(dexseq_gtf, gene, insertions_comb)

        assert all(before['end'] < insertions['position'].min())
        assert all(after['start'] > insertions['position'].max())
        assert dropped_samples == set()

    def test_before(self, dexseq_gtf, gene, insertions):
        """Checks if sample with insertion before gene is dropped."""

        # Add 'before' insertion to insertion frame.
        ins_before = _build_insertions([
            ('1', gene['start'] - 1000, -1, gene.gene_id, 'S3')])
        ins_with_before = pd.concat([insertions, ins_before],
                                    ignore_index=True)

        before, after, dropped_samples = \
            de_exon.split_exons_common(dexseq_gtf, gene, ins_with_before)

        # Check if split is unchanged.
        assert all(before['end'] < insertions['position'].min())
        assert all(after['start'] > insertions['position'].max())

        # Check if sample S3 is correctly dropped.
        assert dropped_samples == {'S3'}

    def test_after(self, dexseq_gtf, gene, insertions):
        """Checks if sample with insertion after gene is dropped."""

        # Add 'after' insertion to insertion frame.
        ins_after = _build_insertions([
            ('1', gene['end'] + 1000, -1, gene.gene_id, 'S3')])
        ins_with_after = pd.concat([insertions, ins_after],
                                   ignore_index=True)

        before, after, dropped_samples = \
            de_exon.split_exons_common(dexseq_gtf, gene, ins_with_after)

        # Check if split is unchanged.
        assert all(before['end'] < insertions['position'].min())
        assert all(after['start'] > insertions['position'].max())

        # Check if sample S3 is correctly dropped.
        assert dropped_samples == {'S3'}

    def test_no_proper(self, dexseq_gtf, gene):
        """Tests if error is raised if no proper samples are given."""

        ins_before = _build_insertions([
            ('1', gene['start'] - 1000, -1, gene.gene_id, 'S3')])

        with pytest.raises(ValueError):
            de_exon.split_exons_common(dexseq_gtf, gene, ins_before)

    def test_multiple_per_sample(self, dexseq_gtf, gene, insertions):
        """Test behaviour with multiple insertions in same sample."""

        extra_ins = _build_insertions([
            ('1', 182431766, 1, 'ENSMUSG00000026510', 'S1')])

        ins_with_extra = pd.concat([insertions, extra_ins],
                                   ignore_index=True)

        before, after, _ = de_exon.split_exons_common(
            dexseq_gtf, gene, ins_with_extra)

        # Check if split is unchanged.
        assert all(before['end'] < insertions['position'].min())
        assert all(after['start'] > insertions['position'].max())


@pytest.fixture()
def exon_counts(dexseq_gtf, gene, insertions):
    # Set seed.
    np.random.seed(1)

    # Create index from exons.
    exons = dexseq_gtf.get_region(
        filters={'gene_id': gene.gene_id,
                 'feature': 'exonic_part'})

    index = exons[['gene_id', 'contig', 'start', 'end', 'strand']]
    index = index.rename(columns={'contig': 'chr'})
    index['chr'] = index['chr'].astype(str)

    # Simulate counts.
    counts = {}
    for sample, grp in insertions.groupby('sample_id'):
        # Find split for sample.
        last_ins = grp['position'].max()
        mask = index['end'] < last_ins

        # Generate count vector.
        vec = np.empty(len(index))
        vec[:sum(mask)] = np.random.randint(-10, 10, sum(mask)) + 20
        vec[sum(mask):] = np.random.randint(-10, 10, sum(~mask)) + 100

        counts[sample] = vec

    # Generate some extra samples.
    for i in range(len(counts), 10):
        counts['S{}'.format(i + 1)] = np.random.randint(-10, 10, len(mask)) + 20

    # Create frame.
    counts = pd.concat([index, pd.DataFrame(counts)], axis=1)
    counts = counts.set_index(list(index.columns))
    counts = counts.astype(int)

    return counts


class TestDeExon(object):

    def test_simple_increase(self, gene, insertions, dexseq_gtf, exon_counts):
        res = de_exon.de_exon(insertions, gene.gene_id, dexseq_gtf, exon_counts)

        assert res.sample_split[0] == set(insertions['sample_id'])

        assert len(res.sample_split[1]) > 0
        assert all(s not in res.sample_split[1]
                   for s in set(insertions['sample_id']))

        assert res.p_value < 0.05
        assert res.direction == 1

    def test_missing_sample(self, gene, insertions, dexseq_gtf, exon_counts):
        exon_counts = exon_counts.drop('S2', axis=1)

        with pytest.raises(ValueError):
            de_exon.de_exon(insertions, gene.gene_id, dexseq_gtf, exon_counts)

    def test_plot_boxplot(self, gene, insertions, dexseq_gtf, exon_counts):
        res = de_exon.de_exon(insertions, gene.gene_id, dexseq_gtf, exon_counts)
        res.plot_boxplot()

    def test_plot_sums(self, gene, insertions, dexseq_gtf, exon_counts):
        res = de_exon.de_exon(insertions, gene.gene_id, dexseq_gtf, exon_counts)
        res.plot_sums()
