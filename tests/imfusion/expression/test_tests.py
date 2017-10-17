# -*- coding: utf-8 -*-
"""Tests for imfusion.expression.test module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import pytest

from imfusion.expression.counts import (ExonExpressionMatrix,
                                        GeneExpressionMatrix)
from imfusion.expression import tests
from imfusion.insertions import Insertion, InsertionSet

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def insertions():
    """Example insertions in gene_a (5 exons)."""
    return InsertionSet.from_tuples([
        # After first exon.
        Insertion(
            id='1',
            chromosome='1',
            position=25,
            strand=1,
            support_junction=1,
            support_spanning=1,
            support=2,
            gene_id='gene_a',
            sample='S1'),

        # After third exon.
        Insertion(
            id='2',
            chromosome='1',
            position=65,
            strand=-1,
            support_junction=1,
            support_spanning=1,
            support=2,
            gene_id='gene_a',
            sample='S2')
    ]) # yapf: disable


@pytest.fixture
def exon_counts():
    """Example exon counts."""

    file_path = pytest.helpers.data_path(
        'exon_counts.txt', relative_to=__file__)
    return ExonExpressionMatrix.from_imf(file_path)


class TestSplitCounts(object):
    """Tests for split_counts function."""

    # TODO: tests for reverse strand.

    def test_proper_example(self, insertions, exon_counts):
        """Tests example with two insertions, after exons 2 and 4."""

        before, after, dropped = tests.split_counts(
            exon_counts, insertions, gene_id='gene_a')

        assert len(before) == 1
        assert len(after) == 2
        assert dropped == set()

    def test_before_gene(self, insertions, exon_counts):
        """Tests dropping sample with an insertion before the gene"""

        insertions.values.loc[0, 'position'] = 5

        before, after, dropped = tests.split_counts(
            exon_counts, insertions, gene_id='gene_a')

        assert len(before) == 3
        assert len(after) == 2
        assert dropped == {'S1'}

    def test_after_gene(self, insertions, exon_counts):
        """Tests dropping sample with an insertion after the gene"""

        insertions.values.loc[1, 'position'] = 110

        before, after, dropped = tests.split_counts(
            exon_counts, insertions, gene_id='gene_a')

        assert len(before) == 1
        assert len(after) == 4
        assert dropped == {'S2'}

    def test_first_in_exon(self, insertions, exon_counts):
        """Tests first insertion in second exon.

        Should only pick first exon in before.
        """

        insertions.values.loc[0, 'position'] = 35

        before, after, dropped = tests.split_counts(
            exon_counts, insertions, gene_id='gene_a')

        assert len(before) == 1
        assert len(after) == 2
        assert dropped == set()

    def test_last_in_exon(self, insertions, exon_counts):
        """Tests second insertion in fourth exon.

        Should only pick last (fifth) exon as after.
        """

        insertions.values.loc[1, 'position'] = 75

        before, after, dropped = tests.split_counts(
            exon_counts, insertions, gene_id='gene_a')

        assert len(before) == 1
        assert len(after) == 1
        assert dropped == set()

    def test_no_valid_insertions(self, insertions, exon_counts):
        """Tests invalid splits."""

        # All insertions outside gene.
        insertions.values['position'] = [110, 120]

        with pytest.raises(ValueError):
            tests.split_counts(exon_counts, insertions, gene_id='gene_a')

        # All insertions in first/last exons.
        insertions.values['position'] = [15, 95]

        with pytest.raises(ValueError):
            tests.split_counts(exon_counts, insertions, gene_id='gene_a')

    def test_no_insertions(self, insertions, exon_counts):
        """Tests split with no insertions in gene."""

        with pytest.raises(ValueError):
            tests.split_counts(exon_counts, insertions, gene_id='gene_b')


@pytest.fixture
def sb_insertions():
    """Example insertions for Trp53bp2 and Nf1."""

    file_path = pytest.helpers.data_path(
        'insertions.txt', relative_to=__file__)
    return InsertionSet.from_csv(file_path, sep='\t')


@pytest.fixture
def sb_exon_counts():
    """Example exon counts for Trp53bp2 and Nf1."""

    file_path = pytest.helpers.data_path(
        'exon_counts.sb.txt', relative_to=__file__)
    #return read_exon_counts(file_path)
    return ExonExpressionMatrix.from_imf(file_path)


class TestTestDeExon(object):
    """Tests for test_de function."""

    def test_pos_example_trp53bp2(self, sb_insertions, sb_exon_counts):
        """Tests positive example of DE in Trp53bp2 in the SB dataset."""

        result = tests.test_de_exon(
            sb_insertions, sb_exon_counts, gene_id='ENSMUSG00000026510')

        assert result.p_value < 0.01
        assert result.direction == 1

    def test_neg_example_nf1(self, sb_insertions, sb_exon_counts):
        """Tests negative example of DE in Nf1 in the SB dataset."""

        result = tests.test_de_exon(
            sb_insertions, sb_exon_counts, gene_id='ENSMUSG00000020716')

        assert result.p_value > 0.05
        assert result.direction == 1

    def test_plot_boxplot(self, sb_insertions, sb_exon_counts):
        """Exercises plotting of boxplot."""

        result = tests.test_de_exon(
            sb_insertions, sb_exon_counts, gene_id='ENSMUSG00000026510')
        result.plot_boxplot()

    def test_plot_sums(self, sb_insertions, sb_exon_counts):
        """Exercises plotting of sums."""

        result = tests.test_de_exon(
            sb_insertions, sb_exon_counts, gene_id='ENSMUSG00000026510')
        result.plot_sums()


class TestDeExonSingle(object):
    """Tests for de_exon_single function."""

    def test_pos_example_trp53bp2(self, sb_insertions, sb_exon_counts):
        """Tests positive example in Trp53bp2."""

        result = tests.test_de_exon_single(
            sb_insertions,
            sb_exon_counts,
            insertion_id='2049_13_11KOU051-R3.INS_2',
            gene_id='ENSMUSG00000026510')

        assert result.p_value < 0.01
        assert result.direction == 1

    def test_neg_example_nf1(self, sb_insertions, sb_exon_counts):
        """Tests negative example in Nf1."""

        result = tests.test_de_exon_single(
            sb_insertions,
            sb_exon_counts,
            insertion_id='2800_69_13SKA003-R3.INS_8',
            gene_id='ENSMUSG00000020716')

        assert result.p_value > 0.05
        assert result.direction == 1


@pytest.fixture
def sb_gene_counts():
    """Example exon counts for Trp53bp2 and Nf1."""

    file_path = pytest.helpers.data_path(
        'gene_counts.sb.txt', relative_to=__file__)
    return GeneExpressionMatrix.from_subread(file_path)


class TestDeGene(object):
    """Tests for test_de_gene function."""

    def test_pos_example_trp53bp2(self, sb_insertions, sb_gene_counts):
        """Tests positive example in Trp53bp2."""

        result = tests.test_de_gene(
            sb_insertions, sb_gene_counts, gene_id='ENSMUSG00000026510')

        assert result.p_value < 0.01
        assert result.direction == 1

    def test_neg_example_nf1(self, sb_insertions, sb_gene_counts):
        """Tests negative example in Nf1."""

        result = tests.test_de_gene(
            sb_insertions, sb_gene_counts, gene_id='ENSMUSG00000020716')

        assert result.p_value > 0.05
        assert result.direction == -1

    def test_boxplot(self, sb_insertions, sb_gene_counts):
        """Tests positive example in Trp53bp2."""

        result = tests.test_de_gene(
            sb_insertions, sb_gene_counts, gene_id='ENSMUSG00000026510')
        result.plot_boxplot(log=True)
