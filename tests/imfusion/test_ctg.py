# -*- coding: utf-8 -*-
"""Tests for imfusion.ctg module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from collections import namedtuple

import pyfaidx
import pytest

from imfusion import ctg
from imfusion.build import Reference
from imfusion.model import Insertion
from imfusion.util.frozendict import frozendict

Sequence = namedtuple('Sequence', ['seq'])
Gene = namedtuple('Gene', ['contig', 'start', 'end', 'strand'])

# pylint: disable=redefined-outer-name,no-self-use,too-few-public-methods
# pylint: disable=protected-access,C0103


def _insertion(id,
               seqname,
               position,
               strand,
               support_junction=1,
               support_spanning=1,
               metadata=None):
    """Helper function for building an Insertion instance."""
    return Insertion(
        id=id,
        seqname=seqname,
        position=position,
        strand=strand,
        support_junction=support_junction,
        support_spanning=support_spanning,
        support=support_junction + support_spanning,
        metadata=frozendict(metadata or {}))


class TestMergeIntervals(object):
    """Tests for the merge_intervals function."""

    def test_merge_intervals(self):
        """Test using example intervals."""
        intervals = [(0, 10), (14, 20), (18, 30), (30, 32), (33, 35)]
        merged = list(ctg.merge_intervals(intervals))
        assert merged == [(0, 10), (14, 32), (33, 35)]

    def test_merge_intervals_empty(self):
        """Test empty interval case."""
        merged = list(ctg.merge_intervals([]))
        assert merged == []


class TestMergeGenomicIntervals(object):
    """Tests for the merge_genomic_intervals function."""

    def test_merge_genomic_intervals(self):
        """Test using example intervals."""

        intervals = [('1', 10, 20), ('1', 22, 30), ('1', 25, 32), ('2', 8, 15),
                     ('3', 25, 40), ('3', 35, 40)]
        merged = list(ctg.merge_genomic_intervals(intervals))

        assert merged == [('1', 10, 20), ('1', 22, 32), ('2', 8, 15),
                          ('3', 25, 40)]

    def test_merge_genomic_intervals_empty(self):
        """Test empty interval case."""
        merged = list(ctg.merge_genomic_intervals([]))
        assert merged == []


@pytest.fixture
def reference_path():
    """Returns path to example reference sequence."""
    return pytest.helpers.data_path('reference.fa', relative_to=__file__)


@pytest.fixture
def reference(reference_path):
    """Returns example reference sequence."""
    return pyfaidx.Fasta(str(reference_path))


class TestCountRegion(object):
    """Tests for count_region function."""

    def test_count_region(self, reference):
        """Test counting without pattern."""
        assert ctg.count_region(reference, region=('1', 10, 30)) == 20

    def test_count_region_pattern(self, reference):
        """Test counting with pattern."""
        assert ctg.count_region(
            reference, region=('1', 10, 30), pattern='AG|GA') == 3


class TestCountTotal(object):
    """Tests for count_total function."""

    def test_count_total(self, reference):
        """Test without intervals, without pattern."""
        assert ctg.count_total(reference) == 280

    def test_count_total_pattern(self, reference):
        """Test without intervals, with pattern."""
        assert ctg.count_total(reference, pattern='TA') == 34

    def test_count_total_intervals(self, reference):
        """Test with intervals, without pattern."""
        intervals = [('1', 10, 20), ('1', 15, 25), ('2', 0, 10)]
        assert ctg.count_total(reference, intervals=intervals) == 25

    def test_count_total_intervals_pattern(self, reference):
        """Test with intervals, with pattern."""
        intervals = [('1', 10, 20), ('1', 15, 25), ('2', 0, 10)]
        assert ctg.count_total(
            reference, pattern='TA', intervals=intervals) == 3


@pytest.fixture
def insertions():
    """Example insertion set."""

    return [
        _insertion(id='1', seqname='1', position=9, strand=1,
                   metadata=frozendict({'gene_id': 'gene_a', 'sample': 'S1'})),
        _insertion(id='2', seqname='1', position=15, strand=-1,
                   metadata=frozendict({'gene_id': 'gene_b', 'sample': 'S2'}))
    ] # yapf: disable


class TestTestRegion(object):
    """Tests for test_region function."""

    def test_test_region(self, insertions, reference):
        """Test enriched region."""
        p_val = ctg.test_region(insertions, reference, region=('1', 5, 20))
        assert p_val < 0.01

    def test_test_region_negative(self, insertions, reference):
        """Test non-enriched region."""
        p_val = ctg.test_region(insertions, reference, region=('1', 10, 30))
        assert p_val > 0.01


class TestApplyWindow(object):
    """Tests for apply_window function."""

    def test_apply_window(self):
        """Tests example on forward strand."""

        new_window = ctg._apply_gene_window(
            Gene('1', 100, 120, '+'), window=(80, 50))

        assert new_window == ('1', 20, 170)

    def test_apply_window_rev(self):
        """Tests example on reverse strand."""

        new_window = ctg._apply_gene_window(
            Gene('1', 100, 120, '-'), window=(80, 50))

        assert new_window == ('1', 50, 200)

    def test_apply_window_none(self):
        """Tests example without window."""

        new_window = ctg._apply_gene_window(
            Gene('1', 100, 120, '-'), window=None)

        assert new_window == ('1', 100, 120)

    def test_apply_window_wrong_strand(self):
        """Tests example without proper strand."""

        with pytest.raises(ValueError):
            ctg._apply_gene_window(Gene('1', 100, 120, None), window=(80, 50))


class TestSubsetToWindows(object):
    """Tests subset_to_windows function."""

    def test_subset_insertions(self, insertions):
        """Test example."""

        windows = {'gene_a': ('1', 8, 12), 'gene_b': ('2', 10, 20)}
        subset = ctg._subset_to_windows(insertions, windows)

        assert len(subset) == 1
        assert subset[0].seqname == '1'

    def test_subset_insertions_no_overlap(self, insertions):
        """Test example with no insertions within windows."""

        windows = {'gene_a': ('1', 100, 120), 'gene_b': ('2', 10, 20)}
        assert len(ctg._subset_to_windows(insertions, windows)) == 0

    def test_subset_insertions_no_seqname(self, insertions):
        """Test example with overlapping position on different sequence."""

        windows = {'gene_a': ('2', 100, 120), 'gene_b': ('2', 10, 20)}
        assert len(ctg._subset_to_windows(insertions, windows)) == 0

    def test_subset_insertions_wrong_gene(self, insertions):
        """Test example."""

        windows = {'gene_a': ('1', 8, 12), 'gene_c': ('1', 10, 20)}
        subset = ctg._subset_to_windows(insertions, windows)

        assert len(subset) == 1
        assert subset[0].seqname == '1'


class TestCollapsePerSample(object):
    """Tests for collapse_per_sample function."""

    def test_example(self, insertions):
        """Tests example with collapsing."""

        insertions[1] = insertions[1]._replace(
            metadata={'gene_id': 'gene_a',
                      'sample': 'S1'})
        merged = list(ctg._collapse_per_sample(insertions))

        assert len(merged) == 1
        assert merged[0].position == 12

    def test_negative_example(self, insertions):
        """Tests example without collapsing."""

        merged = list(ctg._collapse_per_sample(insertions))
        assert merged == insertions


@pytest.fixture
def ctg_insertions():
    """Insertions for test_ctg test case."""
    return [
        _insertion(id='1', seqname='1', position=9, strand=1,
                   metadata={'gene_id': 'gene_a', 'sample': 'S1'}),
        _insertion(id='2', seqname='1', position=9, strand=1,
                   metadata={'gene_id': 'gene_a', 'sample': 'S1'}),
        _insertion(id='3', seqname='1', position=8, strand=-1,
                   metadata={'gene_id': 'gene_a', 'sample': 'S1'}),
        _insertion(id='4', seqname='1', position=8, strand=-1,
                   metadata={'gene_id': 'gene_b', 'sample': 'S2'}),
        _insertion(id='5', seqname='2', position=12, strand=-1,
                   metadata={'gene_id': 'gene_c', 'sample': 'S2'}),
        _insertion(id='6', seqname='1', position=6, strand=-1,
                   metadata={'gene_id': 'gene_a', 'sample': 'S3'})
    ] # yapf: disable


@pytest.fixture
def ctg_reference():
    """Reference for test_ctg test case."""
    return Reference(
        pytest.helpers.data_path(
            'ctg_reference', relative_to=__file__))


class TestTestCtgs(object):
    """Tests for the test_ctgs_function."""

    def test_example(self, ctg_insertions, ctg_reference):
        """Test example with three genes."""

        # TODO: P-values seem higher than in previous tests. Check why.

        # Do CTG test.
        result = ctg.test_ctgs(ctg_insertions, ctg_reference, per_sample=False)
        result = result.set_index('gene_id')

        # Check results.
        assert len(result) == 3
        assert result.ix['gene_a', 'p_value'] < 0.05
        assert result.ix['gene_b', 'p_value'] > 0.05
        assert result.ix['gene_c', 'p_value'] > 0.05

    def test_example_with_collapse(self, ctg_insertions, ctg_reference):
        """Tests if gene_a is no longer significant after collapsing."""

        # Do CTG test.
        result = ctg.test_ctgs(ctg_insertions, ctg_reference, per_sample=True)
        result = result.set_index('gene_id')

        # Check results.
        assert len(result) == 3
        assert result.ix['gene_a', 'p_value'] > 0.05
        assert result.ix['gene_b', 'p_value'] > 0.05
        assert result.ix['gene_c', 'p_value'] > 0.05

    def test_example_with_chromosomes(self, ctg_insertions, ctg_reference):
        """Tests subsetting for specific chromosomes."""

        # Do CTG test.
        result = ctg.test_ctgs(
            ctg_insertions, ctg_reference, chromosomes=['1'], per_sample=True)

        assert len(result) == 2
        assert set(result['gene_id']) == {'gene_a', 'gene_b'}

    def test_example_with_window(self, ctg_insertions, ctg_reference):
        """Tests applying a gene window."""

        # TODO: check the generated windows.

        # Do CTG test.
        result = ctg.test_ctgs(
            ctg_insertions, ctg_reference, window=(4, 0), per_sample=False)
        result = result.set_index('gene_id')

        # Check result.
        assert result.ix['gene_a', 'p_value'] < 0.05
