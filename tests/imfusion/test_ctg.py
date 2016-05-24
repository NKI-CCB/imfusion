# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import re
from collections import namedtuple

import pytest
import mock

import pandas as pd
import pyfaidx

from imfusion import ctg
from imfusion.model import Insertion
from imfusion.util.tabix import GtfFile


Sequence = namedtuple('Sequence', ['seq'])
Gene = namedtuple('Gene', ['contig', 'start', 'end', 'strand'])

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods

class TestCountMatches(object):

    def test_count_matches(self):
        seq = Sequence(seq='ATCGATCGTACG')
        assert ctg.count_matches(seq, regex=re.compile('AT')) == 2


    def test_count_matches_mirrored(self):
        seq = Sequence(seq='ATCGATCGTACG')
        assert ctg.count_matches(seq, regex=re.compile('AT|TA')) == 3


    def test_count_matches_missing(self):
        seq = Sequence(seq='ATCGATCGTACG')
        assert ctg.count_matches(seq, regex=re.compile('BG')) == 0


class TestCountPattern(object):

    def test_count_pattern(self):
        seq = Sequence(seq='ATCGATCGTACG')
        assert ctg.count_pattern(seq, pattern='AT') == 2


class TestMergeIntervals(object):

    def test_merge_intervals(self):
        intervals = [(0, 10), (14, 20), (18, 30), (30, 32), (33, 35)]
        merged = list(ctg.merge_intervals(intervals))
        assert merged == [(0, 10), (14, 32), (33, 35)]


    def test_merge_intervals_empty(self):
        merged = list(ctg.merge_intervals([]))
        assert merged == []


class TestMergeGenomicIntervals(object):

    def test_merge_genomic_intervals(self):
        intervals = [('1', 10, 20), ('1', 22, 30), ('1', 25, 32),
                     ('2', 8, 15), ('3', 25, 40), ('3', 35, 40)]
        merged = list(ctg.merge_genomic_intervals(intervals))

        assert merged == [('1', 10, 20), ('1', 22, 32),
                          ('2', 8, 15), ('3', 25, 40)]


    def test_merge_genomic_intervals_empty(self):
        merged = list(ctg.merge_genomic_intervals([]))
        assert merged == []


@pytest.fixture
def reference_path():
    return pytest.helpers.data_path('reference.fa')


@pytest.fixture
def reference(reference_path):
    return pyfaidx.Fasta(str(reference_path))


class TestCountRegion(object):

    def test_count_region(self, reference):
        assert ctg.count_region('1', 10, 30, reference) == 20


        def test_count_region_pattern(self, reference):
            assert ctg.count_region('1', 10, 30, reference,
                                    pattern='AG|GA') == 3


class TestCountTotal(object):

    def test_count_total(self, reference):
        assert ctg.count_total(reference) == 280


    def test_count_total_pattern(self, reference):
        assert ctg.count_total(reference, pattern='TA') == 21


    def test_count_total_intervals(self, reference):
        intervals = [('1', 10, 20), ('1', 15, 25), ('2', 0, 10)]
        assert ctg.count_total(reference, intervals=intervals) == 25


    def test_count_total_intervals_pattern(self, reference):
        intervals = [('1', 10, 20), ('1', 15, 25), ('2', 0, 10)]
        assert ctg.count_total(reference, pattern='TA',
                               intervals=intervals) == 2


class TestTestRegion(object):

    def test_test_region(self, insertions, reference):
        p_val = ctg.test_region('1', 5, 20, insertions, reference)
        assert p_val < 0.01


    def test_test_region_negative(self, insertions, reference):
        p_val = ctg.test_region('1', 10, 30, insertions, reference)
        assert p_val > 0.01


class TestApplyWindow(object):

    def test_apply_window(self):
        new_window = ctg._apply_gene_window(
            Gene('1', 100, 120, '+'), window=(-80, 50))
        assert new_window == ('1', 20, 170)


    def test_apply_window_rev(self):
        new_window = ctg._apply_gene_window(
            Gene('1', 100, 120, '-'), window=(-80, 50))
        assert new_window == ('1', 50, 200)


    def test_apply_window_none(self):
        new_window = ctg._apply_gene_window(
            Gene('1', 100, 120, '-'), window=None)
        assert new_window == ('1', 100, 120)


    def test_apply_window_wrong_strand(self):
        with pytest.raises(ValueError):
            ctg._apply_gene_window(Gene('1', 100, 120, None),
                                   window=(-80, 50))

@pytest.fixture
def insertions():
    insertions = [
        Insertion(id='1', seqname='1', position=9, strand=1,
                  gene_id='gene_a', sample_id='s1'),
        Insertion(id='2', seqname='1', position=15, strand=-1,
                  gene_id='gene_b', sample_id='s2')]
    return Insertion.to_frame(insertions)


class TestSubsetInsertions(object):

    def test_subset_insertions(self, insertions):
        intervals = [('1', 8, 12), ('2', 10, 20)]

        subset = ctg._subset_to_intervals(insertions, intervals)

        assert len(subset) == 1
        assert subset.ix[0, 'seqname'] == '1'


    def test_subset_insertions_no_overlap(self, insertions):
        intervals = [('1', 100, 120), ('2', 10, 20)]
        assert len(ctg._subset_to_intervals(insertions, intervals)) == 0


    def test_subset_insertions_no_seqname(self, insertions):
        intervals = [('2', 8, 12), ('2', 10, 20)]
        assert len(ctg._subset_to_intervals(insertions, intervals)) == 0

#
# @pytest.fixture()
# def ref_gtf():
#     return GtfFile(pytest.helpers.data_path('mm10.test.gtf.gz'))
#
#
# #    `id`,
# #    `seqname`, `position` and `strand`.
#
#
# @pytest.fixture
# def insertions_ctgs():
#     insertions = [
#         Insertion(id='1', seqname='1', position=9, strand=1,
#                   gene_id='gene_a', sample_id='s1'),
#         Insertion(id='2', seqname='1', position=9, strand=1,
#                   gene_id='gene_a', sample_id='s1'),
#         Insertion(id='3', seqname='1', position=8, strand=-1,
#                   gene_id='gene_a', sample_id='s1'),
#         Insertion(id='4', seqname='1', position=8, strand=-1,
#                   gene_id='gene_b', sample_id='s2'),
#         Insertion(id='5', seqname='2', position=12, strand=-1,
#                   gene_id='gene_c', sample_id='s2'),
#         Insertion(id='6', seqname='1', position=6, strand=-1,
#                   gene_id='gene_a', sample_id='s3')]
#     return Insertion.to_frame(insertions)
#
# @pytest.fixture
# def ref_genes():
#     return pd.DataFrame.from_records(
#         [('1', 8, 10, '+', 'gene_a'),
#          ('1', 30, 40, '-', 'gene_b'),
#          ('2', 10, 20, '+', 'gene_c')],
#         columns=['contig', 'start', 'end', 'strand', 'gene_id'])
#
#
# class TestIdentifyCtgs(object):
#
#     def test_example(self, insertions_ctgs, reference, ref_gtf, ref_genes):
#         with mock.patch('imfusion.util.tabix.GtfFrame.get_region',
#                         return_value=ref_genes):
#             result = ctg.test_ctgs(insertions_ctgs, reference,
#                                    ref_gtf, per_sample=False, threshold=1.0)
#             result.set_index('gene_id', inplace=True)
#
#             assert len(result) == 3
#             assert result.ix['gene_a', 'p_val_corr'] < 0.05
#             assert result.ix['gene_b', 'p_val_corr'] > 0.05
#             assert result.ix['gene_c', 'p_val_corr'] > 0.05
#
#
#     def test_per_sample(self, insertions_ctgs, reference, ref_gtf, ref_genes):
#         with mock.patch('imfusion.util.tabix.GtfFrame',
#                         return_value=ref_genes):
#             result = ctg.test_ctgs(insertions_ctgs, reference,
#                                    ref_gtf, per_sample=True, threshold=1.0)
#             result.set_index('gene_id', inplace=True)
#
#             assert len(result) == 3
#             assert result.ix['gene_a', 'p_val_corr'] > 0.05
#             assert result.ix['gene_b', 'p_val_corr'] > 0.05
#             assert result.ix['gene_c', 'p_val_corr'] > 0.05
#
#     def test_with_chrom(self, insertions_ctgs, reference, ref_gtf, ref_genes):
#         with mock.patch('imfusion.util.tabix.GtfFrame',
#                         return_value=ref_genes):
#             result = ctg.test_ctgs(insertions_ctgs, reference, ref_gtf,
#                                    chromosomes=['1'], per_sample=True,
#                                    threshold=1.0)
#             assert len(result) == 2
#             assert set(result['gene_id']) == {'gene_a', 'gene_b'}
#
#     def test_with_window(self, insertions_ctgs, reference, ref_gtf, ref_genes):
#         with mock.patch('imfusion.util.tabix.GtfFrame',
#                         return_value=ref_genes):
#             result = ctg.test_ctgs(insertions_ctgs, reference, ref_gtf,
#                                    window=(-2, 0), per_sample=False,
#                                    threshold=1.0)
#             result.set_index('gene_id', inplace=True)
#
#             assert result.ix['gene_a', 'p_val_corr'] < 0.05
