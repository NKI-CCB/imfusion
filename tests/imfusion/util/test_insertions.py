# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import pytest

import numpy as np

from imfusion.util import insertions as _mod
from imfusion.model import Insertion

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods

@pytest.fixture
def ori_insertions():
    """Example insertions for orientation test."""

    return [
        # Correct insertion.
        Insertion(id='O1', seqname='1', position=200, strand=1, sample_id='s1',
                  gene_strand=1, feature_name='En2SA', feature_strand=1),
        # Insertion with wrong orientation.
        Insertion(id='O2', seqname='1', position=200, strand=1, sample_id='s1',
                  gene_strand=1, feature_name='SD', feature_strand=-1),
        # Insertion with NA orientation.
        Insertion(id='O3', seqname='1', position=200, strand=1, sample_id='s1',
                  gene_strand=np.nan, feature_name='SD', feature_strand=-1)]


class TestFilterWrongOrientation(object):
    """Tests for filter_wrong_orientation."""

    def test_simple(self, ori_insertions):
        """Tests simple case, keeping na's."""

        filt = list(_mod.filter_wrong_orientation(ori_insertions))

        assert len(filt) == 2
        assert [i.id for i in filt] == ['O1', 'O3']

    def test_na(self, ori_insertions):
        """Tests simple case, dropping na's."""

        filt = list(_mod.filter_wrong_orientation(
            ori_insertions, drop_na=True))

        assert len(filt) == 1
        assert [i.id for i in filt] == ['O1']


@pytest.fixture
def site_insertions():
    """Example insertions for feature name test."""

    return [
        # Correct insertion.
        Insertion(id='S1', seqname='1', position=200, strand=1, sample_id='s1',
                  gene_strand=1, feature_name='En2SA', feature_strand=1),
        # Insertion with wrong feature.
        Insertion(id='S2', seqname='1', position=200, strand=1, sample_id='s1',
                  gene_strand=1, feature_name='MSCV', feature_strand=1)]


class TestFilterUnexpectedSites(object):
    """Tests for filter_unexpected_sites."""

    def test_example(self, site_insertions):
        """Tests example insertions in filter_unexpected_sites."""

        filt = list(_mod.filter_unexpected_sites(site_insertions))

        assert len(filt) == 1
        assert [i.id for i in filt] == ['S1']


@pytest.fixture
def insertions(ori_insertions, site_insertions):
    """Example insertions for combined test."""
    return ori_insertions + site_insertions


class TestFilterInsertions(object):
    """Tests for filter_insertions."""

    def test_example(self, insertions):
        """Tests example insertions in filter_insertions."""

        filt = list(_mod.filter_invalid_insertions(insertions))

        assert len(filt) == 3
        assert [i.id for i in filt] == ['O1', 'O3', 'S1']


@pytest.fixture
def gtf_path():
    """Fixture for gtf file path."""
    return pytest.helpers.data_path('mm10.test.gtf.gz')


@pytest.fixture
def insertions2():
    return [Insertion(id='A', seqname='1', position=10, sample_id='s1',
                      strand=1, gene_id='ENSMUSG00000026510'),
            Insertion(id='B', seqname='1', position=182462432,
                      strand=1, sample_id='s1',),
            Insertion(id='C', seqname='1', position=20, strand=1,
                      sample_id='s1', gene_id='ENSMUSG00000000001')]


class TestFilterBlacklist(object):

    def test_example(self, insertions2, gtf_path):
        # Apply filter.
        filt = list(_mod.filter_blacklist(
            insertions2, gene_ids=['ENSMUSG00000026510'],
            reference_gtf=gtf_path, filter_overlap=True))

        # Check result.
        assert len(filt) == 1
        assert [i.id for i in filt] == ['C']

    def test_no_overlap(self, insertions2, gtf_path):
        # Apply filter.
        filt = list(_mod.filter_blacklist(
            insertions2, gene_ids=['ENSMUSG00000026510'],
            reference_gtf=gtf_path, filter_overlap=False))

        # Check result.
        assert len(filt) == 2
        assert [i.id for i in filt] == ['B', 'C']

    def test_require_gtf(self, insertions2):
        with pytest.raises(ValueError):
            list(_mod.filter_blacklist(
                insertions2, gene_ids=['ENSMUSG00000026510']))
