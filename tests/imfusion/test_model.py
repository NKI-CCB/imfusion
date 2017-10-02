# -*- coding: utf-8 -*-
"""Tests for imfusion.model module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import pytest

from imfusion.model import Fusion, TransposonFusion, Insertion
from imfusion.vendor.frozendict import frozendict

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def fusion():
    """Example fusion."""
    return Fusion(
        chromosome_a='1',
        position_a=300,
        strand_a=1,
        chromosome_b='T2onc',
        position_b=420,
        strand_b=1,
        flank_a=52,
        flank_b=62,
        support_junction=2,
        support_spanning=1,
        sample='S1')


@pytest.fixture
def tr_fusion():
    """Example transposon fusion."""
    return TransposonFusion(
        chromosome='16',
        anchor_genome=52141095,
        anchor_transposon=1541,
        strand_genome=-1,
        strand_transposon=1,
        flank_genome=-78,
        flank_transposon=-76,
        support_junction=380,
        support_spanning=118,
        sample='S1',
        gene_strand=1,
        gene_name='Cblb',
        feature_type='SA',
        feature_strand=-1,
        feature_name='En2SA')


class TestTransposonFusion(object):
    """Tests for TransposonFusion methods/properties."""

    def test_support(self, tr_fusion):
        """Tests support property."""
        expected = tr_fusion.support_junction + tr_fusion.support_spanning
        assert tr_fusion.support == expected

    def test_transposon_region(self, tr_fusion):
        """Tests transposon_region property."""
        assert tr_fusion.transposon_region == (1465, 1541)

    def test_genome_region(self, tr_fusion):
        """Tests genome_region property."""
        assert tr_fusion.genome_region == ('16', 52141017, 52141095)

    def test_from_fusion(self, fusion):
        """Tests from_fusion method."""

        tr_fusion = TransposonFusion.from_fusion(fusion, 'T2onc')

        expected = TransposonFusion(
            chromosome='1',
            anchor_genome=300,
            anchor_transposon=420,
            strand_genome=1,
            strand_transposon=1,
            flank_genome=-52,
            flank_transposon=62,
            support_junction=2,
            support_spanning=1,
            sample='S1')

        assert tr_fusion == expected
        assert isinstance(tr_fusion.metadata, frozendict)


class TestInsertion(object):
    """Tests methods/properties from Insertion class."""

    def test_from_transposon_fusion(self, tr_fusion):
        """Tests from_tranpsoson_fusion for example."""

        insertion = Insertion.from_transposon_fusion(tr_fusion)

        expected = Insertion(
            id=None,
            chromosome='16',
            position=52141095,
            strand=-1,
            support=498,
            sample='S1',
            transposon_anchor=1541,
            orientation='antisense',
            gene_name='Cblb',
            feature_type='SA',
            feature_strand=-1,
            feature_name='En2SA',
            gene_strand=1,
            support_junction=380,
            support_spanning=118)

        assert insertion == expected
        assert isinstance(insertion.metadata, frozendict)

    def test_from_transposon_fusion_id(self, tr_fusion):
        """Tests from_tranpsoson_fusion for example, with id."""

        insertion = Insertion.from_transposon_fusion(tr_fusion, id='a')
        assert insertion.id == 'a'

    def test_from_transposon_fusions(self, tr_fusion):
        """Tests from_tranpsoson_fusions for example."""

        insertions = list(Insertion.from_transposon_fusions([tr_fusion]))

        expected = [
            Insertion(
                id=None,
                chromosome='16',
                position=52141095,
                strand=-1,
                support=498,
                sample='S1',
                transposon_anchor=1541,
                orientation='antisense',
                gene_name='Cblb',
                feature_type='SA',
                feature_strand=-1,
                feature_name='En2SA',
                gene_strand=1,
                support_junction=380,
                support_spanning=118)
        ]

        assert expected == insertions

    def test_from_transposon_fusions_id(self, tr_fusion):
        """Tests from_tranpsoson_fusions for example, with id."""

        insertions = list(
            Insertion.from_transposon_fusions(
                [tr_fusion], id_fmt_str='INS_{}'))

        assert insertions[0].id == 'INS_1'
