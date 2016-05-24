# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import pytest

import pandas as pd
import numpy as np

from imfusion.util import fusions as mod
from imfusion.model import Fusion
from imfusion.util.tabix import GtfFile

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def fusion():
    """Fixture for example fusion."""
    return Fusion(seqname='1', anchor_genome=182453750, anchor_transposon=120,
                  flank_genome=10, flank_transposon=10,
                  strand_genome=1, strand_transposon=-1)

@pytest.fixture
def transposon_features():
    """Fixture for example transposon features."""
    return pd.DataFrame({
        'name': ['En2SA', 'SD', 'SD2'],
        'start': [100, 600, 650],
        'end': [200, 700, 750],
        'strand': [-1, +1, -1],
        'type': ['SA', 'SD', 'SD']
    })


@pytest.fixture
def gtf():
    """Fixture for gtf file path."""
    return GtfFile(pytest.helpers.data_path('mm10.test.gtf.gz'))


class TestAnnotateFusions(object):
    """Tests annotate_fusions function."""

    def test_simple(self, fusion, gtf, transposon_features):
        """Tests behaviour of default annotation."""

        annotated = list(mod.annotate_fusions(
            [fusion], gtf, transposon_features))
        assert len(annotated) == 1

        assert annotated[0].gene_id == 'ENSMUSG00000026510'
        assert annotated[0].gene_name == 'Trp53bp2'
        assert annotated[0].gene_strand == 1

        assert annotated[0].feature_name == 'En2SA'
        assert annotated[0].feature_type == 'SA'
        assert annotated[0].feature_strand == -1


class TestAnnotateWithGene(object):
    """Tests annotate_with_gene function."""

    # TODO: test with flank genome = nan.
    # TODO: test with overlap only from flank.

    def test_single(self, fusion, gtf):
        """Test single overlapping gene."""

        annotated = list(mod.annotate_with_gene([fusion], gtf))

        assert len(annotated) == 1
        assert annotated[0].gene_id == 'ENSMUSG00000026510'
        assert annotated[0].gene_name == 'Trp53bp2'
        assert annotated[0].gene_strand == 1

    def test_different_ref(self, fusion, gtf):
        """Test same position, different chromosome."""

        fusion = fusion._replace(seqname='11')
        annotated = list(mod.annotate_with_gene([fusion], gtf))

        assert len(annotated) == 1
        assert annotated[0].gene_id is None
        assert annotated[0].gene_name is None

    def test_none(self, fusion, gtf):
        """Test case with no overlapping gene."""

        fusion = fusion._replace(anchor_genome=1000)
        annotated = list(mod.annotate_with_gene([fusion], gtf))

        assert len(annotated) == 1
        assert annotated[0].gene_id is None
        assert annotated[0].gene_name is None


class TestAnnotateWithTransposon(object):
    """Tests annotate_with_transposon function."""

    # TODO: test with overlap only from flank.

    def test_single(self, fusion, transposon_features):
        """Tests case where fusion overlaps with single feature."""

        annotated = list(mod.annotate_with_transposon(
            [fusion], transposon_features))

        assert len(annotated) == 1
        assert annotated[0].feature_name == 'En2SA'
        assert annotated[0].feature_type == 'SA'
        assert annotated[0].feature_strand == -1

    def test_multiple(self, fusion, transposon_features):
        """Tests case where fusion overlaps with multiple features."""

        fusion = fusion._replace(anchor_transposon=670)
        annotated = list(mod.annotate_with_transposon(
            [fusion], transposon_features))

        assert len(annotated) == 2

    def test_none(self, fusion, transposon_features):
        """Tests case where fusion does not overlap with feature."""

        fusion = fusion._replace(anchor_transposon=10)
        annotated = list(mod.annotate_with_transposon(
            [fusion], transposon_features))

        assert len(annotated) == 1
        assert annotated[0].feature_name is None
        assert annotated[0].feature_type is None


@pytest.fixture
def fusion2():
    return Fusion(seqname='1', anchor_genome=182436900, anchor_transposon=420,
                  strand_genome=1, strand_transposon=1, flank_genome=35,
                  flank_transposon=35, spanning_reads=15, supporting_mates=0,
                  spanning_mates=0, gene_id='ENSMUSG00000026510',
                  gene_name='Trp53bp2', gene_strand=1, feature_name='En2SA',
                  feature_type='SA', feature_strand=-1)


class TestInsertionPosition(object):
    """Tests for insertion_position."""

    def test_sa(self, fusion2, gtf):
        """Test splice acceptor (should be placed downstream)."""

        position = mod.insertion_position(
            fusion2, gtf, strand=1, orientation='sense')

        # Check if downstream of exon and anchor.
        assert position > 182436999
        assert position > fusion2.anchor_genome

    def test_sd(self, fusion2, gtf):
        """Test splice donor (should be placed upstream)."""

        # Convert to SD.
        fusion2 = fusion2._replace(feature_name='SD')
        position = mod.insertion_position(
            fusion2, gtf, strand=1, orientation='sense')

        # Check if upstream of exon and anchor.
        assert position < 182436897
        assert position < fusion2.anchor_genome

    def test_sa_reverse(self, fusion2, gtf):
        """Test splice acceptor on -1 strand (should be placed downstream)."""

        # Convert to SD.
        fusion2 = fusion2._replace(feature_strand=-1, strand_transposon=-1)
        position = mod.insertion_position(
            fusion2, gtf, strand=-1, orientation='antisense')

        # Check if downstream of exon and anchor.
        assert position > 182436999
        assert position > fusion2.anchor_genome


class TestInsertionOrientation(object):
    """Tests for insertion_orientation."""

    def test_sense(self, fusion2):
        """Test sense insertion."""

        strand, orientation = mod.insertion_orientation(fusion2)

        assert strand == 1
        assert orientation == 'sense'

    def test_antisense(self, fusion2):
        """Test antisense insertion."""

        fusion2 = fusion2._replace(strand_transposon=-1)
        strand, orientation = mod.insertion_orientation(fusion2)

        assert strand == -1
        assert orientation == 'antisense'

    def test_no_gene(self, fusion2):
        """Test sense insertion without gene."""

        fusion2 = fusion2._replace(gene_strand=np.nan)
        strand, orientation = mod.insertion_orientation(fusion2)

        assert strand == 1
        assert orientation is None


class TestPlaceFusions(object):
    """Tests for place_fusions."""

    def test_example(self, fusion2, gtf):
        """Tests conversion of example fusion to insertion."""

        insertions = list(mod.place_fusions(
            [fusion2], sample_id='s1', reference_gtf=gtf))

        # Check if we have our insertion.
        assert len(insertions) == 1

        # Check insertion attributes.
        ins = insertions[0]
        assert ins.seqname == fusion2.seqname
        assert ins.position > fusion2.anchor_genome
        assert ins.strand == 1
        assert ins.sample_id == 's1'
        assert ins.gene_id == fusion2.gene_id
        assert ins.gene_name == fusion2.gene_name
        assert ins.gene_strand == fusion2.gene_strand
        assert ins.orientation == 'sense'
        assert ins.feature_name == 'En2SA'
        assert ins.feature_type == 'SA'
        assert ins.feature_strand == -1
        assert ins.anchor_genome == fusion2.anchor_genome
        assert ins.anchor_transposon == fusion2.anchor_transposon
        assert ins.flank_genome == fusion2.flank_genome
        assert ins.flank_transposon == fusion2.flank_transposon
        assert ins.spanning_reads == fusion2.spanning_reads
        assert ins.supporting_mates == fusion2.supporting_mates
        assert ins.spanning_mates == fusion2.spanning_mates
