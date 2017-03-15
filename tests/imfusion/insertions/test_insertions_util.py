# -*- coding: utf-8 -*-
"""Tests for imfusion.insertions.util module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import pytest

import numpy as np
import toolz

from imfusion.insertions import util
from imfusion.model import TransposonFusion, Insertion
from imfusion.util.frozendict import frozendict

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def fusion():
    """Example fusion."""
    return TransposonFusion(
        seqname='16',
        anchor_genome=52141095,
        anchor_transposon=1541,
        strand_genome=-1,
        strand_transposon=1,
        flank_genome=-78,
        flank_transposon=-76,
        support_junction=380,
        support_spanning=118,
        metadata=frozendict({}))


@pytest.fixture
def features_path():
    """Returns path to example transposon features."""
    return pytest.helpers.data_path('features.txt', relative_to=__file__)


class TestExtractInsertions(object):
    """Tests extract_insertions function."""

    def test_example(self, fusion, gtf_path, features_path):
        """Test example case of Cblb insertion."""

        insertions = util.extract_insertions([fusion], gtf_path, features_path)
        insertions = list(insertions)

        assert len(insertions) == 1
        assert isinstance(insertions[0].metadata, frozendict)

        assert insertions[0].id == 'INS_1'
        assert insertions[0].seqname == '16'
        assert insertions[0].position == 52141095
        assert insertions[0].strand == -1
        assert insertions[0].support_junction == 380
        assert insertions[0].support_spanning == 118
        assert insertions[0].support == 498
        assert insertions[0].metadata['feature_type'] == 'SA'
        assert insertions[0].metadata['feature_name'] == 'En2SA'
        assert insertions[0].metadata['gene_strand'] == 1
        assert insertions[0].metadata['orientation'] == 'antisense'
        assert insertions[0].metadata['feature_strand'] == -1
        assert insertions[0].metadata['gene_name'] == 'Cblb'
        assert insertions[0].metadata['transposon_anchor'] == 1541
        assert insertions[0].metadata['gene_id'] == 'ENSMUSG00000022637'

        assert isinstance(insertions[0].metadata, frozendict)

    def test_assembly_example(self, rgag1_fusion, gtf_path, features_path,
                              assembled_gtf_path):
        """Tests example case with assembled gtf."""

        insertions = util.extract_insertions(
            [rgag1_fusion],
            gtf_path,
            features_path,
            assembled_gtf_path=assembled_gtf_path)
        insertions = list(insertions)

        assert len(insertions) == 1

        assert insertions[0].id == 'INS_1'
        assert insertions[0].seqname == 'X'
        assert insertions[0].position == 143093898
        assert insertions[0].strand == 1
        assert insertions[0].support_junction == 10
        assert insertions[0].support_spanning == 0
        assert insertions[0].support == 10
        assert insertions[0].metadata['feature_type'] == 'SA'
        assert insertions[0].metadata['feature_name'] == 'En2SA'
        assert insertions[0].metadata['gene_strand'] == 1
        assert insertions[0].metadata['orientation'] == 'sense'
        assert insertions[0].metadata['feature_strand'] == -1
        assert insertions[0].metadata['gene_name'] == 'Rgag1'
        assert insertions[0].metadata['transposon_anchor'] == 1541
        assert insertions[0].metadata['gene_id'] == 'ENSMUSG00000085584'
        assert insertions[0].metadata['novel_transcript'] == 'STRG.14160.1'


class TestAnnotateTransposon(object):
    """Tests annotate_transposon function."""

    def test_example(self, fusion, features_path):
        """Tests example case."""

        annotated = next(
            util.annotate_fusions_for_transposon([fusion], features_path))

        assert annotated.metadata['feature_name'] == 'En2SA'
        assert annotated.metadata['feature_type'] == 'SA'
        assert annotated.metadata['feature_strand'] == -1

    def test_fp_example(self, fusion, features_path):
        """Tests transposon that doesn't overlap transposon feature."""

        fusion = fusion._replace(anchor_transposon=100)
        annotated = next(
            util.annotate_fusions_for_transposon([fusion], features_path))

        assert 'feature_name' not in annotated.metadata
        assert fusion == annotated


@pytest.fixture
def gtf_path():
    """Returns path to example gtf."""
    return pytest.helpers.data_path('reference.gtf.gz', relative_to=__file__)


class TestAnnotateGenes(object):
    """Tests annotate_genes function."""

    def test_example(self, fusion, gtf_path):
        """Tests example case."""

        tr_reference = util.TranscriptReference.from_gtf(gtf_path)

        annotated = util.annotate_fusions_for_genes([fusion], tr_reference)
        annotated_fusion = next(annotated)

        assert annotated_fusion.metadata['gene_name'] == 'Cblb'
        assert annotated_fusion.metadata['gene_strand'] == 1

    def test_fp_example(self, fusion, gtf_path):
        """Tests transposon that doesn't overlap gene."""

        fusion = fusion._replace(anchor_genome=100)
        tr_reference = util.TranscriptReference.from_gtf(gtf_path)

        annotated = util.annotate_fusions_for_genes([fusion], tr_reference)
        annotated_fusion = next(annotated)

        assert 'gene_name' not in annotated_fusion.metadata
        assert fusion == annotated_fusion


@pytest.fixture
def insertion():
    return Insertion(
        id='INS_1',
        seqname='16',
        position=52141095,
        strand=-1,
        support_junction=380,
        support_spanning=118,
        support=498,
        metadata=frozendict({
            'feature_type': 'SA',
            'feature_name': 'En2SA',
            'gene_strand': 1,
            'orientation': 'antisense',
            'feature_strand': -1,
            'gene_name': 'Cblb',
            'transposon_anchor': 1541
        }))


class TestFilterUnexpectedFeatures(object):
    """Tests for filter_unexpected_features."""

    def test_positive_example(self, insertion):
        """Tests example that should not be filtered."""

        filt_ins = list(util.filter_unexpected_features([insertion]))
        assert len(filt_ins) == 1

    def test_negative_example(self, insertion):
        """Tests example that should be filtered."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_type': 'LTR'}))
        filt_ins = list(util.filter_unexpected_features([insertion]))

        assert len(filt_ins) == 0


class TestFilterBlacklist(object):
    """Tests for filter_blacklist."""

    def test_positive_example(self, insertion):
        """Tests example that should not be filtered."""

        filt_ins = list(util.filter_blacklist([insertion], {'Fgfr2'}))
        assert len(filt_ins) == 1

    def test_negative_example(self, insertion):
        """Tests example that should be filtered."""

        filt_ins = list(util.filter_blacklist([insertion], {'Cblb'}))
        assert len(filt_ins) == 0


class TestFilterWrongOrientation(object):
    def test_positive_example(self, insertion):
        """Tests example that should not be filtered."""

        filt_ins = list(util.filter_wrong_orientation([insertion]))
        assert len(filt_ins) == 1

    def test_negative_example(self, insertion):
        """Tests example that should be filtered."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_strand': 1}))
        filt_ins = list(util.filter_wrong_orientation([insertion]))

        assert len(filt_ins) == 0

    def test_na_example_true(self, insertion):
        """Tests example with NA with drop_na = False."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_strand': np.nan}))

        filt_ins = list(
            util.filter_wrong_orientation(
                [insertion], drop_na=False))

        assert len(filt_ins) == 1

    def test_negative_na_example(self, insertion):
        """Tests example with NA with drop_na = True."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_strand': np.nan}))

        filt_ins = list(
            util.filter_wrong_orientation(
                [insertion], drop_na=True))

        assert len(filt_ins) == 0


class TestFilterInsertions(object):
    """Tests for filter_insertions."""

    def test_positive_example(self, insertion):
        """Tests example that should not be filtered."""
        filt_ins = util.filter_insertions(
            [insertion], orientation=True, features=True, blacklist={'Fgfr2'})
        assert len(list(filt_ins)) == 1

    def test_negative_blacklist(self, insertion):
        """Tests example filtered by blacklist."""

        filt_ins = util.filter_insertions([insertion], blacklist={'Cblb'})
        assert len(list(filt_ins)) == 0

    def test_negative_orientation_true(self, insertion):
        """Tests example filtered by orientation."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_strand': 1}))
        filt_ins = util.filter_insertions([insertion], orientation=True)

        assert len(list(filt_ins)) == 0

    def test_negative_orientation_false(self, insertion):
        """Tests if example filtered by orientation is kept."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_strand': 1}))
        filt_ins = util.filter_insertions([insertion], orientation=False)

        assert len(list(filt_ins)) == 1

    def test_negative_feature_true(self, insertion):
        """Tests example filtered by features."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_type': 'LTR'}))
        filt_ins = util.filter_insertions([insertion], features=True)

        assert len(list(filt_ins)) == 0

    def test_negative_feature_false(self, insertion):
        """Tests if example filtered by features is kept."""

        insertion = insertion._replace(metadata=toolz.merge(
            insertion.metadata, {'feature_type': 'LTR'}))
        filt_ins = util.filter_insertions([insertion], features=False)

        assert len(list(filt_ins)) == 1


@pytest.fixture
def rgag1_fusion():
    return TransposonFusion(
        seqname='X',
        anchor_genome=143093898,
        anchor_transposon=1541,
        strand_genome=1,
        strand_transposon=1,
        flank_genome=-78,
        flank_transposon=-76,
        support_junction=10,
        support_spanning=0,
        metadata=frozendict({}))


@pytest.fixture
def assembled_gtf_path():
    """Returns path to example assembled gtf."""
    return pytest.helpers.data_path('assembled.gtf.gz', relative_to=__file__)


class TestAnnotateWithAssembly(object):
    """Tests for annotate_with_assembly function."""

    def test_example(self, rgag1_fusion, gtf_path, assembled_gtf_path):

        gtf_ref = util.TranscriptReference.from_gtf(gtf_path)
        assem_ref = util.TranscriptReference.from_gtf(assembled_gtf_path)

        annotated = next(
            util.annotate_fusions_for_assembly([rgag1_fusion], gtf_ref,
                                               assem_ref))

        assert annotated.metadata['gene_name'] == 'Rgag1'
        assert annotated.metadata['gene_strand'] == 1
        assert annotated.metadata['novel_transcript'] == 'STRG.14160.1'

    def test_neg_example(self, rgag1_fusion, gtf_path, assembled_gtf_path):

        gtf_ref = util.TranscriptReference.from_gtf(gtf_path)
        assem_ref = util.TranscriptReference.from_gtf(assembled_gtf_path)

        rgag1_fusion = rgag1_fusion._replace(anchor_genome=100)

        annotated = list(
            util.annotate_fusions_for_assembly([rgag1_fusion], gtf_ref,
                                               assem_ref))

        assert len(annotated) == 1
        assert 'novel_transcript' not in annotated[0].metadata
        assert 'gene_name' not in annotated[0].metadata
