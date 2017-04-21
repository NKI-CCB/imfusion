# -*- coding: utf-8 -*-
"""Tests for imfusion.build.indexers.tophat module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from pathlib2 import Path

import pyfaidx
import pytest

from imfusion.build.indexers import tophat


class TestTophatIndexer(object):
    """Tests for TophatIndexer class."""

    def test_dependencies(self):
        """Check dependencies."""

        indexer = tophat.TophatIndexer()
        assert set(indexer.dependencies) == {'bowtie-build', 'tophat2'}

    def test_build(self, build_kws, mocker):
        """Tests build using example files."""

        # Mock Tophat call.
        mock_bowtie = mocker.patch.object(tophat, 'bowtie_index')
        mock_tophat = mocker.patch.object(tophat,
                                          'tophat2_index_transcriptome')

        # Build reference.
        indexer = tophat.TophatIndexer()
        indexer.build(**build_kws)

        # Check if reference files exist.
        ref = tophat.TophatReference(build_kws['output_dir'])

        assert ref.base_path.exists()
        assert ref.fasta_path.exists()
        assert ref.gtf_path.exists()
        assert ref.indexed_gtf_path.exists()
        # assert ref.index_path.exists()
        assert ref.transposon_name == 'T2onc'
        assert ref.transposon_path.exists()
        assert ref.features_path.exists()

        # Check presence of augmented reference sequences.
        refseq = pyfaidx.Fasta(str(ref.fasta_path))
        assert sorted(refseq.keys()) == ['1', '2', 'T2onc']

        # Check call to Bowtie for building the index.
        mock_bowtie.assert_called_once_with(
            reference_path=ref.fasta_path,
            output_base_path=ref.index_path,
            log_path=build_kws['output_dir'] / 'bowtie.log')

        # Check call to Tophat for building the index.
        mock_tophat.assert_called_once_with(
            bowtie_index_path=ref.index_path,
            gtf_path=ref.gtf_path,
            output_base_path=ref.transcriptome_path,
            log_path=build_kws['output_dir'] / 'transcriptome.log')
