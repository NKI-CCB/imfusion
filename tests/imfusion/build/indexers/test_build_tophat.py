# -*- coding: utf-8 -*-
"""Tests for imfusion.build.indexers.tophat module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

try:
    from pathlib import Path
except ImportError:
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
        mock = mocker.patch.object(tophat.shell, 'run_command')

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
        assert mock.call_args_list[0][1]['args'] == [
            'bowtie-build', str(ref.fasta_path), str(ref.index_path)
        ]

        # Check call to Tophat for building the index.
        assert mock.call_args_list[1][1]['args'] == [
            'tophat2', '--GTF', str(ref.gtf_path), '--transcriptome-index=' +
            str(ref.transcriptome_path), '--bowtie1', '--output-dir',
            pytest.helpers.mock_any(str), str(ref.index_path)
        ]
