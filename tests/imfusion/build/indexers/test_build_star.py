# -*- coding: utf-8 -*-
"""Tests for imfusion.build.indexers.star module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse

import pyfaidx
import pytest

from pathlib2 import Path

from imfusion.build.indexers import star

# pylint: disable=no-self-use,redefined-outer-name


class TestStarIndexer(object):
    """Tests for StarIndexer class."""

    def test_dependencies(self):
        """Check dependencies."""

        indexer = star.StarIndexer()
        assert indexer.dependencies == ['STAR']

    def test_build(self, build_kws, mocker):
        """Tests build using example files."""

        # Mock STAR call.
        mock = mocker.patch.object(star, 'star_index')

        # Build reference.
        indexer = star.StarIndexer()
        indexer.build(**build_kws)

        # Check if reference files exist.
        ref = star.StarReference(build_kws['output_dir'])

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

        # Check call to STAR for building the index.
        mock.assert_called_once_with(
            fasta_path=ref.fasta_path,
            gtf_path=ref.gtf_path,
            output_dir=ref.index_path,
            log_path=build_kws['output_dir'] / 'star.log',
            overhang=100,
            threads=1)

    def test_from_args(self, cmdline_args):
        """Tests creation from command line."""

        # Setup parser.
        parser = argparse.ArgumentParser()
        star.StarIndexer.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        indexer = star.StarIndexer.from_args(args)

        # Check args.
        assert args.reference_seq == Path('/path/to/ref')
        assert args.reference_gtf == Path('/path/to/gtf')
        assert args.transposon_seq == Path('/path/to/tr')
        assert args.transposon_features == Path('/path/to/feat')

        # Check instance.
        assert indexer._overhang == 100
        assert indexer._threads == 1

    def test_from_args_extra(self, cmdline_args):
        """Tests creation from command line using extra args."""

        cmdline_args += ['--star_overhang', '50', '--star_threads', '10']

        # Setup parser.
        parser = argparse.ArgumentParser()
        star.StarIndexer.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        indexer = star.StarIndexer.from_args(args)

        # Check instance.
        assert indexer._overhang == 50
        assert indexer._threads == 10


class TestStarReference(object):
    """Tests for StarReference class."""

    def test_non_existing(self, tmpdir):
        """Tests non-existing reference"""

        with pytest.raises(ValueError):
            star.StarReference(tmpdir / 'test')
