# -*- coding: utf-8 -*-
"""Tests for imfusion.insertions.aligners.tophat module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse
import shutil

from pathlib2 import Path

from future.utils import native_str
import pytest

from imfusion.insertions.aligners import tophat
from imfusion.model import Insertion
from imfusion.util.frozendict import frozendict

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def read_paths():
    """Example read paths."""
    return (Path('a.fastq.gz'), Path('b.fastq.gz'))


@pytest.fixture
def tophat_path():
    """Return path to fusions.out file."""
    return pytest.helpers.data_path('fusions.out', relative_to=__file__)


@pytest.fixture
def tophat_reference():
    """Returns example reference for Tophat2."""
    ref_path = pytest.helpers.data_path('reference', relative_to=__file__)
    return tophat.TophatReference(ref_path)


@pytest.fixture
def tophat_output_dir(tmpdir, tophat_path):
    """Simulates Tophat2 output directory."""

    # Create directories.
    output_dir = Path(native_str(tmpdir / 'out'))

    tophat_dir = output_dir / '_tophat'
    tophat_dir.mkdir(parents=True)

    # Copy / simulate aligner output files.
    shutil.copy(native_str(tophat_path),
                native_str(tophat_dir / 'fusions.out'))  # yapf: disable

    pytest.helpers.touch(tophat_dir / 'Aligned.sortedByCoord.out.bam')

    return output_dir


@pytest.fixture
def cmdline_args(tmpdir):
    """Example command line arguments."""
    return [
        '--fastq', 'a.fastq.gz', '--fastq2', 'b.fastq.gz', '--reference',
        native_str(tmpdir), '--output_dir', '/path/to/out'
    ]


class TestTophatAligner(object):
    """Tests for the TophatAligner class."""

    def test_dependencies(self, tophat_reference):
        """Test dependencies for various configurations."""

        # Basic call.
        aligner = tophat.TophatAligner(tophat_reference)
        assert set(aligner.dependencies) == {'tophat2', 'bowtie'}

        # With assembly.
        aligner = tophat.TophatAligner(tophat_reference, assemble=True)
        assert set(aligner.dependencies) == {'tophat2', 'bowtie', 'stringtie'}

    def test_identify_insertions(self, read_paths, tophat_reference,
                                 tophat_output_dir, mocker):
        """Integration test using data from an example run."""

        # TODO: Check assembly?

        # Mock star_align call.
        tophat_mock = mocker.patch.object(tophat, 'tophat2_align')
        mocker.patch.object(tophat.util, 'count_lines', return_value=8e6)

        # Call identify insertions.
        fastq, fastq2 = read_paths

        aligner = tophat.TophatAligner(tophat_reference)
        ins = list(
            aligner.identify_insertions(
                fastq, tophat_output_dir, fastq2_path=fastq2))

        # Check call to star_align.

        tophat_mock.assert_called_once_with(
            fastq_path=fastq,
            fastq2_path=fastq2,
            output_dir=tophat_output_dir / '_tophat',
            index_path=tophat_reference.index_path,
            extra_args={
                '--fusion-search': (),
                '--transcriptome-index':
                (native_str(tophat_reference.transcriptome_path), ),
                '--num-threads': (1, ),
                '--bowtie1': (),
                '--fusion-anchor-length': (12, )
            })

        # Check result, including specific Cblb insertion.
        assert len(ins) == 7

        assert ins[2].id == 'INS_4'
        assert ins[2].seqname == '16'
        assert ins[2].position == 52141093
        assert ins[2].strand == -1
        assert ins[2].support_junction == 462
        assert ins[2].support_spanning == 103
        assert ins[2].support == 565
        assert ins[2].metadata['gene_id'] == 'ENSMUSG00000022637'
        assert ins[2].metadata['transposon_anchor'] == 1539
        assert ins[2].metadata['feature_name'] == 'En2SA'
        assert ins[2].metadata['gene_name'] == 'Cblb'
        assert ins[2].metadata['feature_type'] == 'SA'
        assert ins[2].metadata['gene_strand'] == 1
        assert ins[2].metadata['orientation'] == 'antisense'
        assert ins[2].metadata['ffpm_junction'] == 231.0
        assert ins[2].metadata['ffpm_spanning'] == 51.5
        assert ins[2].metadata['ffpm'] == 282.5

    def test_from_args_basic(self, cmdline_args):
        """Tests creation with minimal arguments."""

        # Setup parser.
        parser = argparse.ArgumentParser()
        tophat.TophatAligner.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        aligner = tophat.TophatAligner.from_args(args)

        # Check args.
        assert args.fastq == Path('a.fastq.gz')
        assert args.fastq2 == Path('b.fastq.gz')
        assert args.output_dir == Path('/path/to/out')

        # Check aligner.
        # pylint: disable=w0212
        assert not aligner._assemble
        assert aligner._assemble_args == {}
        assert aligner._min_flank == 12
        assert aligner._threads == 1
        assert aligner._extra_args == {}
        assert aligner._filter_features
        assert aligner._filter_orientation
        assert aligner._filter_blacklist is None

    def test_from_args_extra_args(self, cmdline_args):
        """Tests creation with extra options."""

        cmdline_args += [
            '--tophat_threads', '5', '--tophat_min_flank', '20',
            '--tophat_args', "--limitBAMsortRAM 2000", '--assemble',
            '--no_filter_orientation', '--no_filter_feature',
            '--blacklisted_genes', 'En2'
        ] # yapf:disable

        # Setup parser.
        parser = argparse.ArgumentParser()
        tophat.TophatAligner.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        aligner = tophat.TophatAligner.from_args(args)

        # Check args.
        assert args.fastq == Path('a.fastq.gz')
        assert args.fastq2 == Path('b.fastq.gz')
        assert args.output_dir == Path('/path/to/out')

        # Check aligner.
        # pylint: disable=w0212
        assert aligner._assemble
        assert aligner._assemble_args == {}
        assert aligner._min_flank == 20
        assert aligner._threads == 5
        assert aligner._extra_args == {'--limitBAMsortRAM': ('2000', )}
        assert not aligner._filter_features
        assert not aligner._filter_orientation
        assert aligner._filter_blacklist == ['En2']
