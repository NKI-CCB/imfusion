# -*- coding: utf-8 -*-
"""Unit tests for functions calling Tophat externally."""

from pathlib2 import Path

import pytest

from imfusion.external import tophat

# pylint: disable=no-self-use,redefined-outer-name


class TestTophat2IndexTranscriptome(object):
    """Unit tests for the tophat2_index_transcriptome function."""

    def test_call(self, mocker, tmpdir):
        """Tests example call."""

        mock_run = mocker.patch.object(tophat, 'run_command')

        tmp_dir = Path(str(tmpdir / '_tmp'))

        tophat.tophat2_index_transcriptome(
            Path('index'),
            gtf_path=Path('reference.gtf'),
            output_base_path='transcriptome',
            tmp_dir=tmp_dir)

        mock_run.assert_called_once_with(
            args=[
                'tophat2', '--GTF', 'reference.gtf',
                '--transcriptome-index=transcriptome', '--bowtie1',
                '--output-dir', str(tmp_dir), 'index'
            ],
            log_path=None)


@pytest.fixture
def tophat2_align_kws(tmpdir):
    """Example kws for tophat2_align."""

    output_dir = Path(str(tmpdir / 'out'))

    return {
        'fastq_path': Path('sample.R1.fastq.gz'),
        'index_path': Path('index'),
        'output_dir': output_dir
    }


class TestTophat2Align(object):
    """Unit tests for the tophat2_align function."""

    def test_basic_call(self, mocker, tophat2_align_kws):
        """Tests example call."""

        mock_run = mocker.patch.object(tophat, 'run_command')

        tophat.tophat2_align(**tophat2_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'tophat2', '--output-dir',
                str(tophat2_align_kws['output_dir']), 'index',
                'sample.R1.fastq.gz'
            ],
            log_path=None)

    def test_paired_end(self, mocker, tophat2_align_kws):
        """Tests example call with paired-end data."""

        mock_run = mocker.patch.object(tophat, 'run_command')

        tophat2_align_kws['fastq2_path'] = 'sample.R2.fastq.gz'
        tophat.tophat2_align(**tophat2_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'tophat2', '--output-dir',
                str(tophat2_align_kws['output_dir']), 'index',
                'sample.R1.fastq.gz', 'sample.R2.fastq.gz'
            ],
            log_path=None)

    def test_extra_arguments(self, mocker, tophat2_align_kws):
        """Tests example call with extra arguments."""

        mock_run = mocker.patch.object(tophat, 'run_command')

        tophat2_align_kws['extra_args'] = {'--GTF': 'reference.gtf'}
        tophat.tophat2_align(**tophat2_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'tophat2', '--GTF', 'reference.gtf', '--output-dir',
                str(tophat2_align_kws['output_dir']), 'index',
                'sample.R1.fastq.gz'
            ],
            log_path=None)
