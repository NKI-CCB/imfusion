# -*- coding: utf-8 -*-
"""Unit tests for functions calling bowtie externally."""

from pathlib2 import Path

import pytest

from imfusion.external import star

# pylint: disable=no-self-use,redefined-outer-name


class TestStarIndex(object):
    """Unit tests for the star_index function."""

    def test_call(self, mocker, tmpdir):
        """Tests example call."""

        mock_run = mocker.patch.object(star, 'run_command')

        output_dir = Path(str(tmpdir / 'out'))

        star.star_index(
            Path('reference.fa'),
            gtf_path=Path('reference.gtf'),
            output_dir=output_dir,
            overhang=50)

        mock_run.assert_called_once_with(
            args=[
                'STAR', '--runMode', 'genomeGenerate', '--genomeDir',
                str(output_dir), '--genomeFastaFiles', 'reference.fa',
                '--sjdbGTFfile', 'reference.gtf', '--sjdbOverhang', '50',
                '--runThreadN', '1'
            ],
            log_path=None)


@pytest.fixture
def star_align_kws(tmpdir):
    """Example kws for star_align."""

    output_dir = Path(str(tmpdir / 'out'))

    return {
        'fastq_path': Path('sample.R1.fastq.gz'),
        'index_path': Path('index'),
        'output_dir': output_dir
    }


class TestStarAlign(object):
    """Unit tests for the star_align function."""

    def test_basic_call(self, mocker, star_align_kws):
        """Tests example call."""

        mock_run = mocker.patch.object(star, 'run_command')

        star.star_align(**star_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'STAR', '--genomeDir', 'index', '--outFileNamePrefix',
                str(star_align_kws['output_dir']) + '/', '--readFilesIn',
                'sample.R1.fastq.gz', '--readFilesCommand', 'gunzip', '-c'
            ],
            log_path=None)

    def test_basic_without_gzip(self, mocker, star_align_kws):
        """Tests example without gzipped input."""

        mock_run = mocker.patch.object(star, 'run_command')

        star_align_kws['fastq_path'] = Path('sample.R1.fastq')
        star.star_align(**star_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'STAR', '--genomeDir', 'index', '--outFileNamePrefix',
                str(star_align_kws['output_dir']) + '/', '--readFilesIn',
                'sample.R1.fastq'
            ],
            log_path=None)

    def test_paired_end(self, mocker, star_align_kws):
        """Tests example call with paired-end data."""

        mock_run = mocker.patch.object(star, 'run_command')

        star_align_kws['fastq2_path'] = 'sample.R2.fastq.gz'
        star.star_align(**star_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'STAR', '--genomeDir', 'index', '--outFileNamePrefix',
                str(star_align_kws['output_dir']) + '/', '--readFilesIn',
                'sample.R1.fastq.gz', 'sample.R2.fastq.gz',
                '--readFilesCommand', 'gunzip', '-c'
            ],
            log_path=None)

    def test_extra_arguments(self, mocker, star_align_kws):
        """Tests example call with extra arguments."""

        mock_run = mocker.patch.object(star, 'run_command')

        star_align_kws['extra_args'] = {'--outSAMtype': 'Unsorted'}
        star.star_align(**star_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'STAR', '--genomeDir', 'index', '--outFileNamePrefix',
                str(star_align_kws['output_dir']) + '/', '--readFilesIn',
                'sample.R1.fastq.gz', '--outSAMtype', 'Unsorted',
                '--readFilesCommand', 'gunzip', '-c'
            ],
            log_path=None)

    def test_multiple_threads(self, mocker, star_align_kws):
        """Tests example call with extra arguments."""

        mock_run = mocker.patch.object(star, 'run_command')

        star_align_kws['threads'] = 5
        star.star_align(**star_align_kws)

        mock_run.assert_called_once_with(
            args=[
                'STAR', '--genomeDir', 'index', '--outFileNamePrefix',
                str(star_align_kws['output_dir']) + '/', '--readFilesIn',
                'sample.R1.fastq.gz', '--readFilesCommand', 'gunzip', '-c',
                '--runThreadN', '5'
            ],
            log_path=None)
