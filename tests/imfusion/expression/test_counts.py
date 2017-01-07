"""Tests for the imfusion.expression.counts module."""

import shutil

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pytest

from imfusion.expression import counts

# pylint: disable=no-self-use,redefined-outer-name

# Disable E1101 checks which stumble on pandas classes.
# pylint: disable=E1101


@pytest.fixture
def fc_counts(tmpdir):
    """Example exon count output file from featureCounts."""

    # Copy file into place.
    src = pytest.helpers.data_path('fc_counts.txt', relative_to=__file__)
    dest = tmpdir / 'counts.txt'
    shutil.copy(str(src), str(dest))

    # Return copied file path.
    return dest


@pytest.fixture
def count_kws(tmpdir):
    """Basic keyword arguments for generate_exon_counts."""

    return {
        'bam_files': [Path('a.bam'), Path('b.bam')],
        'gtf_path': Path('/path/to/gtf'),
        'names': {
            'a.bam': 'a',
            'b.bam': 'b'
        },
        'tmp_dir': tmpdir
    }


class TestGenerateExonCounts(object):
    """Tests for the generate_exon_counts function."""

    def test_example(self, count_kws, fc_counts, mocker):
        """Tests function using example output file from feature counts."""

        # Mock and copy result into place.
        mock = mocker.patch.object(counts.subprocess, 'check_output')

        # Call function.
        exon_counts = counts.generate_exon_counts(**count_kws)

        # Check call + result.
        mock.assert_called_once_with([
            'featureCounts', '--minOverlap', '1', '-O', '-f', '-t',
            'exonic_part', '-a', '/path/to/gtf', '-o', str(fc_counts), 'a.bam',
            'b.bam'
        ])

        assert exon_counts.shape == (7, 2)
        assert list(exon_counts.columns) == ['a', 'b']


@pytest.fixture
def exon_counts_path():
    """Example exon counts file."""
    return pytest.helpers.data_path('exon_counts.txt', relative_to=__file__)


class TestReadExonCounts(object):
    """Tests for the read_exon_counts function."""

    def test_example(self, exon_counts_path):
        """Tests reading example output from the generate_exon_counts test."""

        exon_counts = counts.read_exon_counts(exon_counts_path)

        assert exon_counts.shape == (7, 2)
        assert list(exon_counts.columns) == ['a', 'b']

    def test_example_with_gene_id(self, exon_counts_path):
        """Tests filtering for a given gene id."""

        exon_counts = counts.read_exon_counts(
            exon_counts_path, gene_id='ENSMUSG00000051951')

        assert exon_counts.shape == (5, 2)
        assert list(exon_counts.columns) == ['a', 'b']
