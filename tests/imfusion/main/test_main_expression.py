# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import shutil

import mock
import pytest

import pandas as pd

from imfusion.main import main
from imfusion.main import expression as expr

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def insertion_dir():
    """Path to insertion directory."""
    return pytest.helpers.data_path('insertions/s1', relative_to=__file__)


@pytest.fixture
def mapping_path():
    """Path to sample mapping."""
    return pytest.helpers.data_path('mapping.txt', relative_to=__file__)


@pytest.fixture
def cmdline_args(insertion_dir, tmpdir):
    """Basic command line args."""

    tmp_path = tmpdir / 's1'
    shutil.copytree(str(insertion_dir), str(tmp_path))

    return ['expression',
            '--sample_dir', str(tmp_path),
            '--exon_gtf', '/path/to/flattened_gtf']


@pytest.fixture
def gene_counts():
    """Gene counts example."""
    counts_path = pytest.helpers.data_path(
        'gene_counts.txt', relative_to=__file__)
    return pd.read_csv(str(counts_path), sep='\t').set_index('gene_id')


@pytest.fixture
def exon_counts():
    """Exon counts example."""
    counts_path = pytest.helpers.data_path(
        'exon_counts.txt', relative_to=__file__)
    index_cols = ['gene_id', 'chr', 'start', 'end', 'strand']
    return pd.read_csv(str(counts_path), sep='\t').set_index(index_cols)


@pytest.fixture
def bam_paths():
    """Bam paths for examples (used to verify call)."""
    template = 'insertions/%s/tophat/accepted_hits.bam'
    bam_paths = {s: Path(pytest.helpers.data_path(template % s,
                                                  relative_to=__file__))
                 for s in {'s1', 's2'}}
    return bam_paths


def test_basic(cmdline_args, gene_counts, exon_counts, tmpdir):
    """Tests basic invocation."""

    with mock.patch.object(expr, 'exon_counts',
                           return_value=exon_counts) as mock_ec:

        main(cmdline_args)

        # Check exon_counts call.
        bam_path = Path(str(tmpdir / 's1' / 'alignment.bam'))

        mock_ec.assert_called_once_with(
            [bam_path],
            Path('/path/to/flattened_gtf'),
            extra_kws={'-T': 1, '-s': 0},
            names={str(bam_path): 's1'})

        # Check ec results.
        ec_path = Path(cmdline_args[2]) / 'exon_counts.txt'
        assert ec_path.exists()

        ec_counts = pd.read_csv(str(ec_path), sep='\t',
                                index_col=[0, 1, 2, 3, 4])
        assert len(ec_counts) == 2
        assert ec_counts.index.names == ['gene_id', 'chr', 'start',
                                         'end', 'strand']
        assert list(ec_counts.columns) == ['s1', 's2']


def test_count_options(cmdline_args, gene_counts, exon_counts):
    """Tests invocation with extra options."""

    cmdline_args += ['--paired', '--stranded', '1', '--threads', '10']

    with mock.patch.object(expr, 'exon_counts',
                           return_value=exon_counts) as mock_ec:

        main(cmdline_args)

        # Check extra_kws.
        _, kwargs = pytest.helpers.mock_call(mock_ec)
        extra_kws = kwargs['extra_kws']

        assert extra_kws['-T'] == 10
        assert extra_kws['-p']
        assert extra_kws['-s'] == 1


def test_extra_kwargs():
    # TODO: write test.
    pass
