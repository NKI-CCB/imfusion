
# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pandas as pd
import pytest

from imfusion.main import main

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def insertion_dir():
    """Path to insertion directory."""
    return pytest.helpers.data_path('insertions_counts', relative_to=__file__)

@pytest.fixture
def mapping_path():
    """Path to sample mapping."""
    return pytest.helpers.data_path('mapping.txt', relative_to=__file__)

@pytest.fixture
def output_base(tmpdir):
    """Example output path."""
    return Path(str(tmpdir / 'merged'))


@pytest.fixture
def cmdline_args(insertion_dir, output_base):
    """Basic command line args."""
    return ['merge',
            '--base_dir', str(insertion_dir),
            '--output_base', str(output_base)]


def test_basic(cmdline_args, output_base):
    """Tests basic invocation."""
    main(cmdline_args)

    # Check written insertions.
    insertions_path = output_base.with_suffix('.insertions.txt')
    insertions = pd.read_csv(str(insertions_path), sep='\t')

    assert len(insertions) == 4
    assert 's1' in set(insertions['sample_id'])
    assert 's2' in set(insertions['sample_id'])

    # Check if ids have been prefixed with sample names.
    for _, row in insertions.iterrows():
        assert row['id'].startswith(row['sample_id'])

    # Check written counts.
    counts_path = output_base.with_suffix('.exon_counts.txt')
    exon_counts = pd.read_csv(str(counts_path), sep='\t')

    assert len(exon_counts) == 4
    assert 's1' in exon_counts.columns
    assert 's2' in exon_counts.columns


def test_subsetting(cmdline_args, output_base):
    """Tests subsetting of samples."""

    cmdline_args += ['--samples', 's1']
    main(cmdline_args)

    # Check insertions.
    insertions_path = output_base.with_suffix('.insertions.txt')
    insertions = pd.read_csv(str(insertions_path), sep='\t')
    assert set(insertions['sample_id']) == {'s1'}

    # Check exon counts.
    counts_path = output_base.with_suffix('.exon_counts.txt')
    counts = pd.read_csv(str(counts_path), sep='\t')
    assert 's1' in counts.columns
    assert 's2' not in counts.columns
