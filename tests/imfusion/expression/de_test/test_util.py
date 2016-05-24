# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import pytest

from imfusion.expression.de_test.util import read_exon_counts

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def exon_counts_path():
    return pytest.helpers.data_path('exon_counts.txt', relative_to=__file__)

def test_read(exon_counts_path):
    """Tests reading of example file."""
    counts = read_exon_counts(exon_counts_path, gene_id='ENSMUSG00000090025')

    assert counts.shape == (1, 2)
    assert list(counts.columns) == ['a', 'b']
    assert counts.index.names == ['gene_id', 'chr', 'start', 'end', 'strand']


def test_read_with_rename(exon_counts_path):
    """Tests reading of example file with renaming of columns."""
    counts = read_exon_counts(exon_counts_path, gene_id='ENSMUSG00000090025',
                              name_map={'a': 'c', 'b': 'd'})

    assert counts.shape == (1, 2)
    assert list(counts.columns) == ['c', 'd']
    assert counts.index.names == ['gene_id', 'chr', 'start', 'end', 'strand']
