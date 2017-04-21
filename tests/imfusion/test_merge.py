# -*- coding: utf-8 -*-
"""Tests for imfusion.merge module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import pytest

from pathlib2 import Path

from imfusion import merge

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def dir_paths():
    """Returns path to directory containing merge example data."""

    base_path = pytest.helpers.data_path('merge', relative_to=__file__)
    return [fp for fp in base_path.iterdir() if fp.is_dir()]


class TestMergeSamples(object):
    """Tests for merge_samples function."""

    def test_merging_insertions(self, dir_paths):
        """Tests merging of insertions."""
        # Merge insertions.
        insertions, _ = merge.merge_samples(dir_paths, with_expression=False)
        insertions = list(insertions)

        # Check result length.
        assert len(insertions) == 18

        # Check if sample names were used.
        assert insertions[0].metadata['sample'] == dir_paths[0].name
        assert insertions[0].id.startswith(dir_paths[0].name)

        assert insertions[9].metadata['sample'] == dir_paths[1].name
        assert insertions[9].id.startswith(dir_paths[1].name)

    def test_merge_ins_with_names(self, dir_paths):
        """Tests merging of insertions with sample names."""
        # Merge insertions.
        insertions, _ = merge.merge_samples(
            dir_paths, sample_names=['a', 'b'], with_expression=False)
        insertions = list(insertions)

        # Check result length.
        assert len(insertions) == 18

        # Check if sample names were used.
        assert insertions[0].metadata['sample'] == 'a'
        assert insertions[0].id.startswith('a')

        assert insertions[9].metadata['sample'] == 'b'
        assert insertions[9].id.startswith('b')

    def test_merge_ins_with_duplicates(self, dir_paths):
        """Tests merging of insertions with duplicate names."""

        with pytest.raises(ValueError):
            list(
                merge.merge_samples(
                    dir_paths, sample_names=['a', 'a'],
                    with_expression=False)[0]) # yapf: disable
