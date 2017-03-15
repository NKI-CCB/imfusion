# -*- coding: utf-8 -*-
"""Tests for imfusion.build.util module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import shutil

import pyfaidx
import pytest

from imfusion.build import util

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def ref_path():
    """Returns path to example reference."""
    return pytest.helpers.data_path('reference.fa', relative_to=__file__)


@pytest.fixture
def gtf_path():
    """Returns path to example gtf."""
    return pytest.helpers.data_path('reference.gtf.gz', relative_to=__file__)


class TestRegionsFromGenes(object):
    """Tests for the regions_from_genes function."""

    def test_example(self, gtf_path):
        """Test example genes."""

        regions = util.regions_from_genes(
            gene_ids=['ENSMUSG00000085584'], gtf_path=gtf_path)

        assert regions == [('X', 143099593, 143104297)]

    def test_unknown_gene(self, gtf_path):
        """Tests gene not in gtf."""

        with pytest.raises(ValueError):
            util.regions_from_genes(
                gene_ids=['ENSMUSG00000095584'], gtf_path=gtf_path)


class TestMaskReference(object):
    """Tests for the mask_reference function."""

    def test_example(self, ref_path, tmpdir):
        """Tests modifying a specific region."""

        # Copy reference to avoid modifying file.
        tmp_ref_path = tmpdir / 'reference.fa'
        shutil.copy(str(ref_path), str(tmp_ref_path))

        # Check original sequence.
        refseq = pyfaidx.Fasta(str(tmp_ref_path))
        assert str(refseq['1'][10:14]) != 'NNNN'

        # Mask.
        util.mask_reference(tmp_ref_path, [('1', 10, 14)])

        # Check modified sequence.
        refseq = pyfaidx.Fasta(str(tmp_ref_path))
        assert str(refseq['1'][10:14]) == 'NNNN'
