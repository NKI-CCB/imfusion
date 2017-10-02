# -*- coding: utf-8 -*-
"""Unit tests for functions calling combinations of tools externally."""

from pathlib2 import Path

import pytest

from imfusion.external import compound


class TestSortBam(object):
    """Unit tests for the sambamba_sort function."""

    def test_sambamba_call(self, mocker):
        """Tests call that uses sambamba."""

        mocker.patch.object(compound, 'which', return_value='sambamba')

        mock_sambamba = mocker.patch.object(compound, 'sambamba_sort')
        compound.sort_bam(Path('test.bam'), Path('sorted.bam'))

        mock_sambamba.assert_called_once_with(
            Path('test.bam'), Path('sorted.bam'), threads=1)

    def test_pysam_call(self, mocker):
        """Tests call that uses pysam (samtools)."""

        mocker.patch.object(compound, 'which', return_value=None)

        mock_pysam = mocker.patch.object(compound.pysam, 'sort')
        compound.sort_bam(Path('test.bam'), Path('sorted.bam'))

        mock_pysam.assert_called_once_with('test.bam', 'sorted.bam')
