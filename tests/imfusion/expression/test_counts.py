# -*- coding: utf-8 -*-
"""Tests for imfusion.expression.counts module."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import tempfile
import shutil

from pathlib2 import Path
import pandas as pd
import pytest

from imfusion.expression import counts
from imfusion.expression.counts import ExonExpressionMatrix

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


class TestExonExpressionMatrix(object):
    """Tests for ExonExpressionMatrix."""

    def test_from_subread(self):
        """Tests reading from subread files."""

        file_path = pytest.helpers.data_path(
            'fc_counts.txt', relative_to=__file__)

        matrix = ExonExpressionMatrix.from_subread(file_path)

        assert list(matrix.values.columns) == ['a.bam', 'b.bam']
        assert list(matrix.values.index.names) == [
            'gene_id', 'chr', 'start', 'end', 'strand'
        ]

    def test_from_imf(self):
        """Tests reading from IM-Fusion expression files."""

        file_path = pytest.helpers.data_path(
            'exon_counts.txt', relative_to=__file__)

        matrix = ExonExpressionMatrix.from_imf(file_path)

        assert list(matrix.values.columns) == ['a', 'b']
        assert list(matrix.values.index.names) == [
            'gene_id', 'chr', 'start', 'end', 'strand'
        ]

    def test_get_exons(self):
        """Tests getting exons."""
        file_path = pytest.helpers.data_path(
            'exon_counts.txt', relative_to=__file__)

        matrix = ExonExpressionMatrix.from_imf(file_path)
        exons = matrix.get_exons(gene_id='gene_a')

        expected = pd.DataFrame(
            {
                'gene_id': 'gene_a',
                'chromosome': '1',
                'start': [10, 30, 50, 70, 90],
                'end': [20, 40, 60, 80, 100],
                'strand': 1
            },
            columns=['gene_id', 'chromosome', 'start', 'end', 'strand'])

        assert all(exons == expected)

    def test_get_exons_invalid(self):
        """Tests getting exons with invalid gene id."""

        file_path = pytest.helpers.data_path(
            'exon_counts.txt', relative_to=__file__)

        matrix = ExonExpressionMatrix.from_imf(file_path)

        with pytest.raises(ValueError):
            matrix.get_exons(gene_id='gene_x')

    def test_from_alignments(self, mocker, tmpdir, fc_counts):
        """Tests generating counts using featureCounts."""

        mocker.patch.object(tempfile, 'mkdtemp', return_value=str(tmpdir))
        mock_fc = mocker.patch.object(counts, 'feature_counts')

        exon_counts = ExonExpressionMatrix.from_alignments(
            file_paths=[Path('a.bam'), Path('b.bam')],
            gtf_path=Path('/path/to/reference.gtf'),
            sample_names=['a', 'b'],
            feature_count_kws={})

        mock_fc.assert_called_once_with(
            bam_files=[Path('a.bam'), Path('b.bam')],
            gtf_path=Path('/path/to/reference.gtf'),
            output_path=tmpdir / 'counts.txt',
            extra_kws={
                '--minOverlap': '1',
                '-O': True,
                '-f': True,
                '-t': 'exonic_part'
            })

        assert len(exon_counts.values) > 0
        assert list(exon_counts.values.columns) == ['a', 'b']
