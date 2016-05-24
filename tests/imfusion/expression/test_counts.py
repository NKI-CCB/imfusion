# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import mock
import pytest

from imfusion.expression import counts

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def fc_kws(tmpdir):
    return {
        'bam_files': [Path('a.bam'), Path('b.bam')],
        'gff_path': Path('/path/to/gff'),
        'names': {'a.bam': 'a', 'b.bam': 'b'},
        'tmp_dir': Path(str(tmpdir / 'counts'))
    }


@pytest.fixture
def gene_counts_path():
    return pytest.helpers.data_path('gene_counts.txt', relative_to=__file__)


@pytest.fixture
def gene_counts(gene_counts_path, fc_kws):
    return counts.read_feature_counts(gene_counts_path, names=fc_kws['names'])


class TestFeatureCounts(object):

    def test_example(self, fc_kws, gene_counts):
        with mock.patch('subprocess.check_output') as call:
            with mock.patch.object(counts, 'read_feature_counts',
                                   return_value=gene_counts):
                fcounts = counts.feature_counts(**fc_kws)

        # Check args.
        expected_args = (['featureCounts',
                          '-a', str(fc_kws['gff_path']),
                          '-o', str(fc_kws['tmp_dir'] / 'counts.txt')] +
                         list(map(str, fc_kws['bam_files'])))
        call.assert_called_once_with(expected_args)

        # Check shape of results.
        assert fcounts.shape[0] > 0
        assert list(fcounts.columns[-2:]) == ['a', 'b']

        # Check if tmpdir was removed.
        assert not fc_kws['tmp_dir'].exists()

    def test_read_counts(self, gene_counts_path):
        fcounts = counts.read_feature_counts(gene_counts_path)

        # Check shape of results.
        assert fcounts.shape[0] > 0
        assert list(fcounts.columns[-2:]) == ['a.bam', 'b.bam']

    def test_read_counts_rename(self, gene_counts_path, fc_kws):
        fcounts = counts.read_feature_counts(
            gene_counts_path, names=fc_kws['names'])

        # Check shape of results.
        assert fcounts.shape[0] > 0
        assert list(fcounts.columns[-2:]) == ['a', 'b']


class TestGeneCounts(object):

    def test_example(self, fc_kws, gene_counts):
        with mock.patch.object(counts, 'feature_counts',
                               return_value=gene_counts) as mock_fc:
            result = counts.gene_counts(**fc_kws)

            mock_fc.assert_called_once_with(
                fc_kws['bam_files'],
                fc_kws['gff_path'],
                names=fc_kws['names'],
                extra_kws={},
                tmp_dir=fc_kws['tmp_dir'])

            # Check result shape.
            assert result.shape == (2, 2)
            assert list(result.columns) == ['a', 'b']
            assert result.index.name == 'gene_id'


@pytest.fixture
def exon_counts_path():
    return pytest.helpers.data_path('exon_counts.txt', relative_to=__file__)


@pytest.fixture
def exon_counts(exon_counts_path, fc_kws):
    return counts.read_feature_counts(exon_counts_path, names=fc_kws['names'])


class TestExonCounts(object):

    def test_example(self, fc_kws, exon_counts):
        with mock.patch.object(counts, 'feature_counts',
                               return_value=exon_counts) as mock_fc:
            result = counts.exon_counts(**fc_kws)

            mock_fc.assert_called_once_with(
                fc_kws['bam_files'],
                fc_kws['gff_path'],
                names=fc_kws['names'],
                tmp_dir=fc_kws['tmp_dir'],
                extra_kws={
                    '-t': 'exonic_part',
                    '-O': True,
                    '-f': True,
                    '--minOverlap': '1'})

            # Check result shape.
            assert result.shape == (2, 2)
            assert list(result.columns) == ['a', 'b']
            assert result.index.names == ['gene_id', 'chr', 'start',
                                          'end', 'strand']
