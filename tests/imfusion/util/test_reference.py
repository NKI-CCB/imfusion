# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import shutil

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pyfaidx
import pytest

from imfusion.util import reference as ref_util
from imfusion.util.tabix import GtfFile


@pytest.fixture
def refseq_path():
    """Fixture for refseq file path."""
    return pytest.helpers.data_path('reference.fa')


@pytest.fixture
def transposon_path():
    """Fixture for refseq file path."""
    return pytest.helpers.data_path('transposon.fa')


class TestConcenateFastas(object):

    def test_concatenation(self, refseq_path, transposon_path, tmpdir):
        # Concatenate fastas.
        output_path = tmpdir / 'result.fa'
        ref_util.concatenate_fastas([refseq_path, transposon_path],
                                    output_path=output_path)

        # Check result exists.
        assert output_path.ensure()

        # Check sequences.
        result = pyfaidx.Fasta(str(output_path))
        assert len(result.keys()) == 3
        assert 'T2onc' in result.keys()


    def test_write_fasta(self):
        """Tests writing of fasta sequence."""

        buf = StringIO()
        ref_util._write_fasta_record(fh=buf, name='TEST',
                                     sequence='ATCGATCGATCG', max_width=10)

        assert buf.getvalue() == '>TEST\nATCGATCGAT\nCG\n'



class TestMaskRegions(object):

    def test_mask(self, refseq_path, tmpdir):
        # Copy test file to temp dir.
        tmp_path = tmpdir / 'result.fa'
        shutil.copy(str(refseq_path), str(tmp_path))

        # Mask region.
        ref_util.mask_regions(tmp_path, [('1', 10, 20)])

        # Check result.
        fasta = pyfaidx.Fasta(str(tmp_path))
        assert str(fasta['1'][10:20]) == 'N' * 10


class TestBlacklistRegions(object):

    def test_example(self):
        """Tests example strings."""

        regions = ['1:20-100', '2:500-600']
        blacklist = ref_util.blacklist_for_regions(regions)

        # Check result.
        assert len(blacklist) == 2
        assert blacklist[0] == ('1', 20, 100)
        assert blacklist[1] == ('2', 500, 600)

    def test_empty(self):
        """Tests empty case."""

        blacklist = ref_util.blacklist_for_regions([])

        # Check result shape.
        assert len(blacklist) == 0


@pytest.fixture
def gtf():
    """Fixture for gtf file path."""
    return GtfFile(pytest.helpers.data_path('mm10.test.gtf.gz'))


class TestBlacklistGenes(object):

    def test_example(self, gtf):
        """Tests example genes."""

        gene_ids = ['ENSMUSG00000026510', 'ENSMUSG00000017146']
        blacklist = ref_util.blacklist_for_genes(gene_ids, gtf)

        # Check result.
        assert len(blacklist) == 2
        assert blacklist[0] == ('1', 182409172, 182462432)
        assert blacklist[1] == ('11', 101488764, 101551955)

    def test_empty(self, gtf):
        """Tests empty case."""

        blacklist = ref_util.blacklist_for_genes([], gtf)

        # Check result shape.
        assert len(blacklist) == 0
