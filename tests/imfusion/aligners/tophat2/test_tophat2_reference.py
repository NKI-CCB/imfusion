# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pyfaidx
import pytest

from imfusion.aligners.tophat2 import reference as ref_th

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def ref_seq_path():
    """Fixture for refseq file path."""
    return pytest.helpers.data_path('reference.fa')


@pytest.fixture
def ref_gtf_path():
    """Fixture for ref gtf file path."""
    return pytest.helpers.data_path('mm10.test.gtf')


@pytest.fixture
def transposon_path():
    """Fixture for refseq file path."""
    return pytest.helpers.data_path('transposon.fa')


@pytest.fixture
def build_kws(ref_seq_path, ref_gtf_path, transposon_path, tmpdir):
    output_path = Path(str(tmpdir / 'result.fa'))
    return dict(ref_seq_path=ref_seq_path,
                ref_gtf_path=ref_gtf_path,
                tr_seq_path=transposon_path,
                output_path=output_path,
                create_transcriptome_index=False)


class TestBuildReference(object):
    """Test the building of an augmented reference for Tophat2.

    Note that we omit the actual building of the transcriptome index for now,
    as our subset gtf file does not seem to be handled properly by tophat.

    """

    def test_example(self, build_kws):
        """Tests building reference."""
        # Build reference.
        ref_th.build_reference(**build_kws)

        # Check sequences.
        output_path = build_kws['output_path']
        result = pyfaidx.Fasta(str(output_path))

        assert len(result.keys()) == 3
        assert 'T2onc' in result.keys()

        # Check if index was built.
        assert output_path.with_suffix('.1.ebwt').exists()

    def test_example_blacklist(self, build_kws, ref_seq_path):
        """Tests building reference with blacklisted regions."""
        # Build reference.
        ref_th.build_reference(blacklist_regions=[('1', 10, 20)], **build_kws)

        # Check sequences.
        output_path = build_kws['output_path']
        result = pyfaidx.Fasta(str(output_path))

        assert len(result.keys()) == 3
        assert str(result['1'][10:20]) == 'N' * 10

        # Check original reference.
        orig = pyfaidx.Fasta(str(ref_seq_path))
        assert str(orig['1'][10:20]) == 'TCTGGAGAGG'

    def test_no_bowtie_index(self, build_kws):
        """Tests building reference without index.."""
        # Build reference.
        ref_th.build_reference(create_index=False, **build_kws)

        # Check that index was not built.
        output_path = build_kws['output_path']
        assert output_path.exists()
        assert not output_path.with_suffix('.1.ebwt').exists()
