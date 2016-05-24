"""Unit tests for insertions main.

Effectively tests link between the command line module and the
identify_insertions function in the insertions module.
"""

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
import mock

from imfusion.main import main
from imfusion.model import Insertion


_patch_func = 'imfusion.main.insertions.identify_insertions'
# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods

@pytest.fixture
def feature_path():
    return pytest.helpers.data_path('transposon.features.txt')

@pytest.fixture
def gtf_path():
    return pytest.helpers.data_path('mm10.test.gtf.gz')

@pytest.fixture
def ins_args(tmpdir, feature_path, gtf_path):
    """Example build args."""
    return ['insertions',
            '--fastq', 'in1.fastq', 'in2.fastq',
            '--reference_index', '/path/to/index',
            '--reference_gtf', str(gtf_path),
            '--transposon_name', 'T2onc',
            '--transposon_features', str(feature_path),
            '--output_dir', str(tmpdir)]

@pytest.fixture
def insertions():
    return [Insertion(id='A', seqname='1', position=10, strand=1,
                      sample_id='s1', gene_id='ENSMUSG00000026510',
                      spanning_reads=10),
            Insertion(id='B', seqname='1', position=20, strand=1,
                      sample_id='s1', gene_id='ENSMUSG00000000001',
                      spanning_reads=5)]


@pytest.fixture
def output_path(tmpdir):
    return Path(str(tmpdir / 'insertions.txt'))


def test_basic(ins_args, insertions, feature_path, gtf_path, output_path):
    """Test basic invocation with required arguments."""

    with mock.patch(_patch_func, return_value=insertions) as func:
        main(ins_args)

        # Check call.
        assert len(func.call_args_list) == 1
        kwargs = func.call_args_list[0][1]

        assert kwargs['fastqs'] == [Path('in1.fastq'), Path('in2.fastq')]
        assert kwargs['index_path'] == Path('/path/to/index')
        assert kwargs['reference_gtf_path'] == Path(gtf_path)
        assert kwargs['transposon_name'] == 'T2onc'
        assert kwargs['work_dir'] == output_path.parent
        assert kwargs['tophat_kws'] is None

        # Check output file.
        assert output_path.exists()

        # Check written frame.
        df = pd.read_csv(str(output_path), sep='\t')
        assert list(df['id']) == ['A', 'B']


def test_transcriptome_index(ins_args, insertions):
    """Tests invocation with transcriptome index."""

    ins_args += ['--transcriptome_index', '/path/to/tr_index']

    with mock.patch(_patch_func, return_value=insertions) as func:
        main(ins_args)

        # Check call.
        assert len(func.call_args_list) == 1
        kwargs = func.call_args_list[0][1]

        assert kwargs['transcriptome_index'] == Path('/path/to/tr_index')



def test_filter_depth(ins_args, insertions, tmpdir, output_path):
    """Test filtering on depth."""

    ins_args += ['--min_support', '6']

    with mock.patch(_patch_func, return_value=insertions):
        main(ins_args)

        # Check output file.
        assert output_path.exists()

        # Check written frame.
        df = pd.read_csv(str(output_path), sep='\t')
        assert list(df['id']) == ['A']


def test_filter_blacklist(ins_args, insertions, tmpdir, output_path):
    """Test filtering of blacklisted genes."""

    ins_args += ['--blacklist', 'ENSMUSG00000026510']

    with mock.patch(_patch_func, return_value=insertions):
        main(ins_args)

        # Check output file.
        assert output_path.exists()

        # Check written frame.
        df = pd.read_csv(str(output_path), sep='\t')
        assert list(df['id']) == ['B']


def test_paired(ins_args, insertions):
    """Test paired-end fastqs."""

    ins_args += ['--fastq2', 'in1.R2.fastq', 'in2.R2.fastq']

    with mock.patch(_patch_func, return_value=insertions) as func:
        main(ins_args)

        # Check call.
        assert len(func.call_args_list) == 1
        kwargs = func.call_args_list[0][1]

        assert kwargs['fastqs'] == [(Path('in1.fastq'), Path('in1.R2.fastq')),
                                    (Path('in2.fastq'), Path('in2.R2.fastq'))]

def test_non_existing_output_dir(ins_args, insertions, tmpdir):
    output_idx = ins_args.index('--output_dir') + 1
    ins_args[output_idx] = str(tmpdir / 'subdir')

    with mock.patch(_patch_func, return_value=insertions) as func:
        main(ins_args)
