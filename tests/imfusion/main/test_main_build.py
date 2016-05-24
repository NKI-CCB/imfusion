# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import pytest

import pyfaidx

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from imfusion.main import main
from imfusion.main import build

# pylint: disable=redefined-outer-name, no-self-use,too-few-public-methods


@pytest.fixture
def build_args(reference_seq_path, transposon_seq_path, gtf_path, output_path):
    """Example build args."""
    return ['build',
            '--reference_seq', str(reference_seq_path),
            '--reference_gtf', str(gtf_path),
            '--transposon_seq', str(transposon_seq_path),
            '--output', str(output_path)]


@pytest.fixture
def build_args_no_idx(build_args):
    """Example build args, omitting index building."""
    build_args += ['--no_index', '--no_transcriptome_index']
    return build_args


@pytest.fixture
def output_path(tmpdir):
    return Path(str(tmpdir / 'augmented.fa'))


def test_basic(build_args, reference_seq_path, reference_gtf_path,
               transposon_seq_path, output_path, mocker):
    """Integration test for basic invocation."""

    # Mock various functions.
    mock_build = mocker.spy(build, 'build_reference')
    mock_idx = mocker.patch('imfusion.aligners.tophat2.'
                            'reference.build_bowtie_index')
    mock_tr_idx = mocker.patch('imfusion.aligners.tophat2.'
                               'reference.build_transcriptome_index')

    # Run command.
    main(build_args)

    # Check created reference.
    assert output_path.exists()
    aug_ref = pyfaidx.Fasta(str(output_path))
    assert list(aug_ref.keys()) == ['1', '2', 'T2onc']

    # Check call to build_reference.
    mock_build.assert_called_once_with(
        reference_seq_path,
        reference_gtf_path,
        transposon_seq_path,
        output_path=output_path,
        blacklist_regions=[],
        create_index=True,
        create_transcriptome_index=True)

    # Check build_bowtie_index.
    mock_idx.assert_called_once_with(output_path)

    # Check build_transcriptome_index.
    mock_tr_idx.assert_called_once_with(
        output_path.with_suffix(''), reference_gtf_path,
        output_base_path=output_path.with_suffix('.transcriptome'))


def test_blacklist_regions(build_args_no_idx, reference_seq_path, output_path):
    """Integration test with blacklist regions."""

    # Add regions.
    build_args_no_idx += ['--blacklist_regions', '1:10-20']

    # Check reference is not yet modified.
    ref_seq = pyfaidx.Fasta(str(reference_seq_path))
    assert str(ref_seq['1'][10:20]) != 'N' * 10

    # Run command.
    main(build_args_no_idx)

    # Check augmented reference is modified.
    aug_seq = pyfaidx.Fasta(str(output_path))
    assert str(aug_seq['1'][10:20]) == 'N' * 10


def test_blacklist_genes(build_args_no_idx, mocker):
    """Unit test with blacklist genes."""

    # Add regions.
    build_args_no_idx += ['--blacklist_genes', 'ENSMUSG00000026510']

    # Mock build function, as gtf doesn't work with this reference.
    mock_build = mocker.patch.object(build, 'build_reference')

    # Run command.
    main(build_args_no_idx)

    # Check call.
    _, kwargs = pytest.helpers.mock_call(mock_build)
    assert kwargs['blacklist_regions'] == [('1', 182409172, 182462432)]


def test_no_index(build_args, mocker):
    """Unit test with no_index flag."""

    # Add flag.
    build_args += ['--no_index']

    # Mock build function, as gtf doesn
    mock_build = mocker.patch.object(build, 'build_reference')

    # Run command.
    main(build_args)

    # Check call.
    _, kwargs = pytest.helpers.mock_call(mock_build)
    assert not kwargs['create_index']


def test_existing(build_args_no_idx, output_path):
    """Unit test with existing output file."""

    # Already create output file.
    output_path.touch()

    with pytest.raises(ValueError):
        main(build_args_no_idx)


def test_existing_override(build_args_no_idx, output_path, mocker):
    """Test existing output file with override flag."""

    # Add overwrite flag.
    build_args_no_idx += ['--force_overwrite']

    # Already create output file.
    output_path.touch()

    # Mock build function
    mock_build = mocker.spy(build, 'build_reference')

    # Run.
    main(build_args_no_idx)

    # Check build is called.
    assert mock_build.call_count == 1


def test_missing_output_dir(build_args_no_idx, output_path):
    """Test existing output file with override flag."""

    # Set output path to subdir.
    new_output_path = output_path.parent / 'subdir' / 'augmented.fa'
    build_args_no_idx[-3] = str(new_output_path)

    # Run.
    main(build_args_no_idx)

    # Check result.
    assert new_output_path.exists()
