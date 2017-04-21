import pytest

from pathlib2 import Path


@pytest.fixture
def data_path():
    """Path to build test data directory."""
    return Path(__file__).parent.parent


@pytest.fixture
def refseq_path(data_path):
    """Path to test reference sequence."""
    return pytest.helpers.data_path('reference.fa', relative_to=data_path)


@pytest.fixture
def gtf_path(data_path):
    """Path to test reference gtf."""
    return pytest.helpers.data_path('reference.gtf', relative_to=data_path)


@pytest.fixture
def transposon_path(data_path):
    """Path to test transposon sequence."""
    return pytest.helpers.data_path('transposon.fa', relative_to=data_path)


@pytest.fixture
def features_path(data_path):
    """Path to test transposon features."""
    return pytest.helpers.data_path('features.txt', relative_to=data_path)


@pytest.fixture
def build_kws(refseq_path, gtf_path, transposon_path, features_path, tmpdir):
    """Test kws for indexer.build."""
    return {
        'refseq_path': refseq_path,
        'gtf_path': gtf_path,
        'transposon_path': transposon_path,
        'transposon_features_path': features_path,
        'output_dir': Path(str(tmpdir / 'ref'))
    }


@pytest.fixture
def cmdline_args():
    """Example command line args."""
    return [
        '--reference_seq', '/path/to/ref',
        '--reference_gtf', '/path/to/gtf',
        '--transposon_seq', '/path/to/tr',
        '--transposon_features', '/path/to/feat',
        '--output_dir', '/path/to/out'
    ] # yapf:disable
