pytest_plugins = ['helpers_namespace']

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pytest


BASE_DIR = Path(__file__).parent


@pytest.helpers.register
def data_path(relative_path, relative_to=None):
    if relative_to is None:
        # Use BASE_DIR as default.
        relative_to = BASE_DIR
    elif not isinstance(relative_to, Path):
        # Ensure relative_to is a Path.
        relative_to = Path(relative_to)

    # If relative_to is not a path, move up one level.
    if not relative_to.is_dir():
        relative_to = relative_to.parent

    return relative_to / 'data' / relative_path


@pytest.helpers.register
def mock_call(mock_obj, index=0):
    args, kwargs = mock_obj.call_args_list[index]
    return args, kwargs


@pytest.fixture
def reference_seq_path():
    return pytest.helpers.data_path('reference.fa')


@pytest.fixture
def transposon_seq_path():
    return pytest.helpers.data_path('transposon.fa')


@pytest.fixture
def reference_gtf_path(gtf_path):
    return gtf_path


@pytest.fixture
def gtf_path():
    return pytest.helpers.data_path('mm10.test.gtf.gz')
