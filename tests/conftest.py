import os

from pathlib2 import Path

pytest_plugins = ['helpers_namespace']

import pytest

BASE_DIR = Path(__file__).parent


@pytest.helpers.register
def data_path(relative_path, relative_to=None):
    """Returns data path to test file."""

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


# For some reason, mock.ANY would not work.
# Using own implementation for now.
@pytest.helpers.register
def mock_any(cls):
    """Returns Any class to use in mock."""

    class Any(cls):
        """Any class used for mocks."""

        def __eq__(self, other):
            return True

    return Any()


@pytest.helpers.register
def touch(path):
    """Touch file."""
    path = str(path)
    with open(path, 'a'):
        os.utime(path, None)
