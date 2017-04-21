# -*- coding: utf-8 -*-
"""Provides functionality for interacting with files and file paths."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import pathlib2 as pathlib


def symlink_relative(src_path, dest_path):
    # type: (pathlib.Path, pathlib.Path) -> None
    """Symlinks file using relative path."""
    dest_path.symlink_to(src_path.relative_to(dest_path.parent))
