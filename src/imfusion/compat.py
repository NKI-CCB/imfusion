# -*- coding: utf-8 -*-
"""Provides several shims for Python 2/3 compatibility."""

import os
import subprocess

# Define FileExistsError for Python 2.7.
try:
    FileExistsError = FileExistsError
except NameError:
    FileExistsError = OSError

# Define FileNotFoundError for Python 2.7.
try:
    FileNotFoundError = FileNotFoundError
except NameError:
    FileNotFoundError = OSError

# Define DEVNULL for both Python 2 and 3.
try:
    DEVNULL = subprocess.DEVNULL  # Python 3
except AttributeError:
    DEVNULL = open(os.devnull, 'wb')  # Python 2
