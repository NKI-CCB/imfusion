# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import subprocess

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path


def get_version(path=None):
    """Get the version of Tophat2 in path.

    Parameters
    ----------
    path : pathlib.Path
        Path to use for the Tophat2 executable.

    Returns
    -------
    str
        Version of the Tophat2 executable in path.

    """

    exec_path = Path(path or '') / 'tophat2'
    output = subprocess.check_output([exec_path, '--version'])

    return output.strip().split()[-1].replace(b'v', b'')


def get_bowtie_version(path=None):
    """Get the version of Bowtie in path.

    Parameters
    ----------
    path : pathlib.Path
        Path to use for the bowtie executable.

    Returns
    -------
    str
        Version of Bowtie executable in path.

    """

    exec_path = Path(path or '') / 'bowtie'

    output = subprocess.check_output([exec_path, '--version'])
    output_lines = output.split(b'\n')

    return output_lines[0].strip().split()[-1]
