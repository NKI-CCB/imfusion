"""Dispatcher script that redirects calls to subcommands.

This script allows users to call subcommands from a single binary
(``imfusion build ...`` instead of ``imfusion-build`` by redirecting calls to
to the appropriate subcommands.

"""

import argparse
import subprocess
import sys

import imfusion

VALID_SUBCOMMANDS = {'build', 'insertions', 'expression', 'merge', 'ctg'}


def main():
    """Main function, parses the subcommand and executes the right script."""

    # Parse subcommand + extra arguments (which are 'unknown' args).
    parser = argparse.ArgumentParser()
    parser.add_argument('subcommand', choices=VALID_SUBCOMMANDS)
    parser.add_argument(
        '--version',
        action='version',
        version='IM-Fusion ' + imfusion.__version__)

    args, extra_args = parser.parse_known_args()

    # Pass arguments to subcommand.
    binary = 'imfusion-' + args.subcommand
    process = subprocess.Popen([binary] + extra_args)
    process.wait()

    # Copy processes returncode.
    sys.exit(process.returncode)


if __name__ == '__main__':
    main()
