# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import argparse
import logging
import pkg_resources
import sys

from imfusion.main import build, insertions, expression, merge, ctg


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='im-fusion')

    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    # Setup sub-command parsers.
    build.setup_parser(subparsers)
    insertions.setup_parser(subparsers)
    expression.setup_parser(subparsers)
    merge.setup_parser(subparsers)
    ctg.setup_parser(subparsers)

    # Parse args.
    args = parser.parse_args(argv)

    # Setup logging.
    _setup_logging(args)

    # Log version.
    logger = logging.getLogger()

    version = pkg_resources.get_distribution('im-fusion').version
    logger.info('--- IM-Fusion v%s ---', version)

    # Dispatch.
    args.main(args)

    logger.info('--- Done! ---')


def _setup_logging(args):
    root = logging.getLogger('')

    # Setup logging format.
    formatter = logging.Formatter(
        '%(asctime)-15s %(levelname)-12s %(message)s',
        datefmt='[%Y-%m-%d %H:%M:%S]')

    # Add console logger.
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    console.setLevel(logging.INFO)

    root.addHandler(console)

    root.setLevel(logging.INFO)


if __name__ == '__main__':
    main()
