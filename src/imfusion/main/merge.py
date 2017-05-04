# -*- coding: utf-8 -*-
"""Script for merging insertion and expression datasets."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse

from pathlib2 import Path

import imfusion
from imfusion.merge import merge_samples
from imfusion.model import Insertion


def main():
    """Main function of im-fusion-merge."""

    args = _parse_args()

    # Use directory names to name samples if no names given.
    names = args.names or [fp.name for fp in args.sample_dirs]

    # Merge samples into single dataset.
    merged_ins, merged_expr = merge_samples(
        args.sample_dirs,
        sample_names=names,
        with_expression=args.output_expression is not None)

    # Write output(s).
    Insertion.to_csv(str(args.output), merged_ins, sep='\t', index=False)

    if args.output_expression is not None:
        merged_expr.to_csv(str(args.output_expression), sep='\t', index=True)


def _parse_args():
    """Parses command-line arguments for im-fusion-merge."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--version',
        action='version',
        version='IM-Fusion ' + imfusion.__version__)

    parser.add_argument(
        '--sample_dirs',
        type=Path,
        nargs='+',
        required=True,
        help='Path to sample directories.')

    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output path for merged insertion file.')

    parser.add_argument(
        '--names',
        nargs='+',
        required=False,
        default=None,
        help='Alternative sample names to use for samples in merged dataset.')

    parser.add_argument(
        '--output_expression',
        type=Path,
        default=None,
        help='Output path for merged expression file.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
