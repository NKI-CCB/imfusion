# -*- coding: utf-8 -*-
"""Script for merging insertion and expression datasets."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse
from collections import Counter
from itertools import chain

import pandas as pd
from pathlib2 import Path

import imfusion
from imfusion.expression import ExonExpressionMatrix
from imfusion.insertions import InsertionSet


def main():
    """Main function for imfusion-merge."""

    args = parse_args()

    # Check if we have the same number of files for
    # insertion and expression data.
    if args.expression and args.output_expression:
        if len(args.expression) != len(args.insertions):
            raise ValueError('Differing number of insertion/expression paths')

    # Merge insertions.
    insertions = [InsertionSet.from_csv(str(fp), sep='\t')
                  for fp in args.insertions]  # yapf: disable

    merged_insertions = InsertionSet.concat(insertions)
    merged_insertions.to_csv(args.output, sep='\t')

    # Merge expression.
    if args.expression and args.output_expression:
        expr_matrices = [ExonExpressionMatrix.from_imf(fp)
                         for fp in args.expression]  # yapf: disable
        merged_expr = ExonExpressionMatrix.concat(expr_matrices, axis=1)

        merged_expr.to_csv(args.output_expression, sep='\t')


def parse_args():
    """Parses command-line arguments for imfusion-merge."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--version',
        action='version',
        version='IM-Fusion ' + imfusion.__version__)

    parser.add_argument(
        '--insertions',
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
        '--expression',
        type=Path,
        nargs='+',
        required=False,
        default=None,
        help='Path to sample directories.')

    parser.add_argument(
        '--output_expression',
        type=Path,
        default=None,
        help='Output path for merged expression file.')

    return parser.parse_args()


if __name__ == '__main__':
    main()
