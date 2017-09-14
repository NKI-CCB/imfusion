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
from imfusion.expression import read_exon_counts
from imfusion.model import Insertion


def main():
    """Main function for imfusion-merge."""

    args = parse_args()

    # Use directory names to name samples if no names given.
    ins_frames = [pd.read_csv(str(fp), sep='\t') for fp in args.insertions]

    # Check for overlapping samples between sets.
    seen = set()
    for ins_frame in ins_frames:
        if len(set(ins_frame['sample']) & seen) > 0:
            raise ValueError('Overlapping samples between inputs')
        seen |= set(ins_frame['sample'])

    # Merge insertions
    ins_merged = pd.concat(ins_frames, axis=0)
    ins_merged.to_csv(str(args.output), sep='\t', index=True)

    if args.expression and args.output_expression:
        if len(args.expression) != len(ins_frames):
            raise ValueError('Differing number of insertion/expression paths')

        # Read expression.
        expr_frames = [read_exon_counts(fp) for fp in args.expression]

        # Check insertion/expression frames have same samples.
        for ins_frame, expr_frame in zip(ins_frames, expr_frames):
            if not set(ins_frame['sample']) == set(expr_frame.columns):
                raise ValueError('Insertion and expression inputs do '
                                 'not have matching samples')

        # Merge expression.
        expr_merged = pd.concat(expr_frames, axis=1)
        expr_merged.to_csv(str(args.output_expression), sep='\t', index=True)


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
