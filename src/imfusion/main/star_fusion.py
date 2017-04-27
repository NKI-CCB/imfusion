# -*- coding: utf-8 -*-
"""Script for identifying endogenous fusions using STAR-Fusion."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse
import logging

from pathlib2 import Path

import imfusion
from imfusion.external.star_fusion import star_fusion
from imfusion.external.util import FileNotFoundError

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    """Main function for imfusion-expression."""

    args = parse_args()

    junction_path = args.sample_dir / '_star' / 'Chimeric.out.junction'

    if not junction_path.exists():
        raise FileNotFoundError('Junction file ({}) could not be found. '
                                'Was the sample analyzed using STAR?'.format(
                                    junction_path))

    star_fusion(
        junction_path,
        reference_path=args.reference,
        output_dir=args.output_dir)


def parse_args():
    """Parses arguments for imfusion-starfusion."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--version',
        action='version',
        version='IM-Fusion ' + imfusion.__version__)

    parser.add_argument(
        '--sample_dir',
        type=Path,
        required=True,
        help=('Path to sample directory (typically the output '
              'of imfusion-insertions).'))

    parser.add_argument(
        '--star_fusion_reference',
        type=Path,
        required=True,
        help='Path to STAR fusion reference.')

    parser.add_argument(
        '--output_dir',
        type=Path,
        default=None,
        help='Path to output directory.')

    return parser.parse_args()
