# -*- coding: utf-8 -*-
"""Script for identifying insertions using available aligners."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse
import logging

from pathlib2 import Path

import imfusion
from imfusion.insertions.aligners import get_aligners
from imfusion.model import Insertion

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    """Main function for imfusion-insertions."""

    args = parse_args()

    # Construct aligner and identify insertions.
    aligner = args.aligner.from_args(args)
    aligner.check_dependencies()

    insertions = aligner.identify_insertions(
        fastq_path=args.fastq,
        output_dir=args.output_dir,
        fastq2_path=args.fastq2)

    # Convert to dataframe.
    insertion_frame = Insertion.to_frame(insertions)
    insertion_frame = insertion_frame.sort_values('support', ascending=False)

    # Write output.
    output_path = args.output_dir / 'insertions.txt'
    insertion_frame.to_csv(str(output_path), sep='\t', index=False)


def parse_args():
    """Parses arguments for imfusion-insertions."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='imfusion-insertions')
    parser.add_argument(
        '--version',
        action='version',
        version='IM-Fusion ' + imfusion.__version__)

    subparsers = parser.add_subparsers(dest='aligner')
    subparsers.required = True

    # Register pipelines.
    for aligner_name, aligner_class in sorted(get_aligners().items()):
        aligner_parser = subparsers.add_parser(aligner_name)
        aligner_class.configure_args(aligner_parser)
        aligner_parser.set_defaults(aligner=aligner_class)

    return parser.parse_args()
