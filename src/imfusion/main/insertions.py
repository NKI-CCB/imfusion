# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import argparse
import logging

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from imfusion.insertions.aligners import get_aligners
from imfusion.model import Insertion

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    """Main function for imfusion-insertions."""

    args = parse_args()

    # Construct aligner.
    aligner = args.aligner.from_args(args)
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
    subparsers = parser.add_subparsers(dest='aligner')
    subparsers.required = True

    # Register pipelines.
    for aligner_name, aligner_class in sorted(get_aligners().items()):
        aligner_parser = subparsers.add_parser(aligner_name)
        aligner_class.configure_args(aligner_parser)
        aligner_parser.set_defaults(aligner=aligner_class)

    return parser.parse_args()
