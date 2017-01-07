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

    # Prepare fastqs (use fastq2 for paired-end).
    if args.fastq2 is not None:
        fastqs = list(zip(args.fastq, args.fastq2))
    else:
        fastqs = args.fastq

    # Construct aligner.
    aligner = args.aligner.from_args(args)
    insertions = aligner.identify_insertions(
        read_paths=fastqs, output_dir=args.output_dir)

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
    aligners = get_aligners()

    for aligner_name in sorted(aligners.keys()):
        aligner_class = aligners[aligner_name]

        aligner_parser = subparsers.add_parser(aligner_name)
        aligner_class.configure_args(aligner_parser)
        aligner_parser.set_defaults(aligner=aligner_class)

    return parser.parse_args()
