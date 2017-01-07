# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import argparse
import logging

from imfusion.build.indexers import get_indexers

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    """Main function for the build command."""

    args = parse_args()

    # Construct aligner.
    indexer = args.indexer.from_args(args)
    indexer.build_index(
        reference_path=args.reference_seq,
        gtf_path=args.reference_gtf,
        transposon_path=args.transposon_seq,
        transposon_features_path=args.transposon_features,
        blacklist_regions=args.blacklist_regions,
        blacklist_genes=args.blacklist_genes,
        output_dir=args.output_dir)


def parse_args():
    """Parses arguments for imfusion-build."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='imfusion-build')
    subparsers = parser.add_subparsers(dest='indexer')
    subparsers.required = True

    # Register pipelines.
    indexers = get_indexers()

    for indexer_name in sorted(indexers.keys()):
        indexer_class = indexers[indexer_name]

        indexer_parser = subparsers.add_parser(indexer_name)
        indexer_class.configure_args(indexer_parser)
        indexer_parser.set_defaults(indexer=indexer_class)

    return parser.parse_args()
