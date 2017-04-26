# -*- coding: utf-8 -*-
"""Script for building references using available indexers."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse
import logging

import imfusion
from imfusion.build.indexers import get_indexers

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    """Main function for the build command."""

    args = parse_args()

    indexer = args.indexer.from_args(args)
    indexer.check_dependencies()

    indexer.build(
        refseq_path=args.reference_seq,
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

    parser.add_argument(
        '--version',
        action='version',
        version='IM-Fusion ' + imfusion.__version__)

    subparsers = parser.add_subparsers(dest='indexer')
    subparsers.required = True

    # Register pipelines.
    for indexer_name, indexer_class in sorted(get_indexers().items()):
        indexer_parser = subparsers.add_parser(indexer_name)
        indexer_class.configure_args(indexer_parser)
        indexer_parser.set_defaults(indexer=indexer_class)

    return parser.parse_args()
