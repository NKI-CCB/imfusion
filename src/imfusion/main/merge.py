# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import logging

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from imfusion.merge import merge_samples


def setup_parser(subparsers):
    subparser = subparsers.add_parser('merge')

    subparser.add_argument('--base_dir', type=Path, required=True,
                           help='Path to a base directory that contains '
                                'outputs for individual samples as ' 'sub-directories.')
    subparser.add_argument('--output_base', type=Path, required=True,
                           help='Base name of the merged output files.')
    subparser.add_argument('--samples', nargs='+', default=None,
                           help='IDs of samples to subset the output to.')
    subparser.add_argument('--no_expression', default=True,
                           dest='with_expression', action='store_false',
                           help='Don\'t include expression counts.')

    subparser.set_defaults(main=main)


def main(args):
    """Merges insertions from different samples into a single file."""

    logger = logging.getLogger()

    # Get sample paths.
    sample_dirs = [fp.parent for fp in args.base_dir.glob('**/insertions.txt')]
    logger.info('Found %s input sample(s)', len(sample_dirs))

    # Read and merge insertions.
    logger.info('Merging samples')
    merged_ins, merged_counts = merge_samples(
        sample_dirs, samples=args.samples,
        with_expression=args.with_expression)

    # Write outputs.
    logger.info('Writing outputs')
    ins_out_path = args.output_base.with_suffix('.insertions.txt')
    merged_ins.to_csv(str(ins_out_path), sep='\t', index=False)

    if args.with_expression:
        exon_out_path = args.output_base.with_suffix('.exon_counts.txt')
        merged_counts.to_csv(str(exon_out_path), sep='\t', index=True)
