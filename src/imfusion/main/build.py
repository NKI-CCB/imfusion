"""Integration and unit tests for the build command."""

# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from imfusion.aligners.tophat2 import build_reference
from imfusion.util.reference import blacklist_for_regions, blacklist_for_genes
from imfusion.util.tabix import GtfFile


def setup_parser(subparsers):
    subparser = subparsers.add_parser('build')

    # Required arguments.
    subparser.add_argument('--reference_seq', type=Path, required=True,
                           help='Path to the reference sequence '
                                '(in Fasta format).')
    subparser.add_argument('--reference_gtf', type=Path, required=True,
                           help='Path to the reference gtf file.')
    subparser.add_argument('--transposon_seq', type=Path, required=True,
                           help='Path to the transposon sequence '
                                '(in Fasta format).')
    subparser.add_argument('--output', type=Path, required=True,
                           help='Output path for the augmented reference.')

    # Optional blacklist arguments.
    subparser.add_argument('--blacklist_regions', nargs='+', default=(),
                           help='Regions of the reference to blacklist. Should '
                                'be specified as \'chromosome:start-end\'.')
    subparser.add_argument('--blacklist_genes', nargs='+', default=(),
                           help='Genes to blacklist. Should be specified as '
                                'gene IDs (Ensembl IDs for Ensembl gtfs).')

    # Optional index arguments.
    subparser.add_argument('--no_index', dest='create_index',
                           default=True, action='store_false',
                           help='Suppresses building of the bowtie index.')
    subparser.add_argument('--no_transcriptome_index', dest='create_tr_index',
                           default=True, action='store_false',
                           help='Suppresses building of the '
                                'transcriptome index.')

    # Optional misc arguments.
    subparser.add_argument('--force_overwrite', default=False,
                           action='store_true',
                           help='Overwrite any existing files.')

    subparser.set_defaults(main=main)


def main(args):
    # Create output directory if needed.
    output_dir = args.output.parent
    if output_dir and not output_dir.exists():
        output_dir.mkdir(parents=True)

    # Check if output file already exists.
    if args.output.exists() and not args.force_overwrite:
        raise ValueError('Output file already exists ({})'
                         .format(args.output))

    # Get blacklist regions.
    blacklist_regions = (
        blacklist_for_regions(args.blacklist_regions) +
        blacklist_for_genes(args.blacklist_genes, GtfFile(args.reference_gtf)))

    # Build new reference with transposon sequence.
    build_reference(
        args.reference_seq,
        args.reference_gtf,
        args.transposon_seq,
        output_path=args.output,
        blacklist_regions=blacklist_regions,
        create_index=args.create_index,
        create_transcriptome_index=args.create_tr_index)
