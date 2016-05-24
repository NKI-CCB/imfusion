# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import logging

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pyfaidx

from imfusion.ctg import test_ctgs, test_de
from imfusion.model import Insertion
from imfusion.util.tabix import GtfFile


def setup_parser(subparsers):
    subparser = subparsers.add_parser('ctg')

    base_group = subparser.add_argument_group('Basic arguments')
    base_group.add_argument('--insertions', type=Path,
                            help='Path to the merged insertions file from '
                                 'im-fusion merge.')
    base_group.add_argument('--reference_seq', required=True,
                            type=pyfaidx.Fasta,
                            help='Path to the reference genome sequence (in '
                                 'fasta format). Can either be the augmented '
                                 'reference genome or the original reference.')
    base_group.add_argument('--reference_gtf', required=True, type=GtfFile,
                            help='Path to the reference gtf file. Typically '
                                 'the same file as used in im-fusion build.')
    base_group.add_argument('--output', type=Path,
                            help='Path for the output CTG file.')

    base_group.add_argument('--threshold', default=0.05, type=float,
                            help='Minimum corrected p-value for CTGs. '
                                 '(Default = 0.05).')
    base_group.add_argument('--pattern', default=None,
                            help='Regular expression reflecting the nucleotide '
                                 'sequence at which the use transposon '
                                 'typically integrates (if any). Used to '
                                 'correct for sequence integration biases '
                                 'along the genome. For example the pattern '
                                 '(AT|TA) is used for the T2onc2 transposon, '
                                 'which is biased towards integrations as '
                                 'TA sites.')
    base_group.add_argument('--window', default=None, type=int, nargs=2,
                            help='Window around the gene within which we test '
                                 'a given gene for enrichment in insertions.')

    # Insertion filtering options.
    ins_group = subparser.add_argument_group('Insertion selection')
    ins_group.add_argument('--chromosomes', default=None, nargs='+',
                           help='Chromosomes to consider. Used to omit '
                                'specific chromosomes from the CTG analysis.')
    ins_group.add_argument('--min_depth', type=int, default=None,
                           help='Minimum supporting number of reads for '
                                'insertions to be included in the CTG '
                                'analysis. Can be used to omit insertions '
                                'with low support for more confidence in '
                                'the analysis.')

    # DE options.
    de_group = subparser.add_argument_group('Differential expression')
    de_group.add_argument('--expression', default=None, type=Path,
                          help='Path to the merged expression file from '
                               'im-fusion merge.')
    de_group.add_argument('--exon_gtf', default=None, type=GtfFile,
                          help='Path to the exon gtf file. Typically the '
                               'same file as used in im-fusion expression.')
    de_group.add_argument('--de_threshold', default=0.05, type=float,
                          help='Minimum p-value for a CTG to be considered '
                               'as differentially expressed.')

    subparser.set_defaults(main=main)


def main(args):
    logger = logging.getLogger()

    # Read insertions and filter for depth.
    ins_frame = Insertion.read_frame(args.insertions)

    if args.min_depth is not None:
        ins_frame = ins_frame.query(
            'spanning_reads >= {0} or spanning_mates >= {0}'
            .format(args.min_depth))

    # Identify CTGs.
    logger.info('Testing for CTGs')

    if args.window is not None:
        logger.info('- Using window (%d, %d)', *args.window)

    ctgs = test_ctgs(ins_frame,
                     reference_seq=args.reference_seq,
                     reference_gtf=args.reference_gtf,
                     chromosomes=args.chromosomes,
                     pattern=args.pattern,
                     window=args.window)

    # If expression is given, test for differential expression.
    if args.expression is not None:
        logger.info('Testing for differential expression')
        ctgs = test_de(ctgs, ins_frame, args.exon_gtf,
                       args.expression, threshold=args.de_threshold)

    # Write outputs.
    logger.info('Writing outputs')
    ctgs.to_csv(str(args.output), sep='\t', index=False)
