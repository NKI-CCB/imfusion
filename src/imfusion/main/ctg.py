# -*- coding: utf-8 -*-
"""Script for testing for CTGs."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse
import logging

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pandas as pd

from imfusion.build import Reference
from imfusion.ctg import test_ctgs
from imfusion.expression.counts import read_exon_counts
from imfusion.expression.test import test_de_genes
from imfusion.model import Insertion

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


# pylint: disable=E1101
def main():
    logger = logging.getLogger()
    args = _parser_args()

    # Read insertions and filter for depth.
    insertions = list(Insertion.from_csv(args.insertions, sep='\t'))

    if args.min_depth is not None:
        insertions = [
            ins for ins in insertions if ins.support > args.min_depth
        ]

    # Identify CTGs.
    logger.info('Testing for CTGs')

    if args.window is not None:
        logger.info('- Using window (%d, %d)', *args.window)

    reference = Reference(args.reference)

    ctgs = test_ctgs(
        insertions,
        reference=reference,
        gene_ids=args.gene_ids,
        chromosomes=args.chromosomes,
        pattern=args.pattern,
        window=args.window)

    # Filter using given threshold.
    if args.threshold is not None:
        ctgs = ctgs.query('q_value <= {}'.format(args.threshold))

    # If expression is given, test for differential expression.
    if args.expression is not None:
        logger.info('Testing for differential expression')

        # Perform DE tests.
        exon_counts = read_exon_counts(args.expression)
        de_results = test_de_genes(
            insertions, exon_counts, gene_ids=ctgs['gene_id'])

        # Combine with CTG result.
        de_results = de_results.rename(
            columns={'direction': 'de_direction',
                     'p_value': 'de_pvalue'})

        ctgs = pd.merge(ctgs, de_results, on='gene_id', how='left')

        if args.de_threshold is not None:
            ctgs = ctgs.query('de_pvalue < {}'.format(args.de_threshold))

    # Write outputs.
    logger.info('Writing outputs')
    ctgs.to_csv(str(args.output), sep='\t', index=False)


def _parser_args():
    parser = argparse.ArgumentParser()

    base_group = parser.add_argument_group('Basic arguments')

    base_group.add_argument(
        '--insertions',
        required=True,
        type=Path,
        help='Path to the merged insertions file from '
        'im-fusion merge.')

    base_group.add_argument(
        '--reference',
        required=True,
        type=Path,
        help='Path to the reference genome sequence (in '
        'fasta format). Can either be the augmented '
        'reference genome or the original reference.')

    base_group.add_argument(
        '--gene_ids', default=None, nargs='+', help='IDs of genes to test.')

    base_group.add_argument(
        '--output',
        type=Path,
        help='Path for the output CTG file.',
        required=True)

    base_group.add_argument(
        '--threshold',
        default=0.05,
        type=float,
        help='Minimum corrected p-value for CTGs. '
        '(Default = 0.01).')

    base_group.add_argument(
        '--pattern',
        default=None,
        help='Regular expression reflecting the nucleotide '
        'sequence at which the use transposon '
        'typically integrates (if any). Used to '
        'correct for sequence integration biases '
        'along the genome. For example the pattern '
        '(AT|TA) is used for the T2onc2 transposon, '
        'which is biased towards integrations as '
        'TA sites.')

    base_group.add_argument(
        '--window',
        default=None,
        type=int,
        nargs=2,
        help='Window around the gene within which we test '
        'a given gene for enrichment in insertions.')

    # Insertion filtering options.
    ins_group = parser.add_argument_group('Insertion selection')

    ins_group.add_argument(
        '--chromosomes',
        default=None,
        nargs='+',
        help='Chromosomes to consider. Used to omit '
        'specific chromosomes from the CTG analysis.')

    ins_group.add_argument(
        '--min_depth',
        type=int,
        default=None,
        help='Minimum supporting number of reads for '
        'insertions to be included in the CTG '
        'analysis. Can be used to omit insertions '
        'with low support for more confidence in '
        'the analysis.')

    # DE options.
    de_group = parser.add_argument_group('Differential expression')

    de_group.add_argument(
        '--expression',
        default=None,
        type=Path,
        help='Path to the merged expression file from '
        'im-fusion merge.')

    de_group.add_argument(
        '--de_threshold',
        default=0.05,
        type=float,
        help='Minimum p-value for a CTG to be considered '
        'as differentially expressed.')

    return parser.parse_args()


def _test_de(insertions, expression_path, dexseq_gtf, gene_ids):
    # Test each gene for differential expression.
    de_results = {}
    for gene_id in gene_ids:
        try:
            de_results[gene_id] = de_exon(insertions, expression_path,
                                          dexseq_gtf, gene_id)
        except ValueError:
            pass

    # Summarize result in a DataFrame.
    de_results = pd.DataFrame(
        ((gene_id, res.p_value, res.direction)
         for gene_id, res in de_results.items()),
        columns=['gene_id', 'de_pvalue', 'de_direction'])

    return de_results


if __name__ == '__main__':
    main()
