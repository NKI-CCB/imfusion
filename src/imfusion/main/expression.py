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

from imfusion.expression.counts import generate_exon_counts

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    """Main function for imfusion-expression."""

    args = parse_args()

    # Assemble featureCounts arguments from options.
    feature_counts_kws = {'-s': args.stranded}

    if args.paired:
        feature_counts_kws['-p'] = True

    # Generate exon counts.
    bam_path = args.sample_dir / 'alignment.bam'

    exon_counts = generate_exon_counts(
        [bam_path],
        args.exon_gtf,
        names={str(bam_path): bam_path.parent.name},
        extra_kws=feature_counts_kws)

    # Write output.
    output = args.output or args.sample_dir / 'expression.txt'
    exon_counts.to_csv(str(output), sep='\t', index=True)


def parse_args():
    """Parses arguments for imfusion-expression."""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--sample_dir',
        type=Path,
        required=True,
        help=('Path to sample directory (typically the output '
              'of imfusion-insertions).'))

    # TODO: Move exon gtf to reference?
    parser.add_argument('--exon_gtf', type=Path, required=True)

    parser.add_argument('--output', type=Path, default=None)

    parser.add_argument(
        '--paired',
        default=False,
        action='store_true',
        help='Generate counts by counting fragments '
        'instead of reads (for paired-end data).')

    parser.add_argument(
        '--stranded',
        default=0,
        type=int,
        choices={0, 1, 2},
        help='Perform strand-specific read counting. '
        'Possible values: 0 (unstranded), 1 (stranded) '
        'and 2 (reversely stranded). 0 by default.')

    return parser.parse_args()
