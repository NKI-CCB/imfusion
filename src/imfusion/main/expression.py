# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import logging

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from imfusion.expression.counts import exon_counts
from imfusion.util.shell import parse_kwarg_str


def setup_parser(subparsers):
    subparser = subparsers.add_parser('expression')

    subparser.add_argument('--sample_dir', type=Path, required=True,
                           help='Path to the sample directory. Effectively '
                                'the same as the output_dir in im-fusion ' 'insertions.')

    # subparser.add_argument('--reference_gtf', type=Path, required=True)
    subparser.add_argument('--exon_gtf', type=Path, required=True,
                           help='Path to the exon gtf file, which contains '
                                'a flattened representation of the exons in '
                                'the previously used reference gtf.')

    subparser.add_argument('--sample_id', default=None,
                           help='Sample id to use in the generated counts '
                                'file. Should reflect the same sample id as '
                                'used by im-fusion insertions. Defaults to the ' 'name of the input directory.')

    subparser.add_argument('--paired', default=False, action='store_true',
                           help='Generate counts by counting fragments '
                                'instead of reads (for paired-end data).')
    subparser.add_argument('--stranded', default=0,
                           type=int, choices={0, 1, 2},
                           help='Perform strand-specific read counting. '
                                'Possible values: 0 (unstranded), 1 (stranded) '
                                'and 2 (reversely stranded). 0 by default.')

    subparser.add_argument('--threads', default=1, type=int,
                           help='Number of threads to use in featureCounts.')
    subparser.add_argument('--extra_kwargs', default='',
                           help='Extra command line options to pass '
                                'to featureCounts.')

    subparser.set_defaults(main=main)


def main(args):
    """Generates gene/exon counts for alignments from im-fusion inputs."""

    logger = logging.getLogger()

    # Check input directory.
    if not (args.sample_dir / 'insertions.txt').exists():
        raise ValueError('Given input directory does not contain output '
                         'from IM-Fusion (no insertions.txt file found)')

    # Get bam path and sample id.
    bam_path = args.sample_dir / 'alignment.bam'
    sample_id = args.sample_id or args.sample_dir.stem

    # Generate gene counts.
    extra_kws = _args_to_extra_kws(args)

    logger.info('Generating exon counts')
    counts_exon = exon_counts(
        [bam_path], args.exon_gtf,
        names={str(bam_path): sample_id}, extra_kws=extra_kws)

    # Write outputs.
    logger.info('Writing output')
    exon_path = args.sample_dir / ('exon_counts.txt')
    counts_exon.to_csv(str(exon_path), sep='\t')


def _args_to_extra_kws(args):
    extra_kws = parse_kwarg_str(args.extra_kwargs)

    if args.paired:
        extra_kws['-p'] = True

    extra_kws['-s'] = args.stranded
    extra_kws['-T'] = args.threads

    return extra_kws
