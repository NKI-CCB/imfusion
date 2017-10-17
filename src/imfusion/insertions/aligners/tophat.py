# -*- coding: utf-8 -*-
"""Implements aligner for building Tophat2 references."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from future.utils import native_str
import pandas as pd
from pathlib2 import Path
import pysam
import toolz

from imfusion.build.indexers.tophat import TophatReference
from imfusion.external.stringtie import stringtie_assemble
from imfusion.external.tophat import tophat2_align
from imfusion.external.util import parse_arguments
from imfusion.insertions import TransposonFusion
from imfusion.util import path, tabix
from imfusion.vendor.frozendict import frozendict

from .base import Aligner, AlignerCommand


class TophatAligner(Aligner):
    """Tophat2 aligner.

    Aligner that identifies insertions using Tophat-Fusion (as implemented in
    Tophat2), which aligns RNA-seq reads and returns a summary file detailing
    the identified fusions. Gene-transposon fusions are identified directly
    from this summary file.

    External dependencies include Tophat2 and Bowtie1, which are both used for
    alignment. Stringtie is required for the detection of novel transcripts.

    Parameters
    ----------
    reference : Reference
        Reference class describing the reference paths.
    assemble : bool
        Whether to perform the assembly of novel transcripts using Stringie.
    assemble_args : Dict[str, Any]
        Extra arguments to pass to Stringie during assembly.
    min_flank : int
        Minimum number of flanking bases required on both sides of the
        breakpoint for a fusion to be considered valid. This means that, for
        a split read overlapping the fusion boundary, at least ``min_flank``
        bases should be on either sides of the fusion.
    threads : int
        Number of threads to use.
    extra_args :
        Extra arguments to pass to Tophat2 for alignment.
    logger : logging.Logger
        Logger to be used for logging messages.
    filter_features : bool
        Whether insertions should be filtered to remove non SA/SD insertions.
    filter_orientation : bool
        Whether insertions should be filtered for fusions in which the
        transposon feature is in the wrong orientation.
    filter_blacklist : List[str]
        List of gene ids to filter insertions for.

    """

    def __init__(self,
                 reference,
                 threads=1,
                 assemble=False,
                 assemble_args=None,
                 min_flank=12,
                 extra_args=None,
                 logger=None,
                 filter_features=True,
                 filter_orientation=True,
                 filter_blacklist=None):

        super().__init__(reference=reference, logger=logger, threads=threads)

        self._assemble = assemble
        self._assemble_args = assemble_args or {}
        self._min_flank = min_flank
        self._threads = threads
        self._extra_args = extra_args or {}

        self._filter_features = filter_features
        self._filter_orientation = filter_orientation
        self._filter_blacklist = filter_blacklist

    @property
    def dependencies(self):
        programs = ['tophat2', 'bowtie']

        if self._assemble:
            programs += ['stringtie']

        return programs

    def identify_insertions(self,
                            fastq_path,
                            output_dir,
                            fastq2_path=None,
                            sample=None):
        """Identifies insertions from given reads."""

        # Perform alignment using STAR.
        alignment_path = output_dir / 'alignment.bam'
        if not alignment_path.exists():
            self._logger.info('Performing alignment using Tophat2')
            self._align(fastq_path, output_dir, fastq2_path=fastq2_path)
            pysam.index(native_str(alignment_path))
        else:
            self._logger.info('Using existing Tophat2 alignment')

        # Assemble transcripts if requested.
        if self._assemble:
            assembled_path = output_dir / 'assembled.gtf.gz'
            if not assembled_path.exists():
                self._logger.info('Assembling transcripts using Stringie')

                # Generate assembled GTF.
                stringtie_out_path = assembled_path.with_suffix('')
                stringtie_assemble(
                    alignment_path,
                    gtf_path=self._reference.gtf_path,
                    output_path=stringtie_out_path)

                # Compress and index.
                tabix.index_gtf(stringtie_out_path, output_path=assembled_path)
                stringtie_out_path.unlink()
            else:
                self._logger.info('Using existing Stringtie assembly')
        else:
            assembled_path = None

        # Extract identified fusions.
        self._logger.info('Extracting gene-transposon fusions')
        fusion_path = output_dir / 'fusions.out'
        fusions = self._extract_fusions(fusion_path, sample_name=sample)

        # Extract insertions.
        self._logger.info('Summarizing insertions')
        insertions = list(
            util.extract_insertions(
                fusions,
                gtf_path=self._reference.indexed_gtf_path,
                features_path=self._reference.features_path,
                assembled_gtf_path=assembled_path,
                ffpm_fastq_path=fastq_path,
                chromosomes=None))

        insertions = util.filter_insertions(
            insertions,
            features=self._filter_features,
            orientation=self._filter_orientation,
            blacklist=self._filter_blacklist)

        for insertion in insertions:
            yield insertion

    def _align(self, fastq_path, output_dir, fastq2_path=None):
        # Setup args.
        transcriptome_path = self._reference.transcriptome_path

        args = {
            '--transcriptome-index': (str(transcriptome_path), ),
            '--fusion-search': (),
            '--fusion-anchor-length': (self._min_flank, ),
            '--bowtie1': (),
            '--num-threads': (self._threads, )
        }

        args = toolz.merge(args, self._extra_args)

        # Run Tophat2.
        tophat_dir = output_dir / '_tophat'

        tophat2_align(
            fastq_path=fastq_path,
            fastq2_path=fastq2_path,
            index_path=self._reference.index_path,
            output_dir=tophat_dir,
            extra_args=args)

        # Symlink alignment into expected location for gene counts.
        path.symlink_relative(
            src_path=tophat_dir / 'accepted_hits.bam',
            dest_path=output_dir / 'alignment.bam')

        path.symlink_relative(
            src_path=tophat_dir / 'fusions.out',
            dest_path=output_dir / 'fusions.out')

    def _extract_fusions(self, fusion_path, sample_name):
        fusion_data = read_fusion_out(fusion_path)

        fusions = extract_transposon_fusions(
            fusion_data,
            transposon_name=self._reference.transposon_name,
            sample_name=sample_name)

        for fusion in fusions:
            yield fusion


class TophatCommand(AlignerCommand):

    name = 'tophat'

    def configure(self, parser):
        super().configure(parser)

        group = parser.add_argument_group('Tophat2 arguments')

        group.add_argument(
            '--tophat_min_flank',
            type=int,
            default=12,
            help=('Minimum mapped length of the two segments '
                  'on each side of the fusion.'))

        group.add_argument(
            '--tophat_args',
            type=parse_arguments,
            default='',
            help='Additional args to pass to Tophat2.')

        assemble_group = parser.add_argument_group('Assembly')
        assemble_group.add_argument(
            '--assemble',
            default=False,
            action='store_true',
            help='Perform de-novo transcript assembly using StringTie.')

        filt_group = parser.add_argument_group('Filtering')
        filt_group.add_argument(
            '--no_filter_orientation',
            dest='filter_orientation',
            default=True,
            action='store_false',
            help=('Don\'t filter fusions with transposon features and genes '
                  'in opposite (incompatible) orientations.'))

        filt_group.add_argument(
            '--no_filter_feature',
            dest='filter_features',
            default=True,
            action='store_false',
            help=('Don\'t filter fusions with non-SA/SD features.'))

        filt_group.add_argument(
            '--blacklisted_genes',
            nargs='+',
            help='Blacklisted genes to filter.')

    @classmethod
    def _build_aligner(cls, args):
        return TophatAligner(
            reference=TophatReference(args.reference),
            threads=args.threads,
            min_flank=args.tophat_min_flank,
            extra_args=args.tophat_args,
            assemble=args.assemble,
            filter_features=args.filter_features,
            filter_orientation=args.filter_orientation,
            filter_blacklist=args.blacklisted_genes)


def read_fusion_out(fusion_path):
    """Reads fusion.out file from Tophat2.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to the Tophat fusion file (fusions.out).

    Returns
    -------
    pandas.DataFrame
        DataFrame containing gene fusions.

    """

    if isinstance(fusion_path, Path):
        fusion_path = str(fusion_path)

    # Read fusions using pandas.
    names = [
        'chromosomes', 'position_a', 'position_b', 'orientation', 'supp_reads',
        'supp_mates', 'supp_spanning_mates', 'contradicting_reads', 'flank_a',
        'flank_b'
    ]

    fusions = pd.read_csv(
        fusion_path, sep='\t', header=None, usecols=range(0, 10), names=names)

    # Split combined entries.
    fusions = _split_fields(fusions)

    # Map orientation to strands.
    fusions['strand_a'] = fusions['strand_a'].map({'f': 1, 'r': -1})
    fusions['strand_b'] = fusions['strand_b'].map({'f': 1, 'r': -1})

    return fusions


def _split_fields(fusions):
    """Splits combined chromosomes/strand entries in tophat fusion frame."""

    columns = [
        'chromosome_a', 'position_a', 'strand_a', 'chromosome_b', 'position_b',
        'strand_b', 'supp_reads', 'supp_mates', 'supp_spanning_mates',
        'contradicting_reads', 'flank_a', 'flank_b'
    ]

    if len(fusions) == 0:
        # Return empty 'split' frame.
        fusions = pd.DataFrame.from_records([], columns=columns)
    else:
        # Split combined entries.
        fusions['chromosome_a'], fusions['chromosome_b'] = \
            zip(*fusions['chromosomes'].str.split('-'))

        fusions['strand_a'], fusions['strand_b'] = \
            zip(*fusions['orientation'].apply(list))

        # Subset/Reorder columns.
        fusions = fusions[columns]

    return fusions


def extract_transposon_fusions(fusion_data, transposon_name, sample_name):
    """
    Extracts gene-transposon fusions from a Tophat fusion.out file.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to the Tophat fusion file (fusions.out).
    transposon_name : str
        Name of the transposon sequence in the augmented reference
        genome that was used for the alignment.
    sample_name : str
        Name of corresponding sample.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing gene-transposon fusions.

    """

    is_paired = (any(fusion_data['supp_mates'] > 0)
                 or any(fusion_data['supp_spanning_mates'] > 0))

    # Select fusions where one chromosome is the transposon and the other isn't.
    fusion_data = fusion_data.loc[(
        (fusion_data['chromosome_a'] == transposon_name) ^
        (fusion_data['chromosome_b'] == transposon_name))]

    # Build frame with candidate insertions.
    for _, row in fusion_data.iterrows():
        yield _to_tr_fusion(row, transposon_name, sample_name, is_paired)


def _to_tr_fusion(row, transposon_name, sample_name, is_paired):
    if row.chromosome_a == transposon_name:
        gen_id, tr_id = 'b', 'a'
        gen_dir, tr_dir = 1, -1
    else:
        gen_id, tr_id = 'a', 'b'
        gen_dir, tr_dir = -1, 1

    strand_genome = row['strand_' + gen_id]
    strand_transposon = row['strand_' + tr_id]

    if is_paired:
        support_junction = row.supp_spanning_mates
        support_spanning = row.supp_mates
    else:
        support_junction = row.supp_reads
        support_spanning = 0

    return TransposonFusion(
        chromosome=row['chromosome_' + gen_id],
        anchor_genome=row['position_' + gen_id],
        anchor_transposon=row['position_' + tr_id],
        flank_genome=row['flank_' + gen_id] * strand_genome * gen_dir,
        flank_transposon=row['flank_' + tr_id] * strand_transposon * tr_dir,
        strand_genome=strand_genome,
        strand_transposon=strand_transposon,
        support_junction=support_junction,
        support_spanning=support_spanning,
        sample=sample_name)
