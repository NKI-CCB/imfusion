# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import logging
import shutil
import subprocess
import tempfile

from imfusion.util.reference import concatenate_fastas, mask_regions


def build_reference(ref_seq_path, ref_gtf_path, tr_seq_path, output_path,
                    blacklist_regions=None, create_index=True,
                    create_transcriptome_index=True):
    """Builds an augmented reference genome for Tophat2.

    This is the main function responsible for building augmented reference
    genomes, which are used as a reference for Tophat2 when identifying
    transposon insertions in RNA-seq samples from an insertional mutagenesis
    screen. An augmented reference contains an extra sequence which reflects
    the sequence of the transposon that was used in the screen.

    Optionally, specific regions of the original reference can be blacklisted if
    these regions contain sequences that are also in the transposon sequence,
    as these may be problematic in the alignment. For example, in our
    Sleeping Beauty screen we masked regions in the Foxf2 and in En2 genes,
    as the T2onc2 transposon we used contains sequences from these genes.
    Sequences in blacklisted regions are replaced by Ns, effectively removing
    the regions from the reference.

    Parameters
    ----------
    refseq_path : pathlib.Path
        Path to the fasta file containing the reference genome.
    ref_gtf_path : pathlib.Path
        Path to the gtf file containing the genes of the reference genome.
        Assumed to have the same structure as Ensembls reference GTF files.
    trseq_path : pathlib.Path
        Path to the fasta file containg the transposon sequence.
    output_path : pathlib.Path
        Output path for the augmented reference genome (which is written
        as a new fasta file).
    blacklist_regions : List[tuple(str, int, int)]
        List of regions that should be blacklisted in the augmented reference
        genome. Regions are specified as a tuple of (chromosome, start_position,
        end_position). For example: ('1', 2000, 2200).
    create_index : bool
        Whether a bowtie index should be created.
    create_transcriptome_index : bool
        Whether a Tophat2 transcriptome index should be created.

    """

    # Create the new reference sequence by concatenating sequences.
    logging.info('Adding transposon to reference')
    concatenate_fastas([ref_seq_path, tr_seq_path], output_path=output_path)

    # Blacklist regions if needed.
    if blacklist_regions is not None and len(blacklist_regions) > 0:
        logging.info('Masking blacklisted regions')
        mask_regions(output_path, blacklist_regions)

    # Create bowtie index unless specified otherwise.
    if create_index:
        logging.info('Building bowtie index')
        build_bowtie_index(output_path)

    # Create transcriptome index unless specified otherwise.
    if create_transcriptome_index:
        logging.info('Building transcriptome index')

        # Check if bowtie index was built.
        if not create_index:
            raise ValueError('Transcriptome index can\'t be built without '
                             'first creating a bowtie index')

        # Build transcriptome index.
        bowtie_idx_path = output_path.with_suffix('')
        tr_idx_path = output_path.with_suffix('.transcriptome')

        build_transcriptome_index(bowtie_idx_path, ref_gtf_path,
                                  output_base_path=tr_idx_path)


def build_bowtie_index(reference_path, output_base_path=None):
    """Builds a bowtie index for a reference genome.

    Parameters
    ----------
    reference_path : pathlib.Path
        Path to the reference genome.
    out_base_path pathlib.Path:
        Base output path for the built index.

    """

    if output_base_path is None:
        output_base_path = reference_path.with_suffix('')

    cmdline_args = ['bowtie-build', str(reference_path),
                    str(output_base_path)]

    with output_base_path.with_suffix('.log').open('w') as log_file:
        subprocess.check_call(cmdline_args, stdout=log_file)


def build_transcriptome_index(bowtie_index_path, reference_gtf_path,
                              output_base_path):
    """Builds a transcriptome index for a reference genome.

    Parameters
    ----------
    reference_path : pathlib.Path
        Path to the reference genome.
    reference_gtf_path : pathlib.Path
        Path to the reference gtf file.
    out_base_path pathlib.Path:
        Base output path for the built index.

    """

    tmp_dir = tempfile.mkdtemp()

    cmdline_args = ['tophat2', '--GTF', str(reference_gtf_path),
                    '--transcriptome-index={}'.format(output_base_path),
                    '--bowtie1', '--output-dir', str(tmp_dir),
                    str(bowtie_index_path)]

    try:
        with output_base_path.with_suffix('.log').open('w') as log_file:
            subprocess.check_call(args=cmdline_args, stderr=log_file)
    finally:
        shutil.rmtree(tmp_dir)
