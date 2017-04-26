"""Module containing functions for calling Tophat."""

from pathlib2 import Path

import tempfile
import shutil

from .util import flatten_arguments, run_command


def tophat2_index_transcriptome(bowtie_index_path,
                                gtf_path,
                                output_base_path,
                                log_path=None,
                                tmp_dir=None):
    """
    Builds a transcriptome index for a reference genome.

    Parameters
    ----------
    bowtie_index_path : pathlib.Path
        Path to the bowtie index.
    gtf_path : pathlib.Path
        Path to the reference gtf file.
    out_base_path pathlib.Path:
        Base output path for the built index.
    log_path : pathlib.Path
        Where to write the log output.
    tmp_dir : pathlib.Path
        Directory to use for temporary outputs.

    """

    if tmp_dir is None:
        tmp_dir = Path(tempfile.mkdtemp())
    else:
        tmp_dir.mkdir(parents=True, exist_ok=True)

    try:
        args = [
            'tophat2', '--GTF', str(gtf_path),
            '--transcriptome-index={}'.format(output_base_path), '--bowtie1',
            '--output-dir', str(tmp_dir), str(bowtie_index_path)
        ]

        run_command(args=args, log_path=log_path)
    finally:
        shutil.rmtree(str(tmp_dir))


def tophat2_align(fastq_path,
                  index_path,
                  output_dir,
                  fastq2_path=None,
                  extra_args=None,
                  log_path=None):
    """Aligns fastq files to a reference genome using Tophat2.

    This function is used to call TopHat2 from Python to perform an
    RNA-seq alignment. As Tophat2 is written in Python 2.7, this function
    cannot be used in Python 3.0+.

    Parameters
    ----------
    fastq_path : pathlib.Path
        Paths to the fastq file that should be used for the Tophat2
        alignment.
    index_path : pathlib.Path
        Path to the bowtie index of the (augmented)
        genome that should be used in the alignment. This index is
        typically generated by the *build_reference* function.
    output_dir : pathlib.Path
        Path to the output directory.
    fastq2_path : pathlib.Path
        Path to the fastq file of the second pair (for paired-end sequencing).
    extra_args : dict
        Dict of extra command line arguments for Tophat2.

    """

    extra_args = extra_args or {}

    # Create output_dir if needed.
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # Inject own arguments.
    extra_args['--output-dir'] = (str(output_dir), )

    # Build command-line arguments.
    if fastq2_path is None:
        fastqs = [str(fastq_path)]
    else:
        fastqs = [str(fastq_path), str(fastq2_path)]

    optional_args = list(flatten_arguments(extra_args))
    positional_args = [str(index_path)] + fastqs

    args = ['tophat2'] + optional_args + positional_args
    args = [str(s) for s in args]

    # Run Tophat2!
    run_command(args=args, log_path=log_path)
