"""Module containing functions for calling STAR."""

from .util import flatten_arguments, run_command


def star_index(fasta_path,
               gtf_path,
               output_dir,
               overhang=100,
               threads=1,
               log_path=None):
    """Builds a bowtie index for a reference genome.

    Parameters
    ----------
    fasta_path : pathlib.Path
        Path to the reference genome.
    gtf_path : pathlib.Path
        Path to the reference gtf file.
    output_dir pathlib.Path:
        Output directory for the built index.
    overhang : int
        Value to use for sjdbOverhang (see STAR manual for more details).
    threads : int
        Number of threads to use.
    log_path : pathlib.Path
        Where to write the log output.

    """

    output_dir.mkdir(exist_ok=False, parents=True)

    args = [
        'STAR', '--runMode', 'genomeGenerate', '--genomeDir', str(output_dir),
        '--genomeFastaFiles', str(fasta_path), '--sjdbGTFfile', str(gtf_path),
        '--sjdbOverhang', str(overhang), '--runThreadN', str(threads)
    ]

    run_command(args=args, log_path=log_path)


def star_align(
        fastq_path,
        index_path,
        output_dir,
        fastq2_path=None,
        extra_args=None,
        threads=1,
        log_path=None, ):
    """Runs STAR in alignment mode for the given fastqs."""

    extra_args = extra_args or {}

    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # Check if files are gzipped.
    if fastq_path.suffixes[-1] == '.gz':
        extra_args['--readFilesCommand'] = ('gunzip', '-c')

    if threads > 1:
        extra_args['--runThreadN'] = threads

    # Assemble arguments.
    if fastq2_path is None:
        fastq_args = [str(fastq_path)]
    else:
        fastq_args = [str(fastq_path), str(fastq2_path)]

    args = [
        'STAR', '--genomeDir', str(index_path), '--outFileNamePrefix',
        str(output_dir) + '/', '--readFilesIn'
    ] + fastq_args

    args += [str(arg) for arg in flatten_arguments(extra_args or {})]

    run_command(args=args, log_path=log_path)
