"""Module containing functions for calling sambamba."""

from .util import run_command


def sambamba_sort(input_bam,
                  output_bam,
                  threads=1,
                  tmp_dir=None,
                  log_path=None):
    """Sorts bam file using sambamba.

    Parameters
    ----------
    input_bam : Path
        Path to input (unsorted) bam file.
    output_bam : Path
        Path for output (sorted) bam file.
    threads : int
        Number of threads to use.
    tmp_dir : Path
        Temporary directory to use.
    log_path : Path
        Path to log file.

    """

    if tmp_dir is None:
        tmp_dir = output_bam.parent / '_tmp'

    args = [
        'sambamba', 'sort', '-o', str(output_bam), '--tmpdir=' + str(tmp_dir),
        '-t', str(threads), str(input_bam)
    ]

    run_command(args=args, log_path=log_path)
