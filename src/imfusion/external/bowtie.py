"""Module containing functions for calling bowtie."""

from .util import run_command


def bowtie_index(reference_path, output_base_path, log_path=None):
    """
    Builds a bowtie index for a reference genome.

    Parameters
    ----------
    reference_path : pathlib.Path
        Path to the reference genome.
    out_base_path : pathlib.Path
        Base output path for the built index.
    log_path : pathlib.Path
        Where to write the log output.

    """

    args = ['bowtie-build', str(reference_path), str(output_base_path)]
    run_command(args=args, log_path=log_path)
