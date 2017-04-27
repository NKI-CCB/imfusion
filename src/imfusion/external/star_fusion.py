"""Module containing functions for calling star-fusion."""

from .util import run_command


def star_fusion(junction_path, reference_path, output_dir=None, log_path=None):
    """Identifies endogenous fusions from an existing STAR alignment.

    Parameters
    ----------
    reference_path : pathlib.Path
        Path to the reference genome.
    out_base_path : pathlib.Path
        Base output path for the built index.
    log_path : pathlib.Path
        Where to write the log output.

    """

    args = [
        'STAR-Fusion', '--genome_lib_dir', str(reference_path), '-J',
        str(junction_path), '--output_dir', str(output_dir)
    ]
    run_command(args=args, log_path=log_path)
