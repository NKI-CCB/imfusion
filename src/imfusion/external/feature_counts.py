"""Module containing functions for calling featureCounts."""

from .util import run_command, flatten_arguments


def feature_counts(
        bam_files,  # type: Iterable[pathlib.Path]
        gtf_path,  # type: pathlib.Path
        output_path,  # type: pathlib.Path
        extra_kws=None,  # type: Dict[str, Any]
        log_path=None  # type: pathlib.Path
):  # type: (...) -> None
    """Generates counts using featureCounts.

    Parameters
    ----------
    bam_files : list[pathlib.Path]
        List of paths to the bam files for which counts should be generated.
    gtf_path : pathlib.Path
        Path to the GTF file containing gene features.
    output_path : pathlib.Path
        Path to output file.
    extra_kws : dict[str, tuple]:
        Dictionary of extra arguments that should be passed to feature counts.
        Keys should correspond to argument names (including dashes),
        values should be tuples containing the argument values.

    Returns
    -------
    pandas.Dataframe
        DataFrame containing feature counts for the given bam files. The rows
        correspond to the counted features, the columns correspond to the
        index values (chomosome, position etc.) and the bam files.

    """

    extra_kws = extra_kws or {}

    # Run feature counts.
    args = (['featureCounts'] + flatten_arguments(extra_kws) +
            ['-a', str(gtf_path), '-o', str(output_path)] +
            [str(bf) for bf in bam_files])

    run_command(args=args, log_path=log_path)
