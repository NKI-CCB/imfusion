"""Module containing functions for calling stringtie."""

from .util import run_command, flatten_arguments


def stringtie_assemble(
        bam_path,  # type: pathlib.Path
        gtf_path,  # type: pathlib.Path
        output_path,  # type: pathlib.Path
        extra_args=None  # type: Dict[str, Iterable[Any]]
):  # type: (...) -> None
    """Runs stringtie to assemble transripts."""

    extra_args = extra_args or {}

    # Inject args.
    extra_args['-G'] = (str(gtf_path), )
    extra_args['-o'] = (str(output_path), )

    # Format args.
    extra_cmdline_args = flatten_arguments(extra_args)
    extra_cmdline_args = [str(arg) for arg in extra_cmdline_args]

    # Assemble full argument list.
    args = ['stringtie', str(bam_path)] + extra_cmdline_args
    run_command(args)
