"""Module containing functions for calling stringtie."""

from pathlib2 import Path

from .util import run_command, flatten_arguments


def stringtie_assemble(
        bam_path,  # type: pathlib2.Path
        gtf_path,  # type: pathlib2.Path
        output_path,  # type: pathlib2.Path
        extra_args=None,  # type: Dict[str, Iterable[Any]]
        log_path=None  # type: pathlib2.Path
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
    run_command(args, log_path=log_path)
