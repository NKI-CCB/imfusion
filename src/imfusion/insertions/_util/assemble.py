# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

from imfusion.util.shell import flatten_arguments, run_command

# def cufflinks_assemble(bam_path, gtf_path, output_dir, extra_args=None):
#     extra_args = extra_args or {}

#     # Inject args and format.
#     extra_args['--GTF-guide'] = (str(gtf_path), )
#     extra_args['--output-dir'] = (str(output_dir), )

#     extra_cmdline_args = flatten_arguments(extra_args)
#     extra_cmdline_args = [str(arg) for arg in extra_cmdline_args]

#     # Assemble full argument list.
#     cmdline_args = ['cufflinks', str(bam_path)] + extra_cmdline_args
#     run_command(cmdline_args)

#     return output_dir / 'transcripts.gtf'


def stringtie_assemble(bam_path, gtf_path, output_path, extra_args=None):
    extra_args = extra_args or {}

    # Inject args.
    extra_args['-G'] = (str(gtf_path), )
    extra_args['-o'] = (str(output_path), )

    # Format args.
    extra_cmdline_args = flatten_arguments(extra_args)
    extra_cmdline_args = [str(arg) for arg in extra_cmdline_args]

    # Assemble full argument list.
    cmdline_args = ['stringtie', str(bam_path)] + extra_cmdline_args
    run_command(args=cmdline_args)

    return output_path
