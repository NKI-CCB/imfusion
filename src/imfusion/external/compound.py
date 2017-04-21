"""Module containing functions for calling multiple external tools."""

import pysam

from .sambamba import sambamba_sort
from .util import which


def sort_bam(input_bam, output_bam, threads=1):
    """Sorts bam file using sambamba or samtools."""

    if which('sambamba') is not None:
        sambamba_sort(input_bam, output_bam, threads=threads)
    else:
        pysam.sort(str(input_bam), str(output_bam))
