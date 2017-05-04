"""Module containing functions for calling multiple external tools."""

from future.utils import native_str
import pysam

from .sambamba import sambamba_sort
from .util import which


def sort_bam(input_bam, output_bam, threads=1):
    """Sorts bam file using sambamba or samtools."""

    if which('sambamba') is not None:
        sambamba_sort(input_bam, output_bam, threads=threads)
    else:
        pysam.sort(native_str(input_bam), native_str(output_bam))
