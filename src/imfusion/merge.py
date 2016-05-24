# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import collections
import pandas as pd

from imfusion.expression.util import read_exon_counts
from imfusion.model import Insertion


def merge_samples(dir_paths, samples=None):
    """Merges samples in dir_paths to a single insertions/exon counts frame.

    Parameters
    ----------
    dir_paths : List[pathlib.Path]
        Paths to the sample directories.
    samples : List[str]
        Samples to subset the results to.

    Returns
    -------
    tuple(pandas.DataFrame, pandas.DataFrame)
        Two DataFrames respectively containing the merged insertions
        and the merged exon counts.

    """

    # Merge insertions.
    ins_paths = [dp / 'insertions.txt' for dp in dir_paths]
    ins_frames = [Insertion.read_frame(ip) for ip in ins_paths]
    merged_insertions = merge_insertion_frames(ins_frames)

    # Merge counts.
    count_paths = [dp / 'exon_counts.txt' for dp in dir_paths]
    count_frames = [read_exon_counts(cp) for cp in count_paths]
    merged_counts = merge_exon_counts(count_frames)

    # Subset to samples if given.
    if samples is not None:
        merged_insertions = merged_insertions.ix[
            merged_insertions['sample_id'].isin(samples)]
        merged_counts = merged_counts[samples]

    return merged_insertions, merged_counts


def merge_insertion_frames(insertion_frames):
    """Merges insertion frames into a single frame."""

    # Check for duplicate samples.
    seen_samples = set()
    for frame in insertion_frames:
        frame_samples = set(frame['sample_id'])
        overlap = frame_samples & seen_samples

        if len(overlap) > 0:
            raise ValueError('Duplicate samples ({})'
                             .format(', '.join(overlap)))

        seen_samples = seen_samples | frame_samples

    # Merge insertions.
    merged = pd.concat(insertion_frames, ignore_index=True)

    # Augment original ids with sample names to keep them unique.
    merged['id'] = (merged['sample_id'] + '.' + merged['id'])

    return merged


def merge_exon_counts(count_frames):
    """Merges expression count frames into a single frame."""

    # Merge counts.
    merged = pd.concat(count_frames, axis=1)

    # Check for duplicate samples.
    col_counts = collections.Counter(merged.columns)
    duplicates = {sample: count for sample, count in col_counts.items()
                  if count > 1}

    if len(duplicates) > 0:
        raise ValueError('Duplicate samples ({})'
                         .format(', '.join(duplicates.keys())))

    return merged
