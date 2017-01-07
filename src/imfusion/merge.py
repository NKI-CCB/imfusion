# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import itertools

import pandas as pd
import toolz

from imfusion.expression import read_exon_counts
from imfusion.model import Insertion


def merge_samples(dir_paths, sample_names=None, with_expression=True):
    """Merges samples in dir_paths to a single insertions/exon counts frame.

    Parameters
    ----------
    dir_paths : List[pathlib.Path]
        Paths to the sample directories.
    samples : List[str]
        Names to use the merged samples.
    with_expression : bool
        Whether to also merge expression counts (assumes that expression.txt
        files have been generated using imfusion-expression).

    Returns
    -------
    tuple(pandas.DataFrame, pandas.DataFrame)
        Two DataFrames respectively containing the merged insertions
        and the merged exon counts. If with_expression is False, the
        merged counts frame is returned as None.

    """

    # Use directory names as sample names if none given.
    if sample_names is None:
        sample_names = [dp.name for dp in dir_paths]

    # Merge insertions.
    ins_paths = [dp / 'insertions.txt' for dp in dir_paths]
    ins_groups = [Insertion.from_csv(fp, sep='\t') for fp in ins_paths]
    merged_insertions = merge_insertions(ins_groups, sample_names)

    # Merge counts.
    if with_expression:
        count_paths = [dp / 'expression.txt' for dp in dir_paths]
        count_frames = [read_exon_counts(cp) for cp in count_paths]
        merged_counts = merge_exon_counts(count_frames, sample_names)
    else:
        merged_counts = None

    return merged_insertions, merged_counts


def merge_insertions(insertion_groups, sample_names):
    """Merges insertions from different samples."""

    # Check for duplicate names.
    if len(sample_names) != len(set(sample_names)):
        raise ValueError('Sample names contains duplicates')

    # Merge renamed insertions into single generator.
    sample_ins = (_apply_sample_name(grp, sample)
                  for (grp, sample) in zip(insertion_groups, sample_names))
    merged = itertools.chain.from_iterable(sample_ins)

    # Yield result.
    for insertion in merged:
        yield insertion


def _apply_sample_name(insertions, sample_name):
    for ins in insertions:
        new_metadata = toolz.merge(ins.metadata, {'sample': sample_name})
        yield ins._replace(id='{}.{}'.format(sample_name, ins.id),
                           metadata=new_metadata) # yapf: disable


def merge_exon_counts(count_frames, sample_names):
    """Merges expression count frames into a single frame."""

    # Check for duplicate names.
    if len(sample_names) != len(set(sample_names)):
        raise ValueError('Sample names contains duplicates')

    # Merge counts.
    merged = pd.concat(count_frames, axis=1)
    merged.columns = sample_names

    return merged
