# -*- coding: utf-8 -*-
"""Contains functions for generating and reading expression counts."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

import itertools
import shutil
import subprocess
import tempfile
from typing import Any, Iterable

import pathlib2 as pathlib

import numpy as np
import pandas as pd
import toolz

from imfusion.external.feature_counts import feature_counts

# Disable E1101 checks which stumble on pandas classes.
# pylint: disable=E1101


def generate_exon_counts(
        bam_files,  # type: Iterable[pathlib.Path]
        gtf_path,  # type: pathlib.Path
        names=None,  # type: List[str]
        extra_kws=None,  # type: Dict[str, Iterable[Any]]
        tmp_dir=None,  # type: pathlib.Path
        keep_tmp=False,  # type: bool
):  # type: (...) -> pd.DataFrame
    """Generates exon counts for given bam files using featureCounts.

    This function is used to generate a m-by-n matrix (m = number of samples,
    n = number of exons) of exon expression counts. This matrix is generated
    using featureCounts, whose results are then parsed and returned.

    Parameters
    ----------
    bam_files : list[pathlib.Path]
        List of paths to the bam files for which counts should be generated.
    gtf_path : pathlib.Path
        Path to the GTF file containing gene features.
    names : dict[str, str]
        Alternative names to use for the given bam
        files. Keys of the dict should correspond to bam file paths, values
        should reflect the sample names that should be used in the
        resulting count matrix.
    extra_kws : dict[str, tuple]:
        Dictionary of extra arguments that should be passed to feature counts.
        Keys should correspond to argument names (including dashes),
        values should be tuples containing the argument values.
    tmp_dir : pathlib.Path
        Temp directory to use for generating counts.
    keep_tmp : bool
        Whether to keep the temp directory (default = False).
    **kwargs
        Any kwargs are passed to `feature_counts`.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing counts. The index of the DataFrame contains
        gene ids corresponding to exons in the gff file, the columns
        correspond to samples/bam files. Column names are either the bam
        file paths, or the alternative sample names if given.

    """

    # Set exon-level parameters for feature counts.
    default_kws = {
        '-f': True,  # Feature level
        '-t': 'exonic_part',  # Use 'exonic_part' features.
        '--minOverlap': '1',  # Minimum overlap with exon.
        '-O': True  # Include if spanning 1+ exons.
    }
    extra_kws = toolz.merge(default_kws, extra_kws or {})

    # Create tmpdir if needed.
    if tmp_dir is None:
        tmp_dir = pathlib.Path(tempfile.mkdtemp())
    elif not tmp_dir.exists():
        tmp_dir.mkdir(parents=True)

    # Run feature counts and read output.
    try:
        output_path = tmp_dir / 'counts.txt'
        feature_counts(
            bam_files=bam_files,
            gtf_path=gtf_path,
            output_path=output_path,
            extra_kws=extra_kws)
        counts = _read_feature_count_output(output_path, names=names)
    finally:
        if not keep_tmp:
            shutil.rmtree(str(tmp_dir))

    # Drop/rename columns.
    counts.drop('Length', axis=1, inplace=True)
    counts.rename(
        columns={
            'Geneid': 'gene_id',
            'Chr': 'chr',
            'Start': 'start',
            'End': 'end',
            'Strand': 'strand'
        },
        inplace=True)

    # Set and sort by index.
    counts.set_index(
        ['gene_id', 'chr', 'start', 'end', 'strand'], inplace=True)
    counts.sort_index(inplace=True)

    return counts


def _read_feature_count_output(file_path, names=None):
    # type (pathlib.Path, Dict[str, str]) -> pd.DataFrame
    """Reads counts from featureCounts output.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to count file.
    names : Optional[Dict[str, str]]
        Optional dictionary that maps featureCount column names, which are
        typically file paths, to more readable sample names.

    """

    # Read counts.
    counts = pd.read_csv(
        str(file_path), sep='\t', comment='#', dtype={'Chr': 'str'})

    # If names are given, rename columns.
    if names is not None:
        for name in names:
            if name not in counts.columns:
                # TODO: Use logging.
                print('Warning: missing sample {} for renaming'.format(name))
        counts = counts.rename(columns=names)

    return counts


def read_exon_counts(file_path, gene_id=None):
    # type: (pathlib.Path, str) -> pd.DataFrame
    """Reads exon counts from file, optionally filtering for given prefix.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to the exon count file.
    gene_id : Optional[str]
        ID of the gene to filter for.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing exon counts. The rows correspond to the
        counted features, the columns correspond to the index values
        (chomosome, position etc.) and the samples.

    """

    return _read_counts(
        file_path,
        prefix=gene_id,
        index_col=['gene_id', 'chr', 'start', 'end', 'strand'])


def _read_counts(file_path, prefix=None, **kwargs):
    # type: (pathlib.Path, str, **Any) -> pd.DataFrame
    """Reads counts from file, optionally filtering for given prefix."""

    # Read counts, optionally filtering for given prefix.
    default_kws = dict(sep='\t', dtype={'chr': 'str'}, comment='#')
    read_csv_kws = toolz.merge(default_kws, kwargs)

    if prefix is None:
        counts = pd.read_csv(str(file_path),
                             **read_csv_kws)  # type: pd.DataFrame
    else:
        counts = _read_csv_startswith(file_path, prefix=prefix, **read_csv_kws)

    return counts


def _read_csv_startswith(file_path, prefix, **kwargs):
    # type: (pathlib.Path, str, **Any) -> pd.DataFrame
    """Reads sorted csv file using only lines that start with given prefix."""

    with open(str(file_path), 'rt') as file_:
        lines = iter(file_)

        # Grab our header and column description.
        header = next(file_)

        sep = kwargs.get('sep', ',')
        columns = [c.strip() for c in header.split(sep)]

        # Drop everything until we find prefix.
        lines = itertools.dropwhile(lambda l: not l.startswith(prefix), lines)

        # Take everything while it still has prefix.
        lines = itertools.takewhile(lambda l: l.startswith(prefix), lines)

        # Read grabbed lines using pandas.
        lines_file = StringIO(u''.join(list(lines)))
        data = pd.read_csv(lines_file, names=columns, **kwargs)

    return data


def normalize_counts(counts):
    # type: (pd.DataFrame) -> pd.DataFrame
    """Normalizes counts for sequencing depth using the median-of-ratios."""

    with np.errstate(divide='ignore'):
        size_factors = estimate_size_factors(counts)
        return counts / size_factors


def estimate_size_factors(counts):
    # type: (pd.DataFrame) -> np.Array
    """Calculates size factors using the median-of-ratios approach."""

    # Convert to float.
    counts = counts.astype(float)

    # Calculate size factors.
    log_geo_means = np.mean(np.log(counts), axis=1)
    size_factors = np.apply_along_axis(
        _estimate_size_factors_col,
        axis=0,
        arr=counts,
        log_geo_means=log_geo_means)

    return size_factors


def _estimate_size_factors_col(counts, log_geo_means):
    # type: (pd.DataFrame, np.Array) -> np.Array
    log_counts = np.log(counts)
    mask = np.isfinite(log_geo_means) & (counts > 0)
    return np.exp(np.median((log_counts - log_geo_means)[mask]))


# def gene_counts(bam_files, gff_path, names=None, extra_kws=None, **kwargs):
#     """Generates gene counts for given bam files using featureCounts.

#     This function is used to generate a m-by-n matrix (m = number of samples,
#     n = number of genes) of gene expression counts. This matrix is generated
#     using featureCounts, whose results are then parsed and returned.

#     Parameters
#     ----------
#     bam_files : list[pathlib.Path]
#         List of paths to the bam files for which counts should be generated.
#     gff_path : pathlib.Path
#         Path to the gene feature file containing gene features.
#     names : dict[str, str]
#         Alternative names to use for the given bam
#         files. Keys of the dict should correspond to bam file paths, values
#         should reflect the sample names that should be used in the
#         resulting count matrix.
#     extra_kws : dict[str, tuple]:
#         Dictionary of extra arguments that should be passed to feature counts.
#         Keys should correspond to argument names (including dashes),
#         values should be tuples containing the argument values.
#     **kwargs
#         Any kwargs are passed to `feature_counts`.

#     Returns
#     -------
#     pandas.DataFrame
#         DataFrame containing counts. The index of the DataFrame contains gene
#         ids corresponding to genes in the gff file, the columns correspond to
#         samples/bam files. Column names are either the bam file paths, or the
#         alternative sample names if given.

#     """

#     extra_kws = extra_kws or {}

#     # Check keyword arguments.
#     if '-f' in extra_kws:
#         raise ValueError('Option -f can\'t be used when '
#                          'counting at gene level')

#     # Run feature counts.
#     counts = feature_counts(
#         bam_files, gff_path, names=names,
#         extra_kws=extra_kws, **kwargs)

#     # Drop extra columns.
#     counts.drop(['Chr', 'Start', 'End',
#                  'Strand', 'Length'], axis=1, inplace=True)

#     # Set and rename index.
#     counts.set_index('Geneid', inplace=True)
#     counts.index.name = 'gene_id'

#     return counts
