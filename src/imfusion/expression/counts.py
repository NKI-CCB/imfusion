"""Module containing functions for generating expression counts."""

# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import shutil
import subprocess
import tempfile

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pandas as pd

from imfusion.util.shell import format_kwargs


def gene_counts(bam_files, gff_path, names=None, extra_kws=None, **kwargs):
    """Generates gene counts for given bam files using featureCounts.

    This function is used to generate a m-by-n matrix (m = number of samples,
    n = number of genes) of gene expression counts. This matrix is generated
    using featureCounts, whose results are then parsed and returned.

    Parameters
    ----------
    bam_files : list[pathlib.Path]
        List of paths to the bam files for which counts should be generated.
    gff_path : pathlib.Path
        Path to the gene feature file containing gene features.
    names : dict[str, str]
        Alternative names to use for the given bam
        files. Keys of the dict should correspond to bam file paths, values
        should reflect the sample names that should be used in the
        resulting count matrix.
    extra_kws : dict[str, any]:
        Dictionary of extra arguments that
        should be passed to feature counts. Keys should correspond
        to argument names (including dashes), values should correspond
        to the argument value. Arguments without values (flags)
        should be given with the boolean value True.
    **kwargs
        Any kwargs are passed to `feature_counts`.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing counts. The index of the DataFrame contains gene
        ids corresponding to genes in the gff file, the columns correspond to
        samples/bam files. Column names are either the bam file paths, or the
        alternative sample names if given.

    """

    extra_kws = extra_kws or {}

    # Check keyword arguments.
    if '-f' in extra_kws:
        raise ValueError('Option -f can\'t be used when '
                         'counting at gene level')

    # Run feature counts.
    counts = feature_counts(
        bam_files, gff_path, names=names,
        extra_kws=extra_kws, **kwargs)

    # Drop extra columns.
    counts.drop(['Chr', 'Start', 'End',
                 'Strand', 'Length'], axis=1, inplace=True)

    # Set and rename index.
    counts.set_index('Geneid', inplace=True)
    counts.index.name = 'gene_id'

    return counts


def exon_counts(bam_files, gff_path, names=None, extra_kws=None, **kwargs):
    """Generates exon counts for given bam files using featureCounts.

    This function is used to generate a m-by-n matrix (m = number of samples,
    n = number of exons) of exon expression counts. This matrix is generated
    using featureCounts, whose results are then parsed and returned.

    Parameters
    ----------
    bam_files : list[pathlib.Path]
        List of paths to the bam files for which counts should be generated.
    gff_path : pathlib.Path
        Path to the gene feature file containing gene features.
    names : dict[str, str]
        Alternative names to use for the given bam
        files. Keys of the dict should correspond to bam file paths, values
        should reflect the sample names that should be used in the
        resulting count matrix.
    extra_kws : dict[str, any]
        Dictionary of extra arguments that
        should be passed to feature counts. Keys should correspond
        to argument names (including dashes), values should correspond
        to the argument value. Arguments without values (flags)
        should be given with the boolean value True.
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

    extra_kws = extra_kws or {}

    # Set exon-level parameters for feature counts.
    extra_kws = dict(extra_kws)
    extra_kws.update({
        '-f': True,             # Feature level
        '-t': 'exonic_part',    # Use 'exonic_part' features.
        '--minOverlap': '1',    # Minimum overlap with exon.
        '-O': True              # Include if spanning 1+ exons.
    })

    # Run feature counts.
    counts = feature_counts(
        bam_files, gff_path, names=names,
        extra_kws=extra_kws, **kwargs)

    # Drop/rename columns.
    counts.drop('Length', axis=1, inplace=True)
    counts.rename(columns={
        'Geneid': 'gene_id',
        'Chr': 'chr',
        'Start': 'start',
        'End': 'end',
        'Strand': 'strand'
    }, inplace=True)

    # Set and sort by index.
    counts.set_index(['gene_id', 'chr', 'start',
                      'end', 'strand'], inplace=True)
    counts.sort_index(inplace=True)

    return counts


def feature_counts(bam_files, gff_path, names=None,
                   extra_kws=None, tmp_dir=None, keep_tmp=False):
    """Runs featureCounts on bam files with given options.

    Main function used to run featureCounts. Used by `gene_counts` and
    `exon_counts` to generate expression counts.

    Parameters
    ----------
    bam_files : list[pathlib.Path]
        List of paths to the bam files for which counts should be generated.
    gff_path : pathlib.Path
        Path to the gff file containing gene features.
    names : dict[str, str]
        Dictionary with sample names, used to rename
        columns from file paths to sample names. Keys of the Dictionary
        should correspond with the bam file paths, values should
        reflect the desired sample name for the respective bam file.
    extra_kws : dict[str, any]
        Dictionary containing extra command line
        arguments that should be passed to featureCounts.
    tmp_dir : pathlib.Path
        Temp directory to use for the generated counts.
    keep_tmp : bool
        Whether to keep the temp directory (default = False).

    Returns
    -------
    pandas.Dataframe
        DataFrame containing feature counts for the given bam files. The rows
        correspond to the counted features, the columns correspond to the
        index values (chomosome, position etc.) and the bam files.

    """

    extra_kws = extra_kws or {}

    # Create tmpdir if needed.
    if tmp_dir is None:
        tmp_dir = Path(tempfile.mkdtemp())
    elif not tmp_dir.exists():
        tmp_dir.mkdir(parents=True)

    try:
        # Output path.
        file_path = tmp_dir / 'counts.txt'

        # Run feature counts.
        args = (['featureCounts'] + list(format_kwargs(extra_kws)) +
                ['-a', str(gff_path), '-o', str(file_path)] +
                [str(bf) for bf in bam_files])

        subprocess.check_output(args)

        # Read counts.
        counts = read_feature_counts(file_path, names=names)
    finally:
        if not keep_tmp:
            shutil.rmtree(str(tmp_dir))

    return counts


def read_feature_counts(file_path, names=None):
    """Reads a feature counts result file, renaming samples if needed.

    Parameters
    ----------
    file_path : pathlib.Path
        Path to the featureCounts output file.
    names : dict[str, str]
        Mapping of file paths to sample names that should be used
        to rename the sample columns.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing feature counts The rows correspond to the
        counted features, the columns correspond to the index values
        (chomosome, position etc.) and the samples.

    """
    # Read resulting counts.
    counts = pd.read_csv(str(file_path), sep='\t',
                         dtype={'Chr': str}, comment='#')

    # If names are given, rename columns.
    if names is not None:
        for name in names:
            if name not in counts.columns:
                print('Warning: missing sample {} for renaming'
                      .format(name))
        counts = counts.rename(columns=names)

    return counts
