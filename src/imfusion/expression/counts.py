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

from future.utils import native, native_str
import pathlib2 as pathlib

import numpy as np
import pandas as pd
import toolz

from imfusion.external.feature_counts import feature_counts

# Disable E1101 checks which stumble on pandas classes.
# pylint: disable=E1101


class _ExpressionMatrix(object):
    """Matrix containing expression values."""

    def __init__(self, values):
        self._values = values

    @property
    def values(self):
        """Matrix values."""
        return self._values

    @classmethod
    def concat(cls, matrices, axis=1, verify_integrity=True):
        """Concatenate expression matrices."""

        merged = pd.concat(
            (mat.values for mat in matrices),
            axis=axis,
            verify_integrity=verify_integrity)

        return cls(merged)

    def to_csv(self, file_path, index=True, **kwargs):
        self._values.to_csv(native_str(file_path), index=index, **kwargs)


class GeneExpressionMatrix(_ExpressionMatrix):
    """Matrix containing gene expression values."""

    def normalize(self, size_factors=None, log2=False):
        """Normalizes expression values for sequencing depth."""

        with np.errstate(divide="ignore"):
            if size_factors is None:
                size_factors = estimate_size_factors(self._values)
            normalized = self._values.divide(size_factors, axis=1)

        if log2:
            normalized = np.log2(normalized + 1)

        return self.__class__(normalized)

    @classmethod
    def from_alignments(cls,
                        file_paths,
                        gtf_path,
                        sample_names=None,
                        feature_count_kws=None):
        """Generates gene expression counts using featureCounts."""

        # Run feature counts and read output.
        tmp_dir = pathlib.Path(tempfile.mkdtemp())
        output_path = tmp_dir / 'counts.txt'

        try:
            feature_counts(
                bam_files=file_paths,
                gtf_path=gtf_path,
                output_path=output_path,
                extra_kws=feature_count_kws or {})

            counts = cls.from_subread(output_path)
        finally:
            shutil.rmtree(str(tmp_dir))

        if sample_names is not None:
            name_map = dict(zip([str(fp) for fp in file_paths], sample_names))
            counts.values.rename(columns=name_map, inplace=True)

        return counts

    @classmethod
    def from_subread(cls, file_path, **kwargs):
        """Reads expression from a subread output file."""

        values = pd.read_csv(
            native_str(file_path),
            sep='\t',
            comment='#',
            index_col=['Geneid'],
            dtype={'Chr': 'str'},
            **kwargs)

        values = values.drop(
            ['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
        values.index.name = 'gene_id'

        return cls(values)

    @classmethod
    def from_csv(cls, file_path, index_col=0, **kwargs):
        """Reads expression from a CSV file."""
        values = pd.read_csv(
            native_str(file_path), index_col=index_col, **kwargs)
        return cls(values)

    @classmethod
    def from_exon_expression(cls, exon_expr):
        """Generates gene expression values from exon expression."""
        gene_expr = exon_expr.values.groupby(level=0).sum()
        return cls(gene_expr)


class ExonExpressionMatrix(_ExpressionMatrix):
    """Matrix containing exon expression values."""

    def get_exons(self, gene_id):
        """Returns exon positions for given gene."""

        values = self.values

        try:
            values = values.loc[[gene_id]]
        except KeyError:
            raise ValueError('Invalid gene_id {!r}'.format(gene_id))

        exons = pd.DataFrame.from_records(
            native(list(values.index.get_values())),
            columns=['gene_id', 'chromosome', 'start', 'end', 'strand'])
        exons['strand'] = exons['strand'].map({'-': -1, '+': 1})

        return exons

    @classmethod
    def from_alignments(cls,
                        file_paths,
                        gtf_path,
                        sample_names=None,
                        feature_count_kws=None):
        """Generates exon expression counts using featureCounts.

        This function is used to generate a m-by-n matrix (m = number of samples,
        n = number of exons) of exon expression counts. This matrix is generated
        using featureCounts, whose results are then parsed and returned.

        Parameters
        ----------
        file_paths : list[pathlib.Path]
            List of paths to the bam files for which counts should be generated.
        gtf_path : pathlib.Path
            Path to the GTF file containing gene features.
        sample_names : List[str]
            Sample names to use for the given bam files.
        feature_count_kws : dict[str, tuple]:
            Dictionary of extra arguments that should be passed to feature
            counts. Keys should correspond to argument names (including dashes),
            values should be tuples containing the argument values.

        Returns
        -------
        ExonExpressionCounts
            Matrix containing expression counts.

        """

        # Set exon-level parameters for feature counts.
        default_fc_kws = {
            '-f': True,  # Feature level
            '-t': 'exonic_part',  # Use 'exonic_part' features.
            '--minOverlap': '1',  # Minimum overlap with exon.
            '-O': True  # Include if spanning 1+ exons.
        }
        feature_count_kws = toolz.merge(default_fc_kws, feature_count_kws
                                        or {})

        # Run feature counts and read output.
        tmp_dir = pathlib.Path(tempfile.mkdtemp())
        output_path = tmp_dir / 'counts.txt'

        try:
            feature_counts(
                bam_files=file_paths,
                gtf_path=gtf_path,
                output_path=output_path,
                extra_kws=feature_count_kws)

            counts = cls.from_subread(output_path)
        finally:
            shutil.rmtree(str(tmp_dir))

        if sample_names is not None:
            name_map = dict(zip([str(fp) for fp in file_paths], sample_names))
            counts.values.rename(columns=name_map, inplace=True)

        return counts

    @classmethod
    def from_subread(cls, file_path, **kwargs):
        """Reads expression from a subread output file."""

        values = pd.read_csv(
            native_str(file_path),
            sep='\t',
            comment='#',
            index_col=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'],
            dtype={'Geneid': 'str',
                   'Chr': 'str'},
            **kwargs)

        values.index = values.index.droplevel(5)
        values.index.names = ['gene_id', 'chr', 'start', 'end', 'strand']

        return cls(values)

    @classmethod
    def from_imf(cls, file_path, **kwargs):
        """Reads expression from IM-Fusion expression file."""

        values = pd.read_csv(
            native_str(file_path),
            sep='\t',
            index_col=['gene_id', 'chr', 'start', 'end', 'strand'],
            dtype={'chr': 'str'},
            **kwargs)

        return cls(values)


def estimate_size_factors(counts):
    # type: (pd.DataFrame) -> np.Array
    """Calculates size factors using median-of-ratios.

    Calculates normalization factors for the median-of-ratios approach for
    normalizing expression counts for differences in sequencing depths between
    samples.

    Parameters
    ----------
    counts : pd.DataFrame
        Matrix of gene counts, with genes along the rows and samples
        along the columns.

    Returns
    -------
    np.Array
        Array of normalization factors, with one entry per sample.

    """

    def _estimate_size_factors_col(counts, log_geo_means):
        # type: (pd.DataFrame, np.Array) -> np.Array
        log_counts = np.log(counts)
        mask = np.isfinite(log_geo_means) & (counts > 0)
        return np.exp(np.median((log_counts - log_geo_means)[mask]))

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
