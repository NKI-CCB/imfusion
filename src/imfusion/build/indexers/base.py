# -*- coding: utf-8 -*-
"""Implements base indexer classes, used for building references."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

# pylint: disable=W0611

import argparse
import logging
import shutil
from typing import Any, Type, Tuple

import pathlib2 as pathlib

import pyfaidx

from imfusion.external.util import check_dependencies
from imfusion.util import tabix

from .. import util as build_util

_indexer_registry = {}


def register_indexer(name, indexer):
    # type: (str, Type[Indexer]) -> None
    """Registers indexer under given name."""
    _indexer_registry[name] = indexer


def get_indexers():
    # type: (...) -> Dict[str, Type[Indexer]]
    """Returns dict of available indexers."""
    return dict(_indexer_registry)


class Indexer(object):
    """Base Indexer class.

    This base class defines the interface for Indexer classes and provides
    basic functionality for building augmented reference genomes. Functionality
    for building the indices for alignment is provided by sub-classes
    (as this is aligner specific), typically by overriding the protected
    ``_build_indices`` method.

    The main interface of the class is formed by the ``check_dependencies``
    and ``build`` methods. The former is used to check whether required
    external dependencies are present in ``$PATH``, whilst the latter is
    called with the required reference files to build the actual reference.
    On completion, the ``build`` method returns a ``Reference`` object, which
    defines paths to the various reference files. The class of this object is
    defined by the ``_reference_class`` method, which may be overridden to
    return aligner-specific ``Reference`` subclasses.

    Additionally, the class provides the methods ``configure_args`` and
    ``from_args``, which are used to configure the argument parser for
    ``imfusion-build`` and to instantate instances from parsed arguments.
    Subclasses with extra options should override the protected ``_parse_args``
    method, which extracts constructor parameters from argparse arguments
    for ``from_args``, allowing for flexible interpretation of these arguments.

    Parameters
    ----------
    logger : logging.Logger
        Logger to be used for logging messages.
    skip_index : bool
        Whether to skip building the index. Mostly used for debugging purposes,
        as this typically results in an incomplete reference.

    """

    def __init__(self, logger=None, skip_index=False):
        # type: (Any) -> None

        self._logger = logger or logging.getLogger()
        self._skip_index = skip_index

    @property
    def _reference_class(self):  # type: (...) -> Type[Reference]
        """Reference class to use for this indexer."""
        return Reference

    @property
    def dependencies(self):  # type: (...) -> List[str]
        """External dependencies required by Indexer."""
        return []

    def check_dependencies(self):
        """Checks if all required external dependencies are in ``$PATH``.

        Raises a ValueError if any dependencies are missing.
        """

        check_dependencies(self.dependencies)

    def build(
            self,
            refseq_path,  # type: pathlib.Path
            gtf_path,  # type: pathlib.Path
            transposon_path,  # type: pathlib.Path
            transposon_features_path,  # type: pathlib.Path
            output_dir,  # type: pathlib.Path
            blacklist_regions=None,  # type: List[Tuple[str, int, int]]
            blacklist_genes=None  # type: List[str]
    ):  # type: (...) -> None
        """Builds an augmented reference containing the transposon sequence.

        The reference is built by first copying the original reference
        sequence into the output directory and then augmenting this reference
        with the transposon sequence. If needed, blacklisted genes/regions are
        then masked within this reference. Next, other supporting files
        (such as the reference GTF) are copied and indexed as needed. As a
        final step, aligner-specific indices are built for the augmented
        reference.

        Parameters
        ----------
        refseq_path : Path
            Path to the reference sequence (in Fasta format).
        gtf_path : Path
            Path to the reference GTF file.
        transposon_path : Path
            Path to the transposon sequence (in Fasta format).
        transposon_features_path : Path
            Path to the transposon features (tab-separated file).
        output_dir : Path
            Path to the output directory. Should not yet exist.
        blacklist_regions: List[Tuple[str, int, int]]
            List or regions that should be blacklisted, specified as a tuple
            of (chromosome, start, end).
        blacklist_genes : List[str]
            List of genes that should be blacklisted. Should correspond to
            ``gene_id`` entries in the passed GTF file.

        Returns
        -------
        Reference
            Returns a ``Reference`` object, which describes the paths to
            various files within the generated reference.

        """

        # Create output directory.
        try:
            output_dir.mkdir(parents=True, exist_ok=False)
        except FileExistsError as ex:
            ex.strerror = 'Output directory already exists'
            raise ex

        # Use dummy Reference instance for paths.
        reference = self._reference_class(output_dir)

        # Copy transposon files.
        self._logger.info('Copying transposon files')
        shutil.copy(str(transposon_path), str(reference.transposon_path))

        build_util.check_feature_file(transposon_features_path)
        shutil.copy(str(transposon_features_path),
                    str(reference.features_path)) # yapf: disable

        # Build indexed reference gtf.
        self._logger.info('Building indexed reference gtf')
        gtf_frame = tabix.read_gtf_frame(gtf_path)
        gtf_frame = tabix.sort_gtf_frame(gtf_frame)

        tabix.write_gtf_frame(gtf_frame, reference.gtf_path)
        tabix.index_gtf(
            reference.gtf_path,
            output_path=reference.indexed_gtf_path,
            sort=False)

        # Build flattened exon gtf.
        self._logger.info('Building flattened exon gtf')
        gtf_frame_flat = tabix.flatten_gtf_frame(gtf_frame)
        tabix.write_gtf_frame(
            gtf_frame_flat, file_path=reference.exon_gtf_path)

        # Build augmented reference.
        self._logger.info('Building augmented reference')

        blacklist_regions = blacklist_regions or []
        blacklist_genes = blacklist_genes or []

        blacklist = (build_util.regions_from_strings(blacklist_regions) +
                     build_util.regions_from_genes(blacklist_genes,
                                                   reference.indexed_gtf_path))

        build_util.build_reference(
            refseq_path,
            transposon_path,
            output_path=reference.fasta_path,
            blacklisted_regions=blacklist)

        # Build any required indices using files.
        if self._skip_index:
            self._logger.warning('Skipping the building of the index. '
                                 'Reference will not be usable for alignment.')
        else:
            self._build_indices(reference)

    def _build_indices(self, reference):
        raise NotImplementedError()

    @classmethod
    def configure_args(cls, parser):
        # type: (argparse.ArgumentParser) -> None
        """Configures an argument parser for the Indexer.

        Used by ``imfusion-build`` to configure the sub-command for
        this indexer (if registered as an Indexer using the
        ``register_indexer`` function).

        Parameters
        ----------
        parser : argparse.ArgumentParser
            Argument parser to configure.

        """

        # Basic arguments.
        base_group = parser.add_argument_group('Basic arguments')

        base_group.add_argument(
            '--reference_seq',
            type=pathlib.Path,
            required=True,
            help='Path to the reference sequence (in Fasta format).')

        base_group.add_argument(
            '--reference_gtf',
            type=pathlib.Path,
            required=True,
            help='Path to the reference gtf file.')

        base_group.add_argument(
            '--transposon_seq',
            type=pathlib.Path,
            required=True,
            help='Path to the transposon sequence (in Fasta format).')

        base_group.add_argument(
            '--transposon_features',
            type=pathlib.Path,
            required=True,
            help='Path to the transposon features (tsv).')

        base_group.add_argument(
            '--output_dir', type=pathlib.Path, required=True)

        # Optional blacklist arguments.
        blacklist_group = parser.add_argument_group('Blacklist arguments')

        blacklist_group.add_argument(
            '--blacklist_regions',
            nargs='+',
            default=(),
            help='Regions of the reference to blacklist. Should '
            'be specified as \'chromosome:start-end\'.')

        blacklist_group.add_argument(
            '--blacklist_genes',
            nargs='+',
            default=(),
            help='Genes to blacklist. Should correspond with '
            'the gene ids used in the reference gtf file.')

        debug_group = parser.add_argument_group('Debugging')
        debug_group.add_argument(
            '--skip_index',
            default=False,
            action='store_true',
            help=('Whether to skip the building of the genome indices. '
                  'Mainly used for debugging purposes.'))

    @classmethod
    def _parse_args(cls, args):
        # type: (argparse.Namespace) -> Dict[str, Any]
        """Parses argparse argument to a dict."""
        return {'skip_index': args.skip_index}

    @classmethod
    def from_args(cls, args):
        # type: (argparse.Namespace) -> Indexer
        """Constructs an Indexer instance from given arguments.

        Instantiates an Indexer instance from the given argparse arguments.
        Uses the ``_parse_args`` method internally, which performs the actual
        argument-to-parameter extraction. As such, the ``_parse_args`` method
        should be overridden in any subclasses with extra parameters.

        Parameters
        ----------
        args : argparse.Namespace
            Arguments to parse.
        """
        return cls(**cls._parse_args(args))


class Reference(object):
    """Base Reference class.

    Reference classes define paths to specific files within the reference. This
    includes, for example, paths to the augmented fasta file and to the r
    eference index. This base class defines paths that should be available
    for all references. Subclasses may define additional paths that are
    specific to the corresponding aligner.

    Parameters
    ----------
    reference_path : Path
        Path to the reference directory.

    """

    def __init__(self, reference_path):
        # type: (pathlib.Path) -> None
        if not reference_path.exists():
            raise ValueError('Reference path does not exist')
        self._reference = reference_path

    @property
    def base_path(self):
        # type: (...) -> pathlib.Path
        """Path to reference base directory."""
        return self._reference

    @property
    def fasta_path(self):
        # type: (...) -> pathlib.Path
        """Path to reference sequence."""
        return self._reference / 'reference.fa'

    @property
    def gtf_path(self):
        # type: (...) -> pathlib.Path
        """Path to reference gtf."""
        return self._reference / 'reference.gtf'

    @property
    def indexed_gtf_path(self):
        # type: (...) -> pathlib.Path
        """Path to reference gtf."""
        return self._reference / 'reference.gtf.gz'

    @property
    def exon_gtf_path(self):
        # type: (...) -> pathlib.Path
        """Path to exon gtf."""
        return self._reference / 'exons.gtf'

    @property
    def index_path(self):
        # type: (...) -> pathlib.Path
        """Path to index."""
        return self._reference / 'index'

    @property
    def transposon_name(self):
        # type: (...) -> str
        """Name of transposon sequence."""
        seqs = pyfaidx.Fasta(str(self.transposon_path)).keys()
        return list(seqs)[0]

    @property
    def transposon_path(self):
        # type: (...) -> pathlib.Path
        """Name of transposon sequence."""
        return self._reference / 'transposon.fa'

    @property
    def features_path(self):
        # type: (...) -> pathlib.Path
        """Path to transposon features."""
        return self._reference / 'features.txt'
