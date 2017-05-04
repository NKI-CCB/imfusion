# -*- coding: utf-8 -*-
"""Implements indexer for building STAR references."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from typing import Any

import toolz

from imfusion.external.star import star_index
from .base import Indexer, Reference, register_indexer


class StarIndexer(Indexer):
    """Indexer that builds references for the STAR aligner.

    Performs the same steps as the base ``Indexer`` class, but additionally
    generates an index for alignment with STAR using STAR's ``genomeGenerate``
    command. Special attention should be paid to the ``overhang`` parameter,
    which defines the overhang used by STAR in the build reference
    (see the ``sjdbOverhang`` parameter in the STAR documentation for more
    details). Ideally, the value for this parameter should be one less than the
    length of the used reads.
    """

    def __init__(self, logger=None, skip_index=False, overhang=100, threads=1):
        # type: (Any, int) -> None
        super().__init__(logger=logger, skip_index=skip_index)
        self._overhang = overhang
        self._threads = threads

    @property
    def _reference_class(self):
        """Reference class to use for this indexer."""
        return StarReference

    @property
    def dependencies(self):
        """External dependencies required by this indexer."""
        return ['STAR']

    def _build_indices(self, reference):
        # type: (StarReference) -> None

        self._logger.info('Building STAR index')

        star_index(
            fasta_path=reference.fasta_path,
            gtf_path=reference.gtf_path,
            output_dir=reference.index_path,
            overhang=self._overhang,
            threads=self._threads,
            log_path=reference.base_path / 'star.log')

    @classmethod
    def configure_args(cls, parser):
        """Configures an argument parser for the Indexer.

        Used by ``imfusion-build`` to configure the sub-command for
        this indexer (if registered as an Indexer using the
        ``register_indexer`` function).

        Parameters
        ----------
        parser : argparse.ArgumentParser
            Argument parser to configure.

        """
        super().configure_args(parser)
        star_group = parser.add_argument_group('STAR arguments')
        star_group.add_argument('--star_overhang', type=int, default=100)
        star_group.add_argument('--star_threads', type=int, default=1)

    @classmethod
    def _parse_args(cls, args):
        super_args = super()._parse_args(args)
        return toolz.merge(super_args, {
            'overhang': args.star_overhang,
            'threads': args.star_threads
        })


class StarReference(Reference):
    """Star Reference class.

    Defines paths to files within the STAR reference. Currently the same
    as the base ``Reference`` class.
    """
    pass


register_indexer('star', StarIndexer)
