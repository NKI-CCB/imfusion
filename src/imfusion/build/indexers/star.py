# -*- coding: utf-8 -*-
"""Implements indexer for building STAR references."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from typing import Any

import toolz

from imfusion.util import shell

from .base import Indexer, Reference, register_indexer


class StarIndexer(Indexer):
    """Indexer that builds indexes for the STAR aligner."""

    def __init__(self, logger=None, overhang=100):
        # type: (Any, int) -> None
        super().__init__(logger=logger)
        self.overhang = overhang

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

        # Create index directory.
        index_dir = reference.index_path
        index_dir.mkdir(exist_ok=False)

        # Run STAR.
        cmdline_args = [
            'STAR', '--runMode', 'genomeGenerate', '--genomeDir',
            str(index_dir), '--genomeFastaFiles', str(reference.fasta_path),
            '--sjdbGTFfile', str(reference.gtf_path), '--sjdbOverhang',
            str(self.overhang)
        ]

        log_path = reference.base_path / 'star.log'
        with log_path.open('w') as log_file:
            shell.run_command(cmdline_args, stdout=log_file)

    @classmethod
    def configure_args(cls, parser):
        super().configure_args(parser)
        star_group = parser.add_argument_group('STAR arguments')
        star_group.add_argument('--overhang', type=int, default=100)

    @classmethod
    def parse_args(cls, args):
        super_args = super().parse_args(args)
        return toolz.merge(super_args, {'overhang': args.overhang})


class StarReference(Reference):
    """Star Reference class."""
    pass


register_indexer('star', StarIndexer)
