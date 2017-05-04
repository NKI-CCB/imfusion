# -*- coding: utf-8 -*-
"""Implements indexer for building Tophat references."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import tempfile
import shutil

import pathlib2 as pathlib

from imfusion.external.bowtie import bowtie_index
from imfusion.external.tophat import tophat2_index_transcriptome

from .base import Indexer, register_indexer, Reference


class TophatIndexer(Indexer):
    """Indexer that builds references for the Tophat-Fusion aligner.

    Performs the same steps as the base ``Indexer`` class, but additionally
    generates an index for alignment with Tophat-Fusion using Bowtie and
    Tophat2.
    """

    def __init__(self, logger=None, skip_index=False):
        super().__init__(logger=logger, skip_index=False)

    @property
    def _reference_class(self):
        """Reference class to use for this indexer."""
        return TophatReference

    @property
    def dependencies(self):
        """External dependencies required by this indexer."""
        return ['bowtie-build', 'tophat2']

    def _build_indices(self, reference):
        # type: (TophatReference) -> None

        # Build bowtie index.
        self._logger.info('Building bowtie index')

        bowtie_idx_path = reference.index_path
        bowtie_log_path = reference.base_path / 'bowtie.log'
        bowtie_index(
            reference_path=reference.fasta_path,
            output_base_path=bowtie_idx_path,
            log_path=bowtie_log_path)

        # Build transcriptome index.
        self._logger.info('Building transcriptome index')
        transcriptome_path = reference.transcriptome_path
        transcriptome_log_path = reference.base_path / 'transcriptome.log'
        tophat2_index_transcriptome(
            bowtie_index_path=bowtie_idx_path,
            gtf_path=reference.gtf_path,
            output_base_path=transcriptome_path,
            log_path=transcriptome_log_path)


register_indexer('tophat', TophatIndexer)


class TophatReference(Reference):
    """Tophat Reference class.

    Defines paths to files within the Tophat-Fusion reference. Compared to the
    base reference, this class adds an additional path to the transcriptome
    index that is used by Tophat during alignment.
    """

    @property
    def transcriptome_path(self):
        # type: (...) -> pathlib.Path
        """Path to transcriptome index."""
        return self._reference / 'transcriptome'

    @property
    def index_path(self):
        # type: (...) -> pathlib.Path
        """Path to index."""
        return self._reference / 'reference'
