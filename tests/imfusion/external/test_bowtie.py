# -*- coding: utf-8 -*-
"""Unit tests for functions calling bowtie externally."""

from pathlib2 import Path

import pytest

from imfusion.external import bowtie


class TestBowtieIndex(object):
    """Unit tests for the bowtie_index function."""

    def test_call(self, mocker):
        """Tests example call with Path paths."""

        mock_run = mocker.patch.object(bowtie, 'run_command')
        bowtie.bowtie_index(Path('reference.fa'), Path('genome'))

        mock_run.assert_called_once_with(
            args=['bowtie-build', 'reference.fa', 'genome'], log_path=None)
