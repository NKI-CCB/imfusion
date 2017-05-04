# -*- coding: utf-8 -*-
"""Unit tests for functions calling sambamba externally."""

from pathlib2 import Path

import pytest

from imfusion.external import sambamba


class TestSambambaSort(object):
    """Unit tests for the sambamba_sort function."""

    def test_call(self, mocker):
        """Tests example call with Path paths."""

        mock_run = mocker.patch.object(sambamba, 'run_command')
        sambamba.sambamba_sort(Path('test.bam'), Path('sorted.bam'))

        mock_run.assert_called_once_with(
            args=[
                'sambamba', 'sort', '-o', 'sorted.bam', '--tmpdir=_tmp', '-t',
                '1', 'test.bam'
            ],
            log_path=None)
