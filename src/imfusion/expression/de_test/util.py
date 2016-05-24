# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401


import itertools

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

import pandas as pd


def read_csv_startswith(file_path, prefix, **kwargs):
    """Reads sorted csv file, using only lines that
       start with the given prefix."""

    with open(str(file_path), 'rt') as file_:
        lines = iter(file_)

        # Grab our header and column description.
        header = next(file_)

        sep = kwargs.get('sep', ',')
        columns = header.split(sep)

        # Drop everything until we find prefix.
        lines = itertools.dropwhile(
            lambda l: not l.startswith(prefix), lines)

        # Take everything while it still has prefix.
        lines = itertools.takewhile(
            lambda l: l.startswith(prefix), lines)

        # Read lines using pandas.
        lines_file = StringIO(u''.join(list(lines)))
        df = pd.read_csv(lines_file, names=columns, **kwargs)

    return df


def read_exon_counts(file_path, gene_id, name_map=None):
    """Reads exon counts for the given gene, mapping sample_ids if needed."""

    # Read counts.
    counts = read_csv_startswith(file_path, sep='\t',
                                 prefix=gene_id, dtype={'chr': str})
    counts.set_index(['gene_id', 'chr', 'start',
                      'end', 'strand'], inplace=True)

    # Remove extra \n from column.
    counts.columns = [c.strip() for c in counts.columns]

    # Map names to sample_ids.
    if name_map is not None:
        sel_cols = [c for c in counts.columns if c in name_map]

        counts = counts[sel_cols]
        counts.columns = [name_map[c] for c in sel_cols]

    return counts
