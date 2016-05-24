import itertools

try:
    from io import StringIO
except ImportError:
    from StringIO import StringIO

import pandas as pd


def read_exon_counts(file_path, gene_id=None):
    """Reads exon counts, optionally only reading a given gene."""

    # Arguments to pass to pandas.read_csv.
    pandas_kws = {'dtype': {'chr': str}}

    if gene_id is None:
        # Read counts using pandas.
        counts = pd.read_csv(str(file_path), sep='\t', **pandas_kws)
    else:
        # Short circuit reading for specific gene.
        counts = read_csv_startswith(file_path, sep='\t',
                                     prefix=gene_id, **pandas_kws)

    # Set index columns.
    counts.set_index(['gene_id', 'chr', 'start',
                      'end', 'strand'], inplace=True)

    return counts


def read_csv_startswith(file_path, prefix, **kwargs):
    """Reads sorted csv file using only lines that start with given prefix."""

    with open(str(file_path), 'rt') as file_:
        lines = iter(file_)

        # Grab our header and column description.
        header = next(file_)

        sep = kwargs.get('sep', ',')
        columns = [c.strip() for c in header.split(sep)]

        # Drop everything until we find prefix.
        lines = itertools.dropwhile(
            lambda l: not l.startswith(prefix), lines)

        # Take everything while it still has prefix.
        lines = itertools.takewhile(
            lambda l: l.startswith(prefix), lines)

        # Read grabbed lines using pandas.
        lines_file = StringIO(u''.join(list(lines)))
        df = pd.read_csv(lines_file, names=columns, **kwargs)

    return df
