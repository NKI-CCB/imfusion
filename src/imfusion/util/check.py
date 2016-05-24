# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401


def check_features(transposon_features):
    """Checks if a transposon feature frame is valid."""

    # Check dataframe columns.
    for req_col in ['name', 'start', 'end', 'strand', 'type']:
        if req_col not in transposon_features.columns:
            raise ValueError('Missing required column '
                             '{!r} in transposon features'.format(req_col))

    # Check if we have any features.
    if len(transposon_features) < 1:
        raise ValueError('Transposon features is empty')

    # Check if we have a valid feature.
    num_valid = transposon_features['type'].isin({'SA', 'SD'}).sum()
    if num_valid < 1:
        raise ValueError('No valid features (with type = SD or SA)'
                         'were found in transposon features')
