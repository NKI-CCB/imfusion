# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import itertools
import re


def format_kwargs(kwargs):
    fmt_kwargs = map(lambda kv: _format_kwarg(*kv), kwargs.items())
    filtered_kwargs = filter(bool, fmt_kwargs)
    flattened_kwargs = itertools.chain.from_iterable(filtered_kwargs)
    return map(str, flattened_kwargs)


def _format_kwarg(key, value, delimiter=' '):
    fmt_kwarg = None
    if type(value) == bool:
        if value:
            fmt_kwarg = (key,)
    elif value is not None:
        if delimiter == ' ':
            fmt_kwarg = key, value
        else:
            fmt_kwarg = (delimiter.join((key, value)), )
    return fmt_kwarg


def parse_kwarg_str(kwarg_str):
    # Pad kwarg_str with space, used to avoid recognizing
    # '-' in between words as a separator.
    kwarg_str = ' ' + kwarg_str

    # Split into separate arg entries.
    arg_list = re.split(r'\s+(-{1,2})', kwarg_str)[1:]
    arg_list = (a + b for a, b in grouped(arg_list, 2))

    # Process entries into k, v entries in dict.
    kwargs = {}
    for arg in arg_list:
        kv = arg.split()

        if len(kv) == 1:
            kwargs[kv[0]] = True
        elif len(kv) == 2:
            kwargs[kv[0]] = kv[1]
        else:
            raise NotImplementedError(
                'Multiple argument values are not yet supported.')

    return kwargs


def grouped(iterable, n):
    """ s -> (s0,s1,s2,...sn-1),
            (sn,sn+1,sn+2,...s2n-1),
            (s2n,s2n+1,s2n+2,...s3n-1), ...
    """
    return zip(*[iter(iterable)]*n)
