# -*- coding: utf-8 -*-
"""Provides functionality for interacting with shell calls and arguments."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import os
import subprocess
from typing import Any, Iterable, Optional

import pyparsing as pp


def run_command(args, stdout=None, stderr=None, logger=None, **kwargs):
    # type(List[str], Any, Any, Any, **Any) -> None
    """Runs command using subprocess.check_call."""

    if logger is not None:
        logger.info('Running command: %s', ' '.join(args))

    try:
        subprocess.check_call(args, stdout=stdout, stderr=stderr, **kwargs)
    except FileNotFoundError as ex:
        extra_msg = ('. Make sure that the required dependencies '
                     'have been installed.')
        ex.strerror += extra_msg
        raise ex


def parse_arguments(arg_str):
    # type: (List[str]) -> Dict[str, Iterable[str]]
    """Parses command line arguments from a string."""

    parser = _setup_parser()
    parsed = parser.parseString(arg_str)

    return dict(parsed.asList())


def _setup_parser():
    arg_prefix = pp.Literal('-')
    arg_prefix2 = pp.Literal('--')

    arg_name = pp.Word(pp.alphanums)

    non_shell_chars = pp.alphanums + ':/*%_'
    arg_value = pp.Word(non_shell_chars)

    argument = (
        pp.Or([arg_prefix, arg_prefix2]) + arg_name + pp.ZeroOrMore(arg_value))
    argument.setParseAction(
        lambda tokens: (tokens[0] + tokens[1], tuple(tokens[2:])))

    argument_list = pp.ZeroOrMore(argument)

    return argument_list


def flatten_arguments(arg_dict):
    # type: (Dict[str, Iterable[Any]]) -> List[str]
    """Flattens a dict of arguments into a list for subprocess."""

    args = []  # type: List[str]

    for arg_name, arg_value in sorted(arg_dict.items()):
        if isinstance(arg_value, bool):
            args += [arg_name]
        elif hasattr(arg_value, '__iter__') and not isinstance(arg_value, str):
            args += [arg_name] + list(arg_value)
        else:
            args += [arg_name, str(arg_value)]

    return args


def which(program):
    # type: (str) -> Optional[str]
    """Determine the location of given executable, if present."""

    def _is_execuable(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if _is_execuable(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if _is_execuable(exe_file):
                return exe_file

    return None


def check_dependencies(programs):
    # type: (Iterable[str]) -> None
    """Checks if listed executables are all available in $PATH."""

    missing = [prog for prog in programs if which(prog) is None]

    if len(missing) > 0:
        raise ValueError('Missing required external dependencies: {}'
                         .format(', '.join(missing)))
