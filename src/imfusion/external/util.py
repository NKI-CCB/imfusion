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

# Define FileNotFoundError for Python 2.7.
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError

# Define DEVNULL.
try:
    DEVNULL = subprocess.DEVNULL  # Python 3
except AttributeError:
    DEVNULL = open(os.devnull, 'wb')  # Python 2


def run_command(args, log_path=None, **kwargs):
    """Runs command using subprocess.check_call."""

    try:
        process = subprocess.Popen(
            args, stdout=DEVNULL, stderr=subprocess.PIPE, **kwargs)
        _, stderr = process.communicate()

        if log_path is not None:
            with log_path.open('w') as log_file:
                print(log_file, file=log_file)

        if process.returncode != 0:
            raise CalledProcessError(
                returncode=process.returncode,
                cmd=args,
                stderr=stderr.decode())

    except FileNotFoundError as ex:
        extra_msg = ('. Make sure that the required dependencies '
                     'have been installed.')
        ex.strerror += extra_msg
        raise ex


class CalledProcessError(subprocess.CalledProcessError):
    """Custom CalledProcessError that includes stderr in message.

    Parameters
    ----------
    *args : Any
        Any args are passed to the CalledProcessError super class.
    stderr : str
        Stderror output to include in message.
    indent : str
        Indentation prefix to use for stderr message.
    **kwargs : Any
        Any other kwargs are passed to the CalledProcessError super class.

    """

    def __init__(self, *args, **kwargs):
        stderr = kwargs.pop('stderr', None)
        indent = kwargs.pop('indent', '    ')

        super().__init__(*args, **kwargs)

        self._stderr = stderr
        self._indent = indent

    def __str__(self):
        message = super().__str__()

        if self._stderr is not None:
            stderr = _indent(_tail(self._stderr), prefix=self._indent)
            message += ('\n\nCaptured the following output '
                        'from stderr:\n\n' + stderr)

        return message


def _indent(message, prefix):
    """Indents message by given prefix.

    Used instead of textwrap.indent for Python 2.7.
    """

    lines = message.split(os.linesep)
    return os.linesep.join(prefix + l for l in lines)


def _tail(message, n=25):
    """Returns tail of message."""

    lines = message.split(os.linesep)
    first = max(0, len(lines) - n)

    return os.linesep.join(lines[first:])


def parse_arguments(arg_str):
    # type: (List[str]) -> Dict[str, Iterable[str]]
    """Parses command line arguments from a string.

    Recognises arguments of the following form:
        - -f
        - --flag
        - -f <value>
        - --flag <value> ...

    """

    # Parse string.
    parser = _setup_parser()
    parsed = parser.parseString(arg_str)

    result = dict(parsed.asList())

    # Check if we were able to parse something.
    if len(arg_str.strip()) > 1 and len(result) == 0:
        raise ValueError('Unable to parse string {!r} into arguments. '
                         'Ensure string only contains properly formatted '
                         'optional arguments, e.g. \'-f\', \'-f <value>\' or '
                         '\'--flag <value>\'.'.format(arg_str))

    return result


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
