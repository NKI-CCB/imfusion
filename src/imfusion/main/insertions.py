# -*- coding: utf-8 -*-
"""Script for identifying insertions using available aligners."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import argparse
import logging

import imfusion
from imfusion.insertions.aligners import AlignerCommand
from imfusion.model import Insertion

FORMAT = "[%(asctime)-15s] %(message)s"
logging.basicConfig(
    format=FORMAT, level=logging.INFO, datefmt="%Y-%m-%d %H:%M:%S")


def main():
    """Main function for imfusion-insertions."""

    args = parse_args()
    args.command.run(args)


def parse_args():
    """Parses arguments for imfusion-insertions."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='imfusion-insertions')
    subparsers = parser.add_subparsers(dest='aligner')
    subparsers.required = True

    # Register pipelines.
    commands = AlignerCommand.available_commands()

    for name, command in commands.items():
        cmd_parser = subparsers.add_parser(name)
        command.configure(cmd_parser)
        cmd_parser.set_defaults(command=command)

    return parser.parse_args()
