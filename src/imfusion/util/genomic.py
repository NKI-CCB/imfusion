# -*- coding: utf-8 -*-
"""Provides functionality for efficiently accessing structured genomic data."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from heapq import heappush, heappop
import itertools
import operator
from typing import Iterable

from intervaltree import IntervalTree


class GenomicIntervalTree(object):
    """Datastructure for efficiently accessing genomic objects by position."""

    def __init__(self, trees):
        # type: (Dict[str, IntervalTree]) -> None
        self._trees = trees

    def __getitem__(self, i):
        # type: (str) -> IntervalTree
        """Returns tree with given chromosome name."""
        return self._trees[i]

    @classmethod
    def from_objects(cls, objects, chrom_attr='chromosome'):
        # type: (Iterable[object], str) -> GenomicIntervalTree
        """Builds an instance using an iterable of objects with
        chromosome, start/end attributes."""

        # Group objects by chromosome.
        keyfunc = operator.attrgetter(chrom_attr)
        objects = sorted(objects, key=keyfunc)
        grouped = itertools.groupby(objects, key=keyfunc)

        # Build trees.
        trees = {}
        for chrom, grp in grouped:
            trees[chrom] = IntervalTree.from_tuples(
                (getattr(obj, 'start'), getattr(obj, 'end'), obj)
                for obj in grp)

        return cls(trees)

    @classmethod
    def from_objects_position(cls, objects, chrom_attr='chromosome'):
        # type: (Iterable[object], str) -> GenomicIntervalTree
        """Builds an instance using an iterable of objects with
        chromosome and position attributes."""

        # Group objects by chromosome.
        keyfunc = operator.attrgetter(chrom_attr)
        objects = sorted(objects, key=keyfunc)
        grouped = itertools.groupby(objects, key=keyfunc)

        # Build trees.
        trees = {}
        for chrom, grp in grouped:
            trees[chrom] = IntervalTree.from_tuples(
                (getattr(obj, 'position'), getattr(obj, 'position') + 1, obj)
                for obj in grp)

        return cls(trees)

    def search(self, chromosome, begin, end=None):
        # type: (str, int, int) -> Iterable[object]
        """Searches the tree for objects within given range."""
        return self._trees[chromosome].search(begin, end)
