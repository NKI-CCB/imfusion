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
from sortedcontainers import SortedList


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


def flatten_intervals(intervals):
    """Flattens a list of (start, end) intervals."""

    intervals = sorted(intervals)

    if len(intervals) > 0:
        # Start first run.
        positions = SortedList(intervals[0])

        for interval in intervals:
            if positions[-1] < interval[0]:
                # We have left the current run, yield everything in positions.
                for new_int in _pairwise(_uniq(positions)):
                    yield new_int

                #Start new run.
                positions = SortedList(interval)
            else:
                # Record current end as part of the current run.
                positions.update(interval)

        # Yield any remaining intervals.
        for new_int in _pairwise(_uniq(positions)):
            yield new_int


def _uniq(iterable):
    """Returns uniq elements for a sorted iterable.

    Example: [1, 1, 2, 2, 2, 3, 3, 1] -> [1, 2, 3, 1

    """

    for k, _ in itertools.groupby(iterable):
        yield k


def _pairwise(iterable):
    """Returns pairs of current, next values from list.

    Example: s -> (s0,s1), (s1,s2), (s2, s3), ...

    """

    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)
