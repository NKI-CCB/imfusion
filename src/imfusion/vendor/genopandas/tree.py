# -*- coding: utf-8 -*-
"""Tree-related functions/classes."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import collections
import itertools
import operator

# pylint: disable=unused-import
from intervaltree import IntervalTree, Interval

# pylint: enable=unused-import


class GenomicIntervalTree(collections.MutableMapping):
    """Datastructure for efficiently accessing genomic objects by position."""

    def __init__(self, *args, **kwargs):
        self._trees = dict()
        self.update(dict(*args, **kwargs))

    def __getitem__(self, key):
        return self._trees[key]

    def __setitem__(self, key, value):
        self._trees[key] = value

    def __delitem__(self, key):
        del self._trees[key]

    def __iter__(self):
        return iter(self._trees)

    def __len__(self):
        return len(self._trees)

    def __repr__(self):
        args = ', '.join('{!r}: {}'.format(k, v) for k, v in self.items())
        return 'GenomicIntervalTree(**{' + args + '})'

    # pylint: disable=too-many-arguments
    def search(self,
               chromosome,
               start,
               end=None,
               strict_left=False,
               strict_right=False):
        """Searches the tree for objects within given range."""

        overlap = self[chromosome].search(start, end)

        if strict_left:
            overlap = (int_ for int_ in overlap if int_.begin >= start)

        if strict_right:
            overlap = (int_ for int_ in overlap if int_.end <= end)

        return set(overlap)

    # pylint: enable=too-many-arguments

    def is_empty(self):
        """Returns True if tree is empty."""

        if len(self) == 0:
            return True

        return all(tree.isempty() for tree in self._trees.values())

    @classmethod
    def from_tuples(cls, tuples):
        """Builds an instance from tuples.

        Assumes tuples are sorted by chromosome.
        """

        # Group by chromosome.
        tuples = sorted(tuples, key=operator.itemgetter(0))
        groups = itertools.groupby(tuples, key=operator.itemgetter(0))

        # Build trees.
        trees = {}
        for chrom, group in groups:
            trees[chrom] = IntervalTree.from_tuples(tup[1:] for tup in group)

        return cls(trees)

    def intersection(self, other):
        """Returns a new tree of all intervals common to both self and other."""

        shared = set(self.keys()) & set(other.keys())
        return self.__class__(
            {k: self[k].intersection(other[k])
             for k in shared})

    def union(self, other):
        """Returns a new tree, comprising all intervals from self and other."""

        merged = {chrom: tree.copy() for chrom, tree in self._trees.items()}

        for chrom, tree in other.items():
            if chrom in merged:
                merged[chrom] = merged[chrom].union(tree)
            else:
                merged[chrom] = tree.copy()

        return self.__class__(merged)

    def difference(self, other):
        """Returns a new tree, comprising all intervals in self but not
           in other."""

        diff = {}

        for chrom, tree in self.items():
            if chrom in other:
                diff[chrom] = tree.difference(other[chrom])
            else:
                diff[chrom] = tree.copy()

        return self.__class__(diff)

    def __or__(self, other):
        return self.union(other)

    def __and__(self, other):
        return self.intersection(other)

    def __sub__(self, other):
        return self.difference(other)
