import itertools
import operator

from intervaltree import IntervalTree


class GenomicIntervalTree(object):
    def __init__(self, trees):
        self._trees = trees

    def __getitem__(self, i):
        return self._trees[i]

    @classmethod
    def from_objects(cls, objects, chrom_attr='chromosome'):
        # Group objects by chromosome.
        keyfunc = operator.attrgetter(chrom_attr)
        objects = sorted(objects, key=keyfunc)
        grouped = itertools.groupby(objects, key=keyfunc)

        # Build trees.
        trees = {}
        for chrom, grp in grouped:
            trees[chrom] = IntervalTree.from_tuples((obj.start, obj.end, obj)
                                                    for obj in objects)

        return cls(trees)

    @classmethod
    def from_objects_position(cls, objects, chrom_attr='chromosome'):
        # Group objects by chromosome.
        keyfunc = operator.attrgetter(chrom_attr)
        objects = sorted(objects, key=keyfunc)
        grouped = itertools.groupby(objects, key=keyfunc)

        # Build trees.
        trees = {}
        for chrom, grp in grouped:
            trees[chrom] = IntervalTree.from_tuples(
                (obj.position, obj.position + 1, obj) for obj in objects)

        return cls(trees)

    def search(self, chromosome, begin, end=None):
        return self._trees[chromosome].search(begin, end)
