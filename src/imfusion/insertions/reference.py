# -*- coding: utf-8 -*-

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from collections import namedtuple
import itertools
import operator

from future.utils import native_str
from intervaltree import IntervalTree
import pysam


class TranscriptReference(object):
    """Reference class, used for efficiently looking up features in
    the reference transciptome."""

    def __init__(self, transcript_trees, exon_trees):
        self._transcripts = transcript_trees
        self._exons = exon_trees

    @classmethod
    def from_gtf(
            cls,
            gtf_path,  # type: pathlib.Path
            chromosomes=None,  # type: List[str]
            record_filter=None  # type: Callable[[Any], bool]
    ):  # type: (...) -> TranscriptReference
        """Builds an Reference instance from the given GTF file."""

        # Open gtf file.
        gtf = pysam.TabixFile(native_str(gtf_path), parser=pysam.asGTF())

        if chromosomes is None:
            chromosomes = gtf.contigs

        # Build the trees.
        transcript_trees = {}
        exon_trees = {}

        for chrom in chromosomes:
            # Collect exons and transcripts.
            transcripts = []
            exons = []

            records = gtf.fetch(reference=chrom)

            if record_filter is not None:
                records = (rec for rec in records if record_filter(rec))

            for record in records:
                if record.feature == 'transcript':
                    transcripts.append(cls._record_to_transcript(record))
                elif record.feature == 'exon':
                    exons.append(cls._record_to_exon(record))

            # Build transcript lookup tree.
            transcript_trees[chrom] = IntervalTree.from_tuples(
                (tr.start, tr.end, tr) for tr in transcripts)

            # Build exon lookup tree.
            keyfunc = lambda rec: rec.transcript_id

            exons = sorted(exons, key=keyfunc)
            grouped = itertools.groupby(exons, key=keyfunc)

            for tr_id, grp in grouped:
                exon_trees[tr_id] = IntervalTree.from_tuples(
                    (exon.start, exon.end, exon) for exon in grp) # yapf: disable

        return cls(transcript_trees, exon_trees)

    @staticmethod
    def _record_to_exon(record):
        attrs = record.asDict()
        return Exon(
            chromosome=record.contig,
            start=record.start,
            end=record.end,
            strand=1 if record.strand == '+' else -1,
            gene_name=attrs.get('gene_name', None),
            gene_id=attrs.get('gene_id', None),
            transcript_id=attrs.get('transcript_id', None))

    @staticmethod
    def _record_to_transcript(record):
        attrs = record.asDict()
        return Transcript(
            id=record['transcript_id'],
            chromosome=record.contig,
            start=record.start,
            end=record.end,
            strand=1 if record.strand == '+' else -1,
            gene_name=attrs.get('gene_name', None),
            gene_id=attrs.get('gene_id', None))

    @staticmethod
    def _lookup_genomic(trees, region):
        chrom, start, end = region
        overlap = trees[chrom][start:end]
        return [tup[2] for tup in overlap]

    @staticmethod
    def _lookup_tree(tree, region):
        start, end = region
        overlap = tree[start:end]
        return [tup[2] for tup in overlap]

    def overlap_transcripts(self, region, strict=True):
        # type: (Tuple[str, int, int], bool) -> List[Transcript]
        """Returns transcripts that overlap with given region."""

        overlap = self._lookup_genomic(self._transcripts, region)

        if strict:
            strict_overlap = []
            for transcript in overlap:
                exons = self._exons[transcript.id]
                exon_overlap = self._lookup_tree(exons, region[1:])

                if len(exon_overlap) > 0:
                    strict_overlap.append(transcript)
            overlap = strict_overlap

        return overlap

    def overlap_genes(self, region, strict=True):
        # type: (Tuple[str, int, int], bool) -> List[Gene]
        """Returns genes that overlap with given region."""

        transcripts = self.overlap_transcripts(region, strict=strict)

        id_func = operator.attrgetter('gene_id')
        transcripts = sorted(transcripts, key=id_func)
        grouped = itertools.groupby(transcripts, key=id_func)

        return [Gene.from_transcripts(list(grp)) for _, grp in grouped]

    def get_exons(self, transcript_id):
        # type: (str) -> List[Exon]
        """Returns exons for given transcript."""
        return [interval[2] for interval in self._exons[transcript_id].items()]


_Exon = namedtuple('Exon', [
    'chromosome', 'start', 'end', 'strand', 'gene_name', 'gene_id',
    'transcript_id'
])


class Exon(_Exon):
    """Exon model class, used in TranscriptReference instances."""
    __slots__ = ()


_Transcript = namedtuple('Transcript', [
    'id', 'chromosome', 'start', 'end', 'strand', 'gene_name', 'gene_id'
])


class Transcript(_Transcript):
    """Transcript model class, used in TranscriptReference instances."""
    __slots__ = ()


_Gene = namedtuple('Gene',
                   ['chromosome', 'start', 'end', 'strand', 'name', 'id'])


class Gene(_Gene):
    """Gene model class, used in TranscriptReference instances."""

    __slots__ = ()

    @classmethod
    def from_transcripts(cls, transcripts):
        # type: (List[Transcript]) -> Gene
        """Builds a gene instance from a collection of transcripts."""

        first = transcripts[0]
        start = min(transcript.start for transcript in transcripts)
        end = max(transcript.end for transcript in transcripts)

        return cls(
            chromosome=first.chromosome,
            start=start,
            end=end,
            strand=first.strand,
            name=first.gene_name,
            id=first.gene_id)
