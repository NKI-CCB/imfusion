# -*- coding: utf-8 -*-
"""Provides functions for annotating gene-transposon fusions."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

from collections import namedtuple
import gzip
import itertools
import operator
from typing import Any, Callable, Iterable, Tuple

import pathlib2 as pathlib

from future.utils import native_str
import numpy as np
import pysam
import toolz

from intervaltree import IntervalTree

from imfusion.model import MetadataFrameMixin, Insertion, Fusion
from imfusion.util import tabix
from imfusion.util.frozendict import frozendict


def extract_insertions(
        fusions,  # type: Iterable[Fusion]
        gtf_path,  # type: pathlib.Path
        features_path,  # type: pathlib.Path
        chromosomes=None,  # type: List[str]
        assembled_gtf_path=None,  # type: pathlib.Path
        ffpm_fastq_path=None  # type: pathlib.Path
):  # type: (...) -> Iterable[Insertion]
    """Extract insertions from gene-transposon fusions."""

    # Annotate for genes.
    gtf_reference = TranscriptReference.from_gtf(
        gtf_path, chromosomes=chromosomes)

    annotated = annotate_fusions_for_genes(fusions, gtf_reference)

    # Annotate for assembly (if given).
    if assembled_gtf_path is not None:
        assem_reference = TranscriptReference.from_gtf(
            assembled_gtf_path, chromosomes=chromosomes)

        annotated = annotate_fusions_for_assembly(annotated, gtf_reference,
                                                  assem_reference)

    # Annotate for transposon.
    annotated = annotate_fusions_for_transposon(annotated, features_path)

    # Drop any fusions without a transposon feature.
    annotated = (fusion for fusion in annotated
                 if 'feature_name' in fusion.metadata)

    # Calculate FFPM scores.
    if ffpm_fastq_path is not None:
        annotated = annotate_ffpm(annotated, fastq_path=ffpm_fastq_path)

    # Convert to insertions.
    insertions = Insertion.from_transposon_fusions(
        annotated, id_fmt_str='INS_{}')

    for insertion in insertions:
        yield insertion


def annotate_fusions_for_genes(fusions, reference):
    # type: (Iterable[Fusion], TranscriptReference) -> Iterable[Fusion]
    """Annotates fusions with genes overlapped by the genomic fusion site.

    Parameters
    ----------
    fusions : iterable[TransposonFusion]
        Fusions to annotate.
    gtf_path : pathlib.Path
        Path to (indexed) gtf file, containing gene exon features.

    Yields
    ------
    Fusion
        Fusions, annotated with overlapped genes.

    """

    for fusion in fusions:
        genes = reference.overlap_genes(fusion.genome_region)

        if len(genes) > 0:
            for gene in genes:
                gene_meta = {
                    'gene_name': gene.name,
                    'gene_strand': gene.strand,
                    'gene_id': gene.id
                }
                merged_meta = toolz.merge(fusion.metadata, gene_meta)
                yield fusion._replace(metadata=frozendict(merged_meta))
        else:
            yield fusion


def annotate_fusions_for_transposon(fusions, feature_path):
    # type: (Iterable[Fusion], pathlib.Path) -> Iterable[Fusion]
    """Annotates fusions with transposon features overlapped by the fusion.

    Parameters
    ----------
    fusions : iterable[TransposonFusion]
        Fusions to annotate.
    feature_path : str or pathlib.Path
        Path to TSV file containing transposon features.

    Yields
    ------
    Fusion
        Fusions, annotated with transposon features.

    """

    def _lookup_tree(tree, region):
        start, end = region
        overlap = tree[start:end]
        return [tup[2] for tup in overlap]

    # Build the feature tree.
    features = TransposonFeature.from_csv(feature_path, sep='\t')
    tree = IntervalTree.from_tuples((feat.start, feat.end, feat)
                                    for feat in features)

    for fusion in fusions:
        overlap = _lookup_tree(tree, fusion.transposon_region)

        if len(overlap) > 0:
            for feature in overlap:
                new_meta = {
                    'feature_name': feature.name,
                    'feature_type': feature.type,
                    'feature_strand': feature.strand
                }
                merged_meta = toolz.merge(fusion.metadata, new_meta)
                yield fusion._replace(metadata=frozendict(merged_meta))
        else:
            yield fusion


def annotate_fusions_for_assembly(
        fusions,  # type: Iterable[Fusion]
        reference,  # type: TranscriptReference
        assembly,  # type: TranscriptReference
        skip_annotated=True  # type: bool
):  # type: (...) -> Iterable[Fusion]
    """Annotates fusions using the assembled GTF."""

    def _exon_region(exon):
        return (exon.chromosome, exon.start, exon.end)

    for fusion in fusions:
        if skip_annotated and 'gene_id' in fusion.metadata:
            # Already annotated
            yield fusion
        else:
            # Identify overlapped transcripts.
            transcripts = assembly.overlap_transcripts(fusion.genome_region)

            if len(transcripts) > 0:
                for transcript in transcripts:
                    # Lookup genes that overlap with exons.
                    exons = assembly.get_exons(transcript.id)

                    genes = set(
                        itertools.chain.from_iterable(
                            reference.overlap_genes(_exon_region(exon))
                            for exon in exons))

                    if len(genes) > 0:
                        for gene in genes:
                            # Yield with information from overlapping genes.
                            new_meta = {
                                'gene_name': gene.name,
                                'gene_strand': gene.strand,
                                'gene_id': gene.id,
                                'novel_transcript': transcript.id
                            }
                            yield fusion._replace(metadata=toolz.merge(
                                fusion.metadata, new_meta))
                    else:
                        # No gene overlap, yield with transcript info.
                        new_meta = {
                            'gene_name': transcript.id,
                            'gene_id': transcript.id,
                            'gene_strand': transcript.strand,
                            'novel_transcript': transcript.id
                        }
                        yield fusion._replace(
                            metadata=toolz.merge(fusion.metadata, new_meta))
            else:
                # No overlap.
                yield fusion


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
                    (exon.start, exon.end, exon) for exon in grp)

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


_Transcript = namedtuple(
    'Transcript',
    ['id', 'chromosome', 'start', 'end', 'strand', 'gene_name', 'gene_id'])


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

        return cls(chromosome=first.chromosome,
                   start=start,
                   end=end,
                   strand=first.strand,
                   name=first.gene_name,
                   id=first.gene_id)


_TransposonFeature = namedtuple(
    'TransposonFeature',
    ['name', 'start', 'end', 'strand', 'type', 'metadata'])


class TransposonFeature(_TransposonFeature, MetadataFrameMixin):
    """Transposon feature model class."""

    __slots__ = ()


def filter_insertions(
        insertions,  # type: Iterable[Insertion]
        features=True,  # type: bool
        orientation=True,  # type: bool
        blacklist=None  # type: Set[str]
):
    """Filters false positive insertions using a common set of filters.

    Parameters
    ----------
    insertions : Iterable[Insertion]
        Insertions to filter.
    features : bool
        Whether to filter insertions that correspond to unexpected
        features of the transposon (non SA/SD) features.
    orientation : bool
        Whether to filter insertions that have the conflicting orientations
        between gene and transposon features.
    blacklist : set[str]
        List of blacklisted genes to filter for.

    Yields
    ------
    Insertion
        Next filtered insertion.

    """

    if features:
        insertions = filter_unexpected_features(insertions)

    if blacklist is not None:
        insertions = filter_blacklist(insertions, blacklist)

    if orientation:
        insertions = filter_wrong_orientation(insertions)

    for insertion in insertions:
        yield insertion


def filter_unexpected_features(insertions):
    # type: (Iterable[Insertion]) -> Iterable[Insertion]
    """Filters insertions that have non splice-acceptor/donor features.

    This filter removes any insertions that splice to tranposon features
    that aren't splice-acceptors or splice-donors. This is based on the
    premise that these other sites are unlikely to be involved in any
    splicing and that therefore these insertions are likely to be false
    positives of the fusion identification.

    Parameters
    ----------
    insertions : iterable[Insertion]
        Insertions to filter.

    Yields
    ------
    Insertion
        Next filtered insertion.

    """

    for ins in insertions:
        if ins.metadata['feature_type'] in {'SA', 'SD'}:
            yield ins


def filter_blacklist(insertions, genes, field='gene_name'):
    # type: (Iterable[Insertion], Set[str], str) -> Iterable[Insertion]
    """Filters insertions for blacklisted genes.

    Parameters
    ----------
    insertions : List[Insertion]
        Insertions to filter.
    genes : set[str]
        Symbols of the blacklisted genes.

    Yields
    ------
    Insertion
        Next filtered insertion.

    """

    for ins in insertions:
        if ins.metadata[field] not in genes:
            yield ins


def filter_wrong_orientation(insertions, drop_na=False):
    # type: (Iterable[Insertion], bool) -> Iterable[Insertion]
    """Filters insertions with wrong feature orientations w.r.t. their genes.

    This filter removes any insertions with a transposon feature that is
    in the wrong orientation with respect to the annotated gene. This is
    based on the premise that, for example, a splice acceptor can only
    splice to a gene that is in the same orientation as the acceptor.

    Parameters
    ----------
    insertions : List[Insertion]
        Insertions to filter.

    Yields
    ------
    Insertion
        Next filtered insertion.

    """

    for ins in insertions:
        feat_strand = ins.metadata.get('feature_strand', np.nan)
        gene_strand = ins.metadata.get('gene_strand', np.nan)

        feat_ori = ins.strand * feat_strand

        if feat_ori == gene_strand:
            # Correct orientation, yield.
            yield ins
        elif not drop_na and (np.isnan(feat_ori) or np.isnan(gene_strand)):
            # Keeping NaN cases, so yield.
            yield ins


def annotate_ffpm(fusions, fastq_path):
    # type: (Iterable[Fusion], pathlib.Path) -> Iterable[Fusion]
    """Annotates fusions with FFPM (Fusion Fragments Per Million) score."""

    # Calculate normalization factor.
    n_reads = count_lines(fastq_path) // 4
    norm_factor = (1.0 / n_reads) * 1e6

    for fusion in fusions:
        ffpm_meta = {
            'ffpm_junction': fusion.support_junction * norm_factor,
            'ffpm_spanning': fusion.support_spanning * norm_factor,
            'ffpm': fusion.support * norm_factor
        }
        merged_meta = toolz.merge(fusion.metadata, ffpm_meta)
        yield fusion._replace(metadata=frozendict(merged_meta))


def count_lines(file_path):
    # type: (pathlib.Path) -> int
    """Counts number of lines in (gzipped) file."""

    if file_path.suffixes[-1] == '.gz':
        with gzip.open(str(file_path)) as file_obj:
            count = _count_lines(file_obj)
    else:
        with file_path.open() as file_obj:
            count = _count_lines(file_obj)
    return count


def _count_lines(file_obj):
    """Counts number of lines in given file."""

    lines = 0
    buf_size = 1024 * 1024
    read_f = file_obj.read  # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)

    return lines
