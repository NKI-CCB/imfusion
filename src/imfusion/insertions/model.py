# -*- coding: utf-8 -*-
"""Module containing model classes for fusions and insertions."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import collections
import gzip
import itertools
import operator
import warnings

from future.utils import native_str
import numpy as np
import pandas as pd
import toolz

from imfusion.util.pandas_ import MetadataRecordSet
from imfusion.vendor.frozendict import frozendict
from imfusion.vendor.genopandas import GenomicDataFrame

from .reference import TranscriptReference

_Fusion = collections.namedtuple('Fusion', [
    'chromosome_a', 'position_a', 'strand_a', 'chromosome_b', 'position_b',
    'strand_b', 'flank_a', 'flank_b', 'support_junction', 'support_spanning',
    'sample'
])


class Fusion(_Fusion):
    """Model Fusion class."""

    __slots__ = ()

    @property
    def support(self):
        """Support score for the fusion."""
        return self.support_junction + self.support_spanning

    def normalize(self, chromosome=None):
        """Normalizes fusions so that the side whose chromosome has the
           lowest lexical ordering is used as donor.
        """

        if chromosome is not None:
            if self.chromosome_a == chromosome:
                is_norm = True
            else:
                if self.chromosome_b != chromosome:
                    raise ValueError(
                        'Fusion does not include given chromosome')
                is_norm = False
        elif self.chromosome_a != self.chromosome_b:
            is_norm = self.chromosome_a < self.chromosome_b
        else:
            is_norm = self.position_a < self.position_b

        if is_norm:
            return self
        else:
            return self._replace(
                chromosome_a=self.chromosome_b,
                position_a=self.position_b,
                strand_a=self.strand_b * -1,
                chromosome_b=self.chromosome_a,
                position_b=self.position_a,
                strand_b=self.strand_a * -1,
                flank_a=self.flank_b,
                flank_b=self.flank_a)

    def distance(self, other):
        """Determine distance to other fusion."""

        if (self.chromosome_a != other.chromosome_a
                or self.chromosome_b != other.chromosome_b):
            raise ValueError('Fusions are on different reference sequences')

        return (abs(self.position_a - other.position_a) +
                abs(self.position_b - other.position_b))

    @classmethod
    def merge(cls, junctions, max_dist):
        """Merges fusion within genomic distance into a single fusion."""
        for grp in cls._groupby_genomic_position(junctions, max_dist):
            if len(grp) == 1:
                yield grp[0]
            else:
                get_support = operator.attrgetter('support_junction')
                sorted_grp = sorted(grp, key=get_support)

                yield sorted_grp[-1]._replace(
                    support_junction=sum(f.support_junction
                                         for f in sorted_grp),
                    support_spanning=sum(f.support_spanning
                                         for f in sorted_grp),
                    flank_a=max(f.flank_a for f in sorted_grp),
                    flank_b=max(f.flank_b for f in sorted_grp))

    @classmethod
    def _groupby_genomic_position(cls, junctions, max_dist):

        # Group junctions by strand and sequence.
        def _keyfunc(fusion):
            return (fusion.chromosome_a, fusion.strand_a, fusion.chromosome_b,
                    fusion.strand_b)

        sorted_juncs = sorted(junctions, key=_keyfunc)
        grouped_juncs = itertools.groupby(sorted_juncs, key=_keyfunc)

        for _, grp in grouped_juncs:
            for grp_a in cls._groupby_position(grp, 'a', max_dist):
                for grp_b in cls._groupby_position(grp_a, 'b', max_dist):
                    yield grp_b

    @staticmethod
    def _groupby_position(junctions, side, max_dist):
        get_loc = operator.attrgetter('position_' + side)
        sorted_juncs = sorted(junctions, key=get_loc)

        grp, prev_pos = [], np.nan
        for junc in sorted_juncs:
            if (get_loc(junc) - prev_pos) > max_dist:
                yield grp
                grp = [junc]
            else:
                grp.append(junc)
            prev_pos = get_loc(junc)

        if len(grp) > 0:
            yield grp


_TransposonFusion = collections.namedtuple('TransposonFusion', [
    'chromosome', 'anchor_genome', 'anchor_transposon', 'strand_genome',
    'strand_transposon', 'flank_genome', 'flank_transposon',
    'support_junction', 'support_spanning', 'sample', 'metadata'
])


class TransposonFusion(_TransposonFusion):
    """Model class representing a gene-transposon fusion."""

    __slots__ = ()

    def __new__(cls, chromosome, anchor_genome, anchor_transposon,
                strand_genome, strand_transposon, flank_genome,
                flank_transposon, support_junction, support_spanning, sample,
                **kwargs):

        return super().__new__(
            cls,
            chromosome=chromosome,
            anchor_genome=int(anchor_genome),
            anchor_transposon=int(anchor_transposon),
            strand_genome=int(strand_genome),
            strand_transposon=int(strand_transposon),
            flank_genome=int(flank_genome),
            flank_transposon=int(flank_transposon),
            support_junction=int(support_junction),
            support_spanning=int(support_spanning),
            sample=sample,
            metadata=frozendict(kwargs))

    def __getattr__(self, name):
        return self.metadata[name]

    @property
    def support(self):
        """Support score for the fusion."""
        return self.support_junction + self.support_spanning

    @property
    def genome_region(self):
        """Returns the spanned region on the genome."""
        if self.flank_genome > 0:
            start = self.anchor_genome
            end = self.anchor_genome + self.flank_genome
        else:
            start = self.anchor_genome + self.flank_genome
            end = self.anchor_genome
        return self.chromosome, start, end

    @property
    def transposon_region(self):
        """Returns the spanned region on the transposon."""
        if self.flank_transposon > 0:
            start = self.anchor_transposon
            end = self.anchor_transposon + self.flank_transposon
        else:
            start = self.anchor_transposon + self.flank_transposon
            end = self.anchor_transposon
        return start, end

    @classmethod
    def from_fusion(cls, fusion, transposon_name, metadata=None):
        """Converts the fusion to a transposon fusion object."""

        if (fusion.chromosome_a == transposon_name
                and fusion.chromosome_b == transposon_name):
            raise ValueError('Fusion does not involve genomic sequence')
        elif (fusion.chromosome_a != transposon_name
              and fusion.chromosome_b != transposon_name):
            raise ValueError('Fusion does not involve transposon')
        elif fusion.chromosome_a == transposon_name:
            tr_key, gen_key = 'a', 'b'
            tr_flank = fusion.flank_a * -fusion.strand_a
            gen_flank = fusion.flank_b * fusion.strand_b
        else:
            tr_key, gen_key = 'b', 'a'
            tr_flank = fusion.flank_b * fusion.strand_b
            gen_flank = fusion.flank_a * -fusion.strand_a

        return cls(
            chromosome=getattr(fusion, 'chromosome_' + gen_key),
            anchor_genome=getattr(fusion, 'position_' + gen_key),
            anchor_transposon=getattr(fusion, 'position_' + tr_key),
            strand_genome=getattr(fusion, 'strand_' + gen_key),
            strand_transposon=getattr(fusion, 'strand_' + tr_key),
            flank_genome=gen_flank,
            flank_transposon=tr_flank,
            support_junction=fusion.support_junction,
            support_spanning=fusion.support_spanning,
            sample=fusion.sample,
            **(metadata or {}))

    @staticmethod
    def annotate_for_genes(fusions, reference):
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

    @staticmethod
    def annotate_for_transposon(fusions, feature_path):
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

        # Build transposon feature frame.
        features = pd.read_csv(native_str(feature_path), sep='\t')
        features['chromosome'] = 'transposon'
        features = features.set_index(['chromosome', 'start', 'end'])

        features = GenomicDataFrame(features)

        # Lookup annotations.
        for fusion in fusions:
            overlap = features.gloc.search('transposon',
                                           *fusion.transposon_region)

            if overlap.shape[0] > 0:
                for feature in overlap.itertuples():
                    new_meta = {
                        'feature_name': feature.name,
                        'feature_type': feature.type,
                        'feature_strand': feature.strand
                    }

                    merged_meta = toolz.merge(fusion.metadata, new_meta)
                    yield fusion._replace(metadata=frozendict(merged_meta))
            else:
                yield fusion

    @staticmethod
    def annotate_for_assembly(fusions,
                              reference,
                              assembly,
                              skip_annotated=True):
        """Annotates fusions using the assembled GTF."""

        def _exon_region(exon):
            return (exon.chromosome, exon.start, exon.end)

        for fusion in fusions:
            if skip_annotated and 'gene_id' in fusion.metadata:
                # Already annotated
                yield fusion
            else:
                # Identify overlapped transcripts.
                transcripts = assembly.overlap_transcripts(
                    fusion.genome_region)

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
                            yield fusion._replace(metadata=toolz.merge(
                                fusion.metadata, new_meta))
                else:
                    # No overlap.
                    yield fusion

    @staticmethod
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

    def _count_lines_obj(file_obj):
        """Counts number of lines in given file object."""

        lines = 0
        buf_size = 1024 * 1024
        read_f = file_obj.read  # loop optimization

        buf = read_f(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = read_f(buf_size)

        return lines

    if file_path.suffixes[-1] == '.gz':
        with gzip.open(str(file_path)) as file_obj:
            count = _count_lines_obj(file_obj)
    else:
        with file_path.open() as file_obj:
            count = _count_lines_obj(file_obj)
    return count


_Insertion = collections.namedtuple('Insertion', [
    'id', 'chromosome', 'position', 'strand', 'support', 'sample', 'metadata'
])


class Insertion(_Insertion):
    """Model class representing an insertion."""

    __slots__ = ()

    def __new__(cls, id, chromosome, position, strand, support, sample,
                **kwargs):

        return super().__new__(
            cls,
            id=id,
            chromosome=chromosome,
            position=int(position),
            strand=int(strand),
            support=int(support),
            sample=sample,
            metadata=frozendict(kwargs))

    def __getattr__(self, name):
        return self.metadata[name]

    @classmethod
    def from_transposon_fusion(cls, fusion, id=None, drop_metadata=None):
        """Converts (annotated) transposon fusion to an insertion.

        Parameters
        ----------
        id : str
            ID to use for the insertion.

        Returns
        -------
        Insertion
            Insertion object derived from the transposon fusion.

        """

        if drop_metadata is None:
            drop_metadata = {'strand_genome', 'strand_transposon'}

        # Fusion should have annotated transposon feature.
        if 'feature_name' not in fusion.metadata:
            raise ValueError(
                'Fusion does not have an annotated transposon feature')

        strand = fusion.strand_genome * fusion.strand_transposon

        if 'gene_strand' in fusion.metadata:
            gene_strand = fusion.metadata['gene_strand']
            orientation = 'sense' if strand == gene_strand else 'antisense'
        else:
            orientation = None

        ins_metadata = toolz.keyfilter(lambda k: k not in drop_metadata,
                                       fusion.metadata)

        ins_metadata['transposon_anchor'] = fusion.anchor_transposon
        ins_metadata['support_junction'] = fusion.support_junction
        ins_metadata['support_spanning'] = fusion.support_spanning

        if orientation is not None:
            ins_metadata['orientation'] = orientation

        return Insertion(
            id=id,
            chromosome=fusion.chromosome,
            position=fusion.anchor_genome,
            strand=strand,
            support=fusion.support,
            sample=fusion.sample,
            **ins_metadata)

    @classmethod
    def from_annotated_transposon_fusions(cls,
                                          fusions,
                                          id_fmt_str=None,
                                          drop_metadata=None):
        """Converts annotated transposon fusions to insertions."""

        if id_fmt_str is not None:
            for i, fusion in enumerate(fusions):
                yield cls.from_transposon_fusion(
                    fusion,
                    id=id_fmt_str.format(num=i + 1, sample=fusion.sample),
                    drop_metadata=drop_metadata)
        else:
            for fusion in fusions:
                yield cls.from_transposon_fusion(
                    fusion, drop_metadata=drop_metadata)

    @classmethod
    def from_transposon_fusions(
            cls,
            fusions,  # type: Iterable[TransposonFusion]
            gtf_path,  # type: pathlib.Path
            features_path,  # type: pathlib.Path
            chromosomes=None,  # type: List[str]
            assembled_gtf_path=None,  # type: pathlib.Path
            ffpm_fastq_path=None,  # type: pathlib.Path,
            id_fmt_str='{sample}.INS_{num}'  # type: str
    ):  # type: (...) -> Iterable[Insertion]
        """Extract insertions from gene-transposon fusions."""

        # Annotate for genes.
        gtf_reference = TranscriptReference.from_gtf(
            gtf_path, chromosomes=chromosomes)

        annotated = TransposonFusion.annotate_for_genes(fusions, gtf_reference)

        # Annotate for assembly (if given).
        if assembled_gtf_path is not None:
            assembled_reference = TranscriptReference.from_gtf(
                assembled_gtf_path, chromosomes=chromosomes)

            annotated = TransposonFusion.annotate_for_assembly(
                annotated, gtf_reference, assembled_reference)

        # Annotate for transposon.
        annotated = TransposonFusion.annotate_for_transposon(
            annotated, features_path)

        # Drop any fusions without a transposon feature.
        annotated = (fusion for fusion in annotated
                     if 'feature_name' in fusion.metadata)

        # Calculate FFPM scores.
        if ffpm_fastq_path is not None:
            annotated = TransposonFusion.annotate_ffpm(
                annotated, fastq_path=ffpm_fastq_path)

        # Convert to insertions.
        insertions = cls.from_annotated_transposon_fusions(
            annotated, id_fmt_str=id_fmt_str)

        for insertion in insertions:
            yield insertion


class InsertionSet(MetadataRecordSet):
    """Class that represents an insertion dataset."""

    @property
    def samples(self):
        return set(self._values['sample'])

    @classmethod
    def _tuple_class(cls):
        return Insertion

    @classmethod
    def from_csv(cls, file_path, **kwargs):
        """Reads a InsertionSet from a csv file using pandas.read_csv."""
        values = pd.read_csv(native_str(file_path), **kwargs)

        if 'chromosome' not in values and 'seqname' in values:
            warnings.warn('Using \'seqname\' to describe insertion positions '
                          'is deprecated, use \'chromosome\' instead',
                          DeprecationWarning)
            values = values.rename(columns={'seqname': 'chromosome'})

        return cls(values)

    @classmethod
    def concat(cls, insertion_sets, check_overlap=True):
        """Concatenates multiple records sets into a single set."""

        seen = set()
        if check_overlap:
            for ins_set in insertion_sets:
                if len(ins_set) > 0 and (seen & ins_set.samples):
                    raise ValueError('Duplicate samples in insertion sets')
                seen |= ins_set.samples

        return super().concat(insertion_sets)
