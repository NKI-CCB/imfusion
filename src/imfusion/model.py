"""Module containing model classes for fusions and insertions."""

# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import collections

import pandas as pd
import numpy as np
import toolz


_Fusion = collections.namedtuple(
    'Fusion', ['seqname', 'anchor_genome', 'anchor_transposon',
               'strand_genome', 'strand_transposon', 'flank_genome',
               'flank_transposon', 'spanning_reads', 'supporting_mates',
               'spanning_mates', 'gene_id', 'gene_name', 'gene_strand',
               'feature_name', 'feature_type', 'feature_strand'])


class Fusion(_Fusion):
    """Class representing a gene-transposon fusion.

    Used by fusion identification tools (such as Tophat2) to return the
    fusions that are identified. Not all fields are required if these
    are not available. However, the following fields should at least
    be provided: `seqname`, `anchor_genome`, `anchor_transposon`,
    `strand_genome`, `strand_transposon`.

    Attributes:
        seqname (str): Chromosome involved in the fusion.
        anchor_genome (int): Genomic fusion breakpoint.
        anchor_transposon (int): Transposon fusion breakpoint.
        strand_genome (int): Strand of fusion in genome (-1 or 1)
        strand_transposon (int): Strand of fusion in transposon (-1 or 1)
        flank_genome (int): Size of flanking region on genome.
        flank_transposon (int): Size of flanking region in transposon.
        gene_id (str): ID of affected gene.
        gene_name (str): Name of affected gene.
        gene_strand (int): Strand of affected gene.
        feature_name (str): Name of affected transposon feature.
        feature_type (str): Feature type (SD/SA).
        feature_strand (int): Strand of affected transposon feature.
        spanning_reads (int): Number of supporting single-end reads.
        supporting_mates (int): Number of mate pairs that support the fusion,
            but do not span the breakpoint with either mate.
        supporting_spanning_mates (int): Number of mate pairs that support
            the fusion and have at least one mate spanning the breakpoint.

    """
    __slots__ = ()

    def __new__(cls, *args, **kwargs):
        defaults = {'flank_genome': np.nan, 'flank_transposon': np.nan,
                    'spanning_reads': np.nan, 'supporting_mates': np.nan,
                    'spanning_mates': np.nan, 'gene_id': None,
                    'gene_name': None, 'gene_strand': np.nan,
                    'feature_name': None, 'feature_type': None,
                    'feature_strand': np.nan}
        return super().__new__(cls, *args, **toolz.merge(defaults, kwargs))



_Insertion = collections.namedtuple(
    'Insertion', ['id', 'seqname', 'position', 'strand', 'sample_id',
                  'gene_id', 'gene_name', 'gene_strand', 'orientation',
                  'feature_name', 'feature_type', 'feature_strand',
                  'anchor_genome', 'anchor_transposon',
                  'flank_genome', 'flank_transposon',
                  'spanning_reads', 'supporting_mates', 'spanning_mates'])


class Insertion(_Insertion):
    """Class respresenting a RNA-seq transposon insertion site.

    Used to represent insertions derived from RNA-seq fusions. Not all fields
    are required, though at least the following should be specified: `id`,
    `seqname`, `position`, `strand` and `sample`.

    Attributes:
        id (str): ID of the insertion.
        seqname (str): Chromosome of the insertion.
        position (int): Genomic position of the insertion.
        strand (int): Strand of the insertion (-1 or 1).
        sample_id (str): Sample in which the insertion was identified.
        gene_id (str): ID of the gene involved in the fusion.
        gene_name (str): Name of the gene involved in the fusion.
        gene_strand (int): Strand of the gene involved in the fusion.
        orientation (str): Relative orientation of the insertion.
        feature_name (str): Name of the transposon feature involved
            in the fusion.
        feature_strand (int): Strand of transposon feature involved
            in the fusion.
        anchor_genome (int): Genomic fusion breakpoint.
        anchor_transposon (int): Transposon fusion breakpoint.
        flank_genome (int): Size of flanking region on genome.
        flank_transposon (int): Size of flanking region in transposon.
        spanning_reads (int): Number of supporting single-end reads.
        supporting_mates (int): Number of mate pairs that support the fusion,
            but do not span the breakpoint with either mate.
        spanning_mates (int): Number of mate pairs that support
            the fusion and have at least one mate spanning the breakpoint.

    """

    __slots__ = ()

    def __new__(cls, *args, **kwargs):
        defaults = {'anchor_genome': np.nan, 'anchor_transposon': np.nan,
                    'flank_genome': np.nan, 'flank_transposon': np.nan,
                    'spanning_reads': np.nan, 'supporting_mates': np.nan,
                    'spanning_mates': np.nan, 'gene_id': None,
                    'gene_name': None, 'gene_strand': np.nan,
                    'orientation': None, 'feature_name': None,
                    'feature_type': None, 'feature_strand': np.nan}
        return super().__new__(cls, *args, **toolz.merge(defaults, kwargs))

    @classmethod
    def read(cls, file_path):
        frame = cls.read_frame(file_path)
        for ins in cls.from_frame(frame):
            yield ins

    @classmethod
    def read_frame(cls, file_path):
        return pd.read_csv(str(file_path), sep='\t',
                           dtype={'seqname': str})

    @classmethod
    def from_frame(cls, ins_frame):
        for _, row in ins_frame.iterrows():
            yield cls(**row.to_dict())

    @classmethod
    def to_frame(cls, insertions):
        return pd.DataFrame.from_records(insertions, columns=cls._fields)
