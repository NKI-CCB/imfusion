# -*- coding: utf-8 -*-
"""Module containing model classes for fusions and insertions."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import collections
import itertools
import operator

import numpy as np
import pandas as pd
import toolz

from imfusion.util.frozendict import frozendict


class FrameMixin(object):
    @classmethod
    def _get_columns(cls):
        return cls._fields

    @classmethod
    def _to_record(cls, obj):
        return obj._asdict()

    @classmethod
    def _to_obj(cls, record):
        return cls(**record._asdict())

    @classmethod
    def read_csv(cls, file_path, **kwargs):
        """Reads objects from a csv file into a pandas DataFrame."""
        df = pd.read_csv(str(file_path), **kwargs)
        return cls.format_frame(df)

    @classmethod
    def to_csv(cls, file_path, objects, index=False, **kwargs):
        """Writes objects to a csv file."""
        df = cls.to_frame(objects)
        df.to_csv(str(file_path), index=index, **kwargs)

    @classmethod
    def from_csv(cls, file_path, **kwargs):
        """Reads objects from a csv file."""

        df = cls.read_csv(file_path, **kwargs)

        for obj in cls.from_frame(df):
            yield obj

    @classmethod
    def from_frame(cls, df):
        """Converts dataframe into an interable of objects."""
        for tup in df.itertuples():
            yield cls._to_obj(tup)

    @classmethod
    def format_frame(cls, df):
        """Formats dataframe into a cohorent format."""
        return cls._reorder_columns(df, order=cls._get_columns())

    @classmethod
    def to_frame(cls, objects):
        """Converts list of objects to a dataframe representation."""

        # Check if insertions is empty.
        is_empty, objects = cls._is_empty(objects)

        if is_empty:
            df = pd.DataFrame.from_records([], columns=cls._get_columns())
        else:
            rows = (cls._to_record(obj) for obj in objects)
            df = pd.DataFrame.from_records(rows)
            df = cls.format_frame(df)

        return df

    @staticmethod
    def _is_empty(iterable):
        try:
            _, iterable = toolz.peek(iterable)
            empty = False
        except StopIteration:
            empty = True

        return empty, iterable

    @classmethod
    def _reorder_columns(cls, df, order):
        extra_cols = set(df.columns) - set(order)
        col_order = list(order) + sorted(extra_cols)
        return df[col_order]


class MetadataFrameMixin(FrameMixin):
    """Mixin class adding namedtuple/frame conversion support."""

    @classmethod
    def _get_columns(cls):
        fields = list(cls._fields)
        del fields[fields.index('metadata')]
        return fields

    @classmethod
    def _to_record(cls, obj):
        obj_data = obj._asdict()
        metadata = obj_data.pop('metadata')
        return toolz.merge(metadata, obj_data)

    @classmethod
    def _to_obj(cls, record):
        record_dict = record._asdict()

        metadata_fields = [
            k for k in record_dict.keys() if k not in set(cls._get_columns())
        ]
        metadata = {k: record_dict.pop(k) for k in metadata_fields}

        metadata.pop('Index', None)

        return cls(metadata=frozendict(metadata), **record_dict)

    def __getattr__(self, name):
        if name in self.metadata:
            return self.metadata[name]
        else:
            raise AttributeError


_Fusion = collections.namedtuple('Fusion', [
    'seqname_a', 'location_a', 'strand_a', 'seqname_b', 'location_b',
    'strand_b', 'flank_a', 'flank_b', 'support_junction', 'support_spanning'
])


class Fusion(FrameMixin, _Fusion):
    """Model Fusion class."""

    __slots__ = ()

    @property
    def support(self):
        """Support score for the fusion."""
        return self.support_junction + self.support_spanning

    def normalize(self, seqname=None):
        """Normalizes fusions so that the side whose seqname has the
           lowest lexical ordering is used as donor.
        """

        if seqname is not None:
            if self.seqname_a == seqname:
                is_norm = True
            else:
                if self.seqname_b != seqname:
                    raise ValueError('Fusion does not include given seqname')
                is_norm = False
        elif self.seqname_a != self.seqname_b:
            is_norm = self.seqname_a < self.seqname_b
        else:
            is_norm = self.location_a < self.location_b

        if is_norm:
            return self
        else:
            return self._replace(
                seqname_a=self.seqname_b,
                location_a=self.location_b,
                strand_a=self.strand_b * -1,
                seqname_b=self.seqname_a,
                location_b=self.location_a,
                strand_b=self.strand_a * -1,
                flank_a=self.flank_b,
                flank_b=self.flank_a)

    def distance(self, other):
        """Determine distance to other fusion."""

        if (self.seqname_a != other.seqname_a or
                self.seqname_b != other.seqname_b):
            raise ValueError('Fusions are on different reference sequences')

        return (abs(self.location_a - other.location_a) +
                abs(self.location_b - other.location_b))

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
            return (fusion.seqname_a, fusion.strand_a, fusion.seqname_b,
                    fusion.strand_b)

        sorted_juncs = sorted(junctions, key=_keyfunc)
        grouped_juncs = itertools.groupby(sorted_juncs, key=_keyfunc)

        for _, grp in grouped_juncs:
            for grp_a in cls._groupby_position(grp, 'a', max_dist):
                for grp_b in cls._groupby_position(grp_a, 'b', max_dist):
                    yield grp_b

    @staticmethod
    def _groupby_position(junctions, side, max_dist):
        get_loc = operator.attrgetter('location_' + side)
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
    'seqname', 'anchor_genome', 'anchor_transposon', 'strand_genome',
    'strand_transposon', 'flank_genome', 'flank_transposon',
    'support_junction', 'support_spanning', 'metadata'
])


class TransposonFusion(MetadataFrameMixin, _TransposonFusion):
    """Model class representing a gene-transposon fusion."""

    __slots__ = ()

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
        return self.seqname, start, end

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

        if (fusion.seqname_a == transposon_name and
                fusion.seqname_b == transposon_name):
            raise ValueError('Fusion does not involve genomic sequence')
        elif (fusion.seqname_a != transposon_name and
              fusion.seqname_b != transposon_name):
            raise ValueError('Fusion does not involve transposon')
        elif fusion.seqname_a == transposon_name:
            tr_key, gen_key = 'a', 'b'
            tr_flank = fusion.flank_a * -fusion.strand_a
            gen_flank = fusion.flank_b * fusion.strand_b
        else:
            tr_key, gen_key = 'b', 'a'
            tr_flank = fusion.flank_b * fusion.strand_b
            gen_flank = fusion.flank_a * -fusion.strand_a

        return TransposonFusion(
            seqname=getattr(fusion, 'seqname_' + gen_key),
            anchor_genome=getattr(fusion, 'location_' + gen_key),
            anchor_transposon=getattr(fusion, 'location_' + tr_key),
            strand_genome=getattr(fusion, 'strand_' + gen_key),
            strand_transposon=getattr(fusion, 'strand_' + tr_key),
            flank_genome=gen_flank,
            flank_transposon=tr_flank,
            support_junction=fusion.support_junction,
            support_spanning=fusion.support_spanning,
            metadata=frozendict(metadata or {}))


_Insertion = collections.namedtuple('Insertion', [
    'id', 'seqname', 'position', 'strand', 'support_junction',
    'support_spanning', 'support', 'metadata'
])


class Insertion(MetadataFrameMixin, _Insertion):
    """Model class representing an insertion."""

    __slots__ = ()

    @classmethod
    def from_transposon_fusion(cls, fusion, id_=None, drop_metadata=None):
        """Converts (annotated) transposon fusion to an insertion.

        Requires

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

        if orientation is not None:
            ins_metadata['orientation'] = orientation

        return Insertion(
            id=id_,
            seqname=fusion.seqname,
            position=fusion.anchor_genome,
            strand=strand,
            support_junction=fusion.support_junction,
            support_spanning=fusion.support_spanning,
            support=fusion.support,
            metadata=frozendict(ins_metadata))

    @classmethod
    def from_transposon_fusions(cls,
                                fusions,
                                id_fmt_str=None,
                                drop_metadata=None):
        """Converts annotated transposon fusions to insertions."""

        if id_fmt_str is not None:
            insertions = (cls.from_transposon_fusion(
                fusion,
                id_=id_fmt_str.format(i + 1),
                drop_metadata=drop_metadata)
                          for i, fusion in enumerate(fusions))
        else:
            insertions = (cls.from_transposon_fusion(
                fusion, drop_metadata=drop_metadata) for fusion in fusions)

        for insertion in insertions:
            yield insertion
