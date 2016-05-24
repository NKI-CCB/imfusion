"""Module containing functions for filtering insertions."""

# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import numpy as np

from imfusion.util.tabix import GtfFile


def filter_invalid_insertions(insertions):
    """Filters invalid insertions.

    Main function for filtering invalid insertions. Effectively applies
    both the `filter_wrong_orientation` and `filter_unexpected_sites` filters
    to filter insertions that have the wrong orientation (w.r.t the transposon
    feature and the annotated gene) or involve features of the transposon
    that we are not interested in (typically non-splice acceptor/donor
    features).

    Parameters
    ----------
    insertions : List[Insertion]
        Insertions to filter.

    Yields
    ------
    Insertion
        Next filtered insertion.

    """
    # Filter fusions based on characteristics that
    # reflect implausible or unwanted situations.
    insertions = filter_wrong_orientation(insertions)
    insertions = filter_unexpected_sites(insertions)

    # Yield remaining fusions.
    for ins in insertions:
        yield ins


def filter_wrong_orientation(insertions, drop_na=False):
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
        feat_ori = ins.strand * ins.feature_strand

        if feat_ori == ins.gene_strand:
            # Correct orientation, yield.
            yield ins
        elif not drop_na and (np.isnan(feat_ori) or np.isnan(ins.gene_strand)):
            # Keeping NaN cases, so yield.
            yield ins


def filter_unexpected_sites(insertions):
    """Filters insertions that have non splice-acceptor/donor features.

    This filter removes any insertions that splice to tranposon features
    that aren't splice-acceptors or splice-donors. This is based on the
    premise that these other sites are unlikely to be involved in any
    splicing and that therefore these insertions are likely to be false
    positives of the fusion identification.

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
        if ins.feature_name.endswith('SA') or \
                ins.feature_name.endswith('SD'):
            yield ins


def filter_blacklist(insertions, gene_ids, reference_gtf=None,
                     filter_overlap=True):
    """Filters insertions for blacklisted genes.

    Parameters
    ----------
    insertions : List[Insertion]
        Insertions to filter.
    gene_ids : set[str]
        IDs of the blacklisted genes.
    reference_gtf : GtfFile
        GtfFile instance containing reference gene features.
        (Only needed if filter_overlap is True).
    filter_overlap : bool
        Whether to filter any insertions overlapping with the listed genes.
        If False (default), only genes explicitly splicing to the
        gene are filtered.

    Yields
    ------
    Insertion
        Next filtered insertion.

    """

    if reference_gtf is not None and not isinstance(reference_gtf, GtfFile):
        reference_gtf = GtfFile(reference_gtf)

    # Filter insertions that are directly splicing.
    filtered = (ins for ins in insertions
                if ins.gene_id is None or
                ins.gene_id not in gene_ids)

    if filter_overlap:
        # Check if gtf is given.
        if reference_gtf is None:
            raise ValueError('Gtf must be given if filter_overlap is True')

        # Filter overlapping insertions.
        for gene_id in gene_ids:
            gene = reference_gtf.get_gene(gene_id)
            filtered = _filter_overlap(filtered, gene=gene)

    # Yield remaining insertions.
    for ins in filtered:
        yield ins


def _filter_overlap(insertions, gene):
    """Filters insertions overlapping gene."""

    for ins in insertions:
        if not (ins.seqname == gene.contig and
                ins.position >= gene.start and
                ins.position <= gene.end):
            yield ins
