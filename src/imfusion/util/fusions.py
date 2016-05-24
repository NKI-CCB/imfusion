"""Module containing functions for annotating gene-transposon fusions."""

# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import numpy as np

from imfusion.model import Insertion
from imfusion.util.tabix import GtfFile


def annotate_fusions(fusions, reference_gtf, transposon_features):
    """Annotates fusions with gene and transposon features.

    Main function for annotating identified gene-transposon fusions.
    Adds the following annotations to fusions: gene features, transposon
    features and the relative orientation of the corresponding insertion
    with respect to the identified gene.

    Parameters
    ----------
    fusions : List[Fusion]
        Fusions to annotate.
    reference_gtf : GtfFile
        GtfFile instance containing reference gene features.
    transposon_features : pandas.DataFrame
        Dataframe containing positions for the features present in
        the transposon sequence. Used to identify transposon features
        (such as splice acceptors or donors) that are involved in the
        identified fusions.

    Yields
    -------
    Fusion
        Next fusion, annotated with gene/transposon features.

    """

    # Annotate fusions with features.
    fusions = annotate_with_gene(fusions, reference_gtf)
    fusions = annotate_with_transposon(fusions, transposon_features)

    for fusion in fusions:
        yield fusion


def annotate_with_gene(fusions, reference_gtf):
    """Annotates fusions with genes overlapped by the genomic fusion site.

    Parameters
    ----------
    fusions : List[Fusion]
        Fusions to annotate.
    reference_gtf : GtfFile
        GtfFile instance containing reference gene features.

    Yields
    -------
    Fusion
        Next fusion, annotated with gene features (includes `gene_id`,
        `gene_name` and `gene_strand`).

    """

    # Get gtf file.
    # gtf = GtfFile(reference_gtf)

    # Yield annotated fusions.
    for fusion in fusions:
        for annotated_fusion in _annotate_genes(fusion, reference_gtf):
            yield annotated_fusion


def _annotate_genes(fusion, ref_gtf):
    """Annotates a fusion with overlapping genes."""

    # Define region for feature search.
    if fusion.strand_genome == 1:
        start = fusion.anchor_genome - fusion.flank_genome
        end = fusion.anchor_genome
    else:
        start = fusion.anchor_genome
        end = fusion.anchor_genome + fusion.flank_genome

    # Identify overlapped genes.
    exons = ref_gtf.get_region(fusion.seqname, start, end,
                               filters={'feature': 'exon'})

    if len(exons) > 0:
        # Get corresponding gene details.
        genes = set(zip(exons['gene_id'], exons['gene_name'], exons['strand']))

        # Return annotated fusions.
        for gene_id, gene_name, gene_strand in genes:
            yield fusion._replace(gene_id=gene_id, gene_name=gene_name,
                                  gene_strand=_numeric_strand(gene_strand))
    else:
        # Return original, unannotated fusion.
        yield fusion



def _numeric_strand(strand):
    if strand == '+':
        return 1
    elif strand == '-':
        return -1
    else:
        return np.nan


def annotate_with_transposon(fusions, transposon_features):
    """Annotates fusions with transposon features overlapped by the fusion.

    Parameters
    ----------
    fusions : List[Fusion]
        Fusions to annotate, in DataFrame format.
    transposon_features : pandas.DataFrame
        Dataframe containing positions for the features present in
        the transposon sequence. Used to identify transposon features
        (such as splice acceptors or donors) that are involved in the
        identified fusions.

    Yields
    -------
    Fusion
        Next fusion, annotated with transposon features.

    """

    for fusion in fusions:
        for annotated_fusion in _annotate_features(fusion, transposon_features):
            yield annotated_fusion


def _annotate_features(fusion, transposon_features):
    # Identify overlapped features.
    features = transposon_features.ix[
        (transposon_features['start'] < fusion.anchor_transposon) &
        (transposon_features['end'] > fusion.anchor_transposon)]

    if len(features) > 0:
        for _, feature in features.iterrows():
            yield fusion._replace(feature_name=feature['name'],
                                  feature_type=feature['type'],
                                  feature_strand=feature['strand'])
    else:
        yield fusion


def place_fusions(fusions, sample_id, reference_gtf, offset=20, max_dist=5000):
    """Derives insertions by placing fusions at approximate genomic locations.

    Main function for deriving insertions from annotated gene-transposon
    fusions. Derives insertions by determining an approximate genomic
    location that is compatible with the gene/transposon feature annotations
    of the fusions. Fusions are therefore expected to be properly annotated
    for gene/transposon features.

    An insertion is essentially 'placed' by looking for the first genomic
    position that does not overlap with a reference feature, in the direction
    that is compatible with the insertions orientation w.r.t. its target gene.

    Parameters
    ----------
    fusions : List[Fusion]
        List of fusions to convert.
    sample_id : str
        Sample id that should be used for the insertions.
    reference_gtf : GtfFile
        GtfFile containing the reference features. Expected to conform
        to the Ensembl reference gtf format.
    offset : int
        Minimum offset of the transposon to the closest reference gene feature.
    max_dist : int
        Maximum distance that an insertion may be placed from the
        genomic anchor of the fusion.

    Yields
    -------
    Insertion
        Next insertion derived from the given fusions.

    """

    for i, fusion in enumerate(fusions):
        yield place_fusion(fusion, sample_id, reference_gtf, offset,
                           max_dist, id_='INS_{}'.format(i + 1))


def place_fusion(fusion, sample_id, reference_gtf,
                 offset=20, max_dist=5000, id_=None):
    """Derives an insertion from a given fusion.

    Parameters
    ----------
    fusion : Fusion
        Fusions to derive an insertion for.
    sample_id : str
        Sample id that should be used for the insertion.
    reference_gtf : GtfFile
        GtfFile containing the reference features. Expected to conform
        to the Ensembl reference gtf format.
    offset : int
        Minimum offset of the transposon to the closest reference gene feature.
    max_dist : int
        Maximum distance that an insertion may be placed from the
        genomic anchor of the fusion.
    id_: str
        ID to use for the insertion.

    Returns
    -------
    Insertion
        The insertion derived from the given fusion.
    """

    strand, orientation = insertion_orientation(fusion)
    position = insertion_position(fusion, reference_gtf, strand, orientation,
                                  offset=offset, max_dist=max_dist)

    insertion = Insertion(
        id=id_,
        seqname=fusion.seqname,
        position=position,
        strand=strand,
        sample_id=sample_id,
        gene_id=fusion.gene_id,
        gene_name=fusion.gene_name,
        gene_strand=fusion.gene_strand,
        orientation=orientation,
        feature_name=fusion.feature_name,
        feature_type=fusion.feature_type,
        feature_strand=fusion.feature_strand,
        anchor_genome=fusion.anchor_genome,
        anchor_transposon=fusion.anchor_transposon,
        flank_genome=fusion.flank_genome,
        flank_transposon=fusion.flank_transposon,
        spanning_reads=fusion.spanning_reads,
        supporting_mates=fusion.supporting_mates,
        spanning_mates=fusion.spanning_mates)

    return insertion


def insertion_position(fusion, reference_gtf, strand, orientation,
                       offset=20, max_dist=5000):
    """Determines an approximate genomic location for an RNA-seq insertion.

    This position is determined by looking for the first genomic
    position that does not overlap with a reference feature, in the direction
    that is compatible with the insertions orientation w.r.t. its target gene.

    Parameters
    ----------
    fusions : Fusion
        Fusion to find a position fore.
    reference_gtf : GtfFile
        GtfFile containing the reference features. Expected to conform
        to the Ensembl reference gtf format.
    strand : int
        Strand of the corresponding insertion (1 or -1).
    orientation : str
        Orientation of the corresponding insertion with respect to the
        annotated gene (sense or antisense).
    offset : int
        Minimum offset of the transposon to the closest reference gene feature.
    max_dist : int
        Maximum distance that an insertion may be placed from the
        genomic anchor of the fusion.

    Returns
    -------
    int
        Approximate position for the insertion.

    """

    # Determine SD/SA type.
    splice_type = _get_splice_type(fusion)

    # Determine direction in which to place insertion. Here the orientation
    # and insertion strand is used to first determine the orientation of
    # the gene, and the splice type then determines the actual direction.
    direction = ((-1 if splice_type == 'SD' else 1) *
                 (1 if orientation == 'sense' else -1) * strand)

    # Search for a suitable position, starting from the initial anchor.
    anchor_location = fusion.anchor_genome
    location = anchor_location + (direction * offset)

    # Perform the search, within a specified maximum distance.
    converged = False

    while not converged and abs(location - anchor_location) < max_dist:
        # Check if no features are within 'offset' distance.
        bounds = location, location + ((direction * -1) * (offset - 1))
        start, end = min(bounds), max(bounds)

        overlap = reference_gtf.get_region(fusion.seqname, start, end,
                                           filters=dict(feature='exon'))

        if len(overlap) > 0:
            # Advance our position to avoid the encountered overlap.
            if direction == 1:
                # Move in forward direction.
                location = overlap['end'].max() + offset
            else:
                # Move in reverse direction.
                location = overlap['start'].min() - offset
        else:
            # We have found our position, so stop search.
            converged = True

    if not converged:
        location = np.nan

    return location


def _get_splice_type(fusion):
    """Determines the splice type (SA/SD) of a fusion."""

    feature = fusion.feature_name

    if feature is None:
        raise ValueError('Fusion does not have an annotated transposon feature')
    elif feature.endswith('SA'):
        return 'SA'
    elif feature.endswith('SD'):
        return 'SD'
    else:
        raise ValueError('Unknown splice type for feature {}'.format(feature))


def insertion_orientation(fusion):
    """Identifies the strand and orientation of the corresponding insertion.

    Parameters
    ----------
    fusion : Fusion
        Fusion to identify the strand and orientation for.

    Returns
    -------
    tuple(int, str)
        Strand and orientation of the insertion corresponding
        with the given fusion.

    """

    strand = fusion.strand_genome * fusion.strand_transposon

    if np.isnan(strand) or np.isnan(fusion.gene_strand):
        orientation = None
    elif strand == fusion.gene_strand:
        orientation = 'sense'
    else:
        orientation = 'antisense'

    return strand, orientation
