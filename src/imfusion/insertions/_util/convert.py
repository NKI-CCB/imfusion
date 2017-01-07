from imfusion.model import Insertion
from imfusion.util.tabix import GtfFile


def convert_to_insertions(fusions, id_format='INS_{}'):
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

    dropped_metadata = {'strand_genome', 'strand_transposon'}

    for i, fusion in enumerate(fusions):
        # Ignore fusions not from a recognised transposon feature.
        if 'feature_name' in fusion.metadata:
            strand = fusion.strand_genome * fusion.strand_transposon

            if 'gene_strand' in fusion.metadata:
                gene_strand = fusion.metadata['gene_strand']
                orientation = 'sense' if strand == gene_strand else 'antisense'
            else:
                orientation = None

            ins_metadata = {
                k: v
                for k, v in fusion.metadata.items()
                if k not in dropped_metadata
            }
            ins_metadata['transposon_anchor'] = fusion.anchor_transposon

            if orientation is not None:
                ins_metadata['orientation'] = orientation

            yield Insertion(
                id=id_format.format(i + 1),
                seqname=fusion.seqname,
                position=fusion.anchor_genome,
                strand=strand,
                support=fusion.support,
                metadata=ins_metadata)
