# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import pyfaidx

from imfusion.util.tabix import GtfFile


def concatenate_fastas(fasta_paths, output_path):
    """Concatenates multiple fasta files into a single file.

    This function combines multiple fasta files into a single output
    file. It is mainly used to generate the combined reference genome
    that contains both the reference genome sequence and the
    transposon sequence.

    Parameters
    ----------
    fasta_paths : List[pathlib.Path]
        Paths to fasta files that should be concatentated.
    output_path : pathlib.Paths
        Path for the combined output file.

    """

    with output_path.open('wt') as fh:
        for fasta_path in fasta_paths:
            # Open fasta file.
            fasta = pyfaidx.Fasta(str(fasta_path))

            # Read/write sequences.
            for record in fasta:
                _write_fasta_record(fh, name=record.long_name,
                                    sequence=str(record))


def _write_fasta_record(fh, name, sequence, max_width=80):
    """Writes a fastq record to an opened file."""

    # Group sequence into chunks of max_width.
    chunked = '\n'.join((sequence[i:i+max_width] for i in
                        range(0, len(sequence), max_width)))

    # Write record to file.
    fh.write(u'>{}\n{}\n'.format(name, chunked))


def mask_regions(refseq_path, blacklist_regions):
    """Masks blacklisted regions in a given reference sequence.

    This function removes blacklisted regions from a given reference genome.
    Blacklisted regions are removed by replacing their original sequence
    with a sequences of 'N' nucleotides.

    Parameters
    ----------
    refseq_path : pathlib.Path
        Path to the reference sequence in fasta format.
    blacklist_regions : List[tuple(str, int, int)]:
        List of regions that should be blacklisted in the reference
        sequence. Regions should be specified as a tuple of (chromosome,
        start_position, end_position). For example: ('1', 2000, 2200).

    Returns
    -------
    pathlib.Path
        Path to edited reference sequence.

    """

    # Open the reference as a mutable sequence.
    ref = pyfaidx.Fasta(str(refseq_path), mutable=True)

    # Mask regions.
    for seqname, start, end in blacklist_regions:
        ref[seqname][start:end] = 'N' * (end - start)

    # Close file.
    ref.close()

    return refseq_path


def blacklist_for_regions(region_strs):
    """Builds blacklist region list for region strings.

    Parses region strings into a list of blacklist region tuples.
    Region strings should be provided in the following format:
    'chromosome:start-end'. For example, 'X:1000-2000' denotes a
    a region on chromosome X from position 1000 to 2000.

    Parameters
    ----------
    region_strs : List[str]
        List of region strings.

    Returns
    -------
    tuple(str, int, int)
        List of blacklist region tuples.

    """

    return [_parse_region_str(s) for s in region_strs]


def _parse_region_str(region_str):
    seqname, range_str = region_str.split(':')
    start, end = tuple(int(e) for e in range_str.split('-'))
    return seqname, start, end


def blacklist_for_genes(gene_ids, reference_gtf):
    """Builds blacklist frame for given genes.

    Returns a list of blacklist regions encompassing the regions
    spanned by the genes corresponding to the given gene ids.

    Parameters
    ----------
    gene_ids : List[str]
        List of (Ensembl) gene ids.
    reference_gtf : GtfFile
        GtfFile instance containing reference gene features.

    Returns
    -------
    tuple(str, int, int)
        List of blacklist region tuples.

    """

    return [_get_gene_region(gene_id, reference_gtf)
            for gene_id in gene_ids]


def _get_gene_region(gene_id, reference_gtf):
    gene = reference_gtf.get_gene(gene_id)
    return (gene.contig, gene.start, gene.end)
