# -*- coding: utf-8 -*-
"""Contains utility functions for building references."""

# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import re
from typing import Iterable, Tuple

import pathlib2 as pathlib

import pandas as pd
import pyfaidx

from imfusion.util.tabix import GtfIterator

REGION_REGEX = re.compile(r'(?P<chromosome>\w+):(?P<start>\d+)-(?P<end>\d+)')


def build_reference(
        reference_path,  # type: pathlib.Path
        transposon_path,  # type: pathlib.Path
        output_path,  # type: pathlib.Path
        blacklisted_regions=None  # type: Iterable[Tuple[str, int, int]]
):  # type: (...) -> None
    """Builds an augmented reference genome containing the transposon.

    Optionally, specific regions of the original reference can be blacklisted if
    these regions contain sequences that are also in the transposon sequence,
    as these may be problematic in the alignment. For example, in our
    Sleeping Beauty screen we masked regions in the Foxf2 and in En2 genes,
    as the T2onc2 transposon we used contains sequences from these genes.
    Sequences in blacklisted regions are replaced by Ns, effectively removing
    the regions from the reference.

    Parameters
    ----------
    reference_path : pathlib.Path
        Path to the fasta file containing the reference genome.
    transposon_path : pathlib.Path
        Path to the fasta file containg the transposon sequence.
    output_path : pathlib.Path
        Output path for the augmented reference genome (which is written
        as a new fasta file).
    blacklisted_regions : List[(chrom, start, ends)]
        List of regions that should be blacklisted in the augmented reference
        genome. Regions are specified as a tuple of (chr, start, end).
    """

    # Combine reference and transposon sequences.
    _concatenate_fastas(
        [reference_path, transposon_path], output_path=output_path)

    # Mask any blacklisted regions.
    if len(blacklisted_regions) > 0:
        mask_reference(output_path, blacklisted_regions)


def _concatenate_fastas(fasta_paths, output_path):
    # type: (Iterable[pathlib.Path], pathlib.Path) -> None
    """Concatenates multiple fasta files into a single file.

    This function combines multiple fasta files into a single output
    file. It is mainly used to generate the combined reference genome
    that contains both the reference genome sequence and the
    transposon sequence.

    Parameters
    ----------
    fasta_paths : Iterable[pathlib.Path]
        Paths to fasta files that should be concatentated.
    output_path : pathlib.Paths
        Path for the combined output file.

    """

    with output_path.open('wt') as fileobj:
        for fasta_path in fasta_paths:
            # Open fasta file.
            fasta = _open_fasta_file(fasta_path)

            # Read/write sequences.
            for record in fasta:
                _write_fasta_record(
                    fileobj, name=record.long_name, sequence=str(record))


def _open_fasta_file(file_path, **kwargs):
    """Opens a fasta file using pyfaidx."""

    try:
        fasta = pyfaidx.Fasta(str(file_path), **kwargs)
    except pyfaidx.FastaIndexingError:
        raise ValueError('Failed to read fasta file ({}). Ensure that '
                         'the file is a valid fasta file.'.format(
                             str(file_path)))
    return fasta


def _write_fasta_record(fileobj, name, sequence, max_width=80):
    """Writes a fastq record to an opened file."""

    # Group sequence into chunks of max_width.
    chunked = '\n'.join((sequence[i:i + max_width]
                         for i in range(0, len(sequence), max_width)))

    # Write record to file.
    fileobj.write(u'>{}\n{}\n'.format(name, chunked))


def mask_reference(reference_path, blacklist_regions):
    # type: (pathlib.Path, Iterable[Tuple[str, int, int]]) -> None
    """Masks blacklisted regions in a given reference sequence.

    This function removes blacklisted regions from a given reference genome.
    Blacklisted regions are removed by replacing their original sequence
    with a sequences of 'N' nucleotides.

    Parameters
    ----------
    refseq_path : pathlib.Path
        Path to the reference sequence in fasta format.
    blacklist_regions : Iterable[tuple(str, int, int)]:
        List of regions that should be blacklisted in the reference
        sequence. Regions should be specified as a tuple of (chromosome,
        start_position, end_position). For example: ('1', 2000, 2200).

    """

    # Open the reference as a mutable sequence.
    ref = _open_fasta_file(reference_path, mutable=True)

    # Mask regions.
    for seqname, start, end in blacklist_regions:
        ref[seqname][start:end] = 'N' * (end - start)

    # Close file.
    ref.close()


def regions_from_strings(region_strs):
    # type: (Iterable[str]) -> List[Tuple[str, int, int]]
    """Builds blacklist region list for region strings.

    Parses region strings into a list of blacklist region tuples.
    Region strings should be provided in the following format:
    'chromosome:start-end'. For example, 'X:1000-2000' denotes a
    a region on chromosome X from position 1000 to 2000.

    Parameters
    ----------
    region_strs : Iterable[str]
        List of region strings.

    Returns
    -------
    Listp[Tuple(str, int, int]]
        List of blacklist region tuples.

    """

    return [_parse_region_str(s) for s in region_strs]


def _parse_region_str(region_str):
    # type: (str) -> Tuple[str, int, int]

    match = REGION_REGEX.search(region_str)

    if match is None:
        raise ValueError('Unable to parse region {!r}. Ensure that region '
                         'specifications are in the required format '
                         '(chromosome:start-end).'.format(region_str))
    else:
        groups = match.groupdict()
        return groups['chromosome'], int(groups['start']), int(groups['end'])


def regions_from_genes(gene_ids, gtf_path):
    # type: (List[str], pathlib.Path) -> List[Tuple[str, int, int]]
    """Builds blacklist frame for given genes.

    Returns a list of blacklist regions encompassing the regions
    spanned by the genes corresponding to the given gene ids.

    Parameters
    ----------
    gene_ids : List[str]
        List of (Ensembl) gene ids.
    gtf_path : pathlib.Path
        Path to gtf file.

    Returns
    -------
    Tuple[str, int, int]
        List of blacklist region tuples.

    """

    # Fetch genes.
    gtf_iter = GtfIterator(gtf_path)
    gene_records = list(gtf_iter.fetch_genes(gene_ids))

    # Check for missing genes.
    for gene_id, record in zip(gene_ids, gene_records):
        if record is None:
            raise ValueError(
                'Gene {!r} not found in reference gtf. Make sure that the '
                'correct gene ID is being used (corresponding with the '
                'gene_id annotation in the GTF file).'.format(gene_id))

    return [(gene.contig, gene.start, gene.end) for gene in gene_records]


def check_feature_file(feature_path):
    """Checks if a transposon feature frame is valid."""

    transposon_features = pd.read_csv(str(feature_path), sep='\t')

    # Check dataframe columns.
    req_cols = {'name', 'start', 'end', 'strand', 'type'}
    missing_cols = req_cols - set(transposon_features.columns)

    if len(missing_cols) > 0:
        raise ValueError('Missing required columns {!r} in the transposon '
                         'file. Make sure that the file is a properly '
                         'formatted tab-separated file (see documentation for '
                         'more details).'.format(missing_cols))

    # Check if we have any features.
    if len(transposon_features) < 1:
        raise ValueError('Transposon feature file is empty.')

    # Check if we have a valid feature.
    num_valid = transposon_features['type'].isin({'SA', 'SD'}).sum()
    if num_valid < 1:
        raise ValueError('No valid features (with type = SD or SA)'
                         'were found in the transposon feature file.')
