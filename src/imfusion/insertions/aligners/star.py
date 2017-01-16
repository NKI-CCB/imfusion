# pylint: disable=wildcard-import,redefined-builtin,unused-wildcard-import
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=wildcard-import,redefined-builtin,unused-wildcard-import

import collections
import operator
import itertools
import sys
from typing import Any
import re

from intervaltree import IntervalTree
import numpy as np
import pandas as pd
import toolz

from imfusion.build.indexers.star import StarReference
from imfusion.model import Fusion, TransposonFusion
from imfusion.util import shell, tabix, path

from .base import Aligner, register_aligner
from .. import util

CIGAR_MATCH_REGEX = re.compile(r'(\d+)M')


class StarAligner(Aligner):
    """STAR aligner."""

    def __init__(
            self,
            reference,  # type: StarReference
            assemble=False,  # type: bool
            assemble_args=None,  # type: Dict[str, Any]
            min_flank=12,  # type: int
            threads=1,  # type: int
            external_sort=False,  # type: bool
            extra_args=None,  # type: Dict[str, Any]
            logger=None,  # type: Any
            merge_junction_dist=10,  # type: int
            max_spanning_dist=300,  # type: int
            max_junction_dist=10000,  # type: int
            filter_features=True,  # type: bool
            filter_orientation=True,  # type: bool
            filter_blacklist=None  # type: List[str]
    ):  # type: (...) -> None

        super().__init__(reference=reference, logger=logger)

        self._assemble = assemble
        self._assemble_args = assemble_args or {}
        self._min_flank = min_flank
        self._threads = threads
        self._external_sort = external_sort
        self._extra_args = extra_args or {}

        self._merge_junction_dist = merge_junction_dist
        self._max_spanning_dist = max_spanning_dist
        self._max_junction_dist = max_junction_dist

        self._filter_features = filter_features
        self._filter_orientation = filter_orientation
        self._filter_blacklist = filter_blacklist

    @property
    def dependencies(self):
        programs = ['STAR']

        if self._external_sort:
            if shell.which('sambamba') is None:
                # Fall back to samtools if sambamba is not available.
                programs += ['samtools']
            else:
                programs += ['sambamba']

        if self._assemble:
            programs += ['stringtie']

        return programs

    def identify_insertions(self, fastq_path, output_dir, fastq2_path=None):
        """Identifies insertions from given reads."""

        output_dir.mkdir(parents=True, exist_ok=True)

        # Perform alignment using STAR.
        alignment_path = output_dir / 'alignment.bam'
        if not alignment_path.exists():
            self._logger.info('Performing alignment')
            self._align(fastq_path, output_dir, fastq2_path=fastq2_path)
        else:
            self._logger.info('Using existing alignment')

        # Assemble transcripts if requested.
        if self._assemble:
            assembled_path = output_dir / 'assembled.gtf.gz'
            if not assembled_path.exists():
                self._logger.info('Assembling transcripts')

                # Generate assembled GTF.
                stringtie_out_path = assembled_path.with_suffix('')
                util.stringtie_assemble(
                    alignment_path,
                    gtf_path=self._reference.gtf_path,
                    output_path=stringtie_out_path)

                # Compress and index.
                tabix.index_gtf(stringtie_out_path, output_path=assembled_path)
                stringtie_out_path.unlink()
            else:
                self._logger.info('Using existing assembly')
        else:
            assembled_path = None

        # Extract identified fusions.
        self._logger.info('Extracting fusions')
        fusion_path = output_dir / 'fusions.out'
        fusions = list(self._extract_fusions(fusion_path))

        # Extract insertions.
        self._logger.info('Extracting insertions')
        insertions = list(
            util.extract_insertions(
                fusions,
                gtf_path=self._reference.indexed_gtf_path,
                features_path=self._reference.features_path,
                assembled_gtf_path=assembled_path,
                ffpm_fastq_path=fastq_path,
                chromosomes=None))

        self._logger.info('Filtering insertions')
        insertions = util.filter_insertions(
            insertions,
            features=self._filter_features,
            orientation=self._filter_orientation,
            blacklist=self._filter_blacklist)

        for insertion in insertions:
            yield insertion

    def _align(self, fastq_path, output_dir, fastq2_path=None):
        # Put output into subdirectory, we will symlink the
        # expected outputs from STAR later.
        star_dir = output_dir / '_star'

        # Gather default arguments.
        sort_type = 'Unsorted' if self._external_sort else 'SortedByCoordinate'

        args = {
            '--twopassMode': ('Basic', ),
            '--outReadsUnmapped': ('None', ),
            '--outSAMtype': ('BAM', sort_type),
            '--runThreadN': (self._threads, ),
            '--chimSegmentMin': (self._min_flank, ),
            '--outSAMstrandField': ('intronMotif', )  # XS field for Stringtie.
        }

        star_align(
            fastq_path=fastq_path,
            fastq2_path=fastq2_path,
            index_path=self._reference.index_path,
            output_dir=star_dir,
            extra_args=toolz.merge(args, self._extra_args),
            stdout=sys.stderr,
            logger=self._logger)

        # If not yet sorted, sort bam file using samtools/sambamba.
        sorted_bam_path = star_dir / 'Aligned.sortedByCoord.out.bam'

        if sort_type == 'Unsorted':
            unsorted_bam_path = star_dir / 'Aligned.out.bam'
            _sort_bam(
                unsorted_bam_path, sorted_bam_path, threads=self._threads)
            unsorted_bam_path.unlink()

        # Symlink fusions and alignment into expected location.
        dest_bam_path = output_dir / 'alignment.bam'
        path.symlink_relative(
            src_path=sorted_bam_path, dest_path=dest_bam_path)

        path.symlink_relative(
            src_path=star_dir / 'Chimeric.out.junction',
            dest_path=output_dir / 'fusions.out')

    def _extract_fusions(self, fusion_path):
        # Read chimeric junction data.
        chimeric_data = read_chimeric_junctions(fusion_path)

        # Extract transposon fusions.
        fusions = extract_transposon_fusions(
            chimeric_data,
            self._reference.transposon_name,
            merge_junction_dist=self._merge_junction_dist,
            max_spanning_dist=self._max_spanning_dist,
            max_junction_dist=self._max_junction_dist)

        for fusion in fusions:
            yield fusion

    @classmethod
    def configure_args(cls, parser):
        super().configure_args(parser)

        star_group = parser.add_argument_group('STAR arguments')
        star_group.add_argument('--star_threads', type=int, default=1)
        star_group.add_argument('--star_min_flank', type=int, default=12)
        star_group.add_argument(
            '--star_external_sort', default=False, action='store_true')
        star_group.add_argument(
            '--star_args', type=shell.parse_arguments, default='')

        star_group.add_argument('--merge_junction_dist', default=10, type=int)
        star_group.add_argument('--max_spanning_dist', default=300, type=int)
        star_group.add_argument('--max_junction_dist', default=10000, type=int)

        assemble_group = parser.add_argument_group('Assembly')
        assemble_group.add_argument(
            '--assemble', default=False, action='store_true')

        filt_group = parser.add_argument_group('Filtering')
        filt_group.add_argument(
            '--no_filter_orientation',
            dest='filter_orientation',
            default=True,
            action='store_false')
        filt_group.add_argument(
            '--no_filter_feature',
            dest='filter_features',
            default=True,
            action='store_false')
        filt_group.add_argument('--blacklisted_genes', nargs='+')

    @classmethod
    def parse_args(cls, args):
        kws = dict(
            reference=StarReference(args.reference),
            min_flank=args.star_min_flank,
            threads=args.star_threads,
            extra_args=args.star_args,
            external_sort=args.star_external_sort,
            assemble=args.assemble,
            merge_junction_dist=args.merge_junction_dist,
            max_spanning_dist=args.max_spanning_dist,
            max_junction_dist=args.max_junction_dist,
            filter_features=args.filter_features,
            filter_orientation=args.filter_orientation,
            filter_blacklist=args.blacklisted_genes)

        return kws


register_aligner('star', StarAligner)


def _sort_bam(input_bam, output_bam, threads=1):
    if shell.which('sambamba') is not None:
        tmp_dir = output_bam.parent / '_tmp'
        args = [
            'sambamba', 'sort', '-o', str(output_bam),
            '--tmpdir={}'.format(str(tmp_dir)), '-t', str(threads),
            str(input_bam)
        ]
    elif shell.which('samtools') is not None:
        args = ['samtools', 'sort', '-o', str(output_bam), str(input_bam)]
    else:
        raise ValueError('Samtools or Sambamba must be installed for '
                         'external sorting')

    shell.run_command(args=args)


def build_star_index(reference_path,
                     gtf_path,
                     output_base_path,
                     overhang=100,
                     stdout=None,
                     stderr=None):
    """Runs STAR in build mode to generate a reference index."""

    if not output_base_path.exists():
        output_base_path.mkdir(parents=True)

    cmdline_args = [
        'STAR', '--runMode', 'genomeGenerate', '--genomeDir',
        str(output_base_path), '--genomeFastaFiles', str(reference_path),
        '--sjdbGTFfile', str(gtf_path), '--sjdbOverhang', str(overhang)
    ]

    shell.run_command(args=cmdline_args, stdout=stdout, stderr=stderr)


def star_align(fastq_path,
               index_path,
               output_dir,
               fastq2_path=None,
               extra_args=None,
               stdout=None,
               stderr=None,
               logger=None):
    """Runs STAR in alignment mode for the given fastqs."""

    extra_args = extra_args or {}

    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # Check if files are gzipped.
    if fastq_path.suffixes[-1] == '.gz':
        extra_args['--readFilesCommand'] = ('gunzip', '-c')

    # Assemble arguments.
    if fastq2_path is None:
        fastq_args = [str(fastq_path)]
    else:
        fastq_args = [str(fastq_path), str(fastq2_path)]

    cmdline_args = [
        'STAR', '--genomeDir', str(index_path), '--outFileNamePrefix',
        str(output_dir) + '/', '--readFilesIn'
    ] + fastq_args
    cmdline_args += [
        str(arg) for arg in shell.flatten_arguments(extra_args or {})
    ]

    shell.run_command(
        args=cmdline_args, stdout=stdout, stderr=stderr, logger=logger)


def read_chimeric_junctions(chimeric_path):
    """Reads junctions from STARs Chimeric.out.junction output file.

    Parameters
    ----------
    pathlib.Path
        Path to the Chimeric.out.junction file.

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame containing the junction data.

    """

    # File structure is as follows.
    #
    # http://labshare.cshl.edu/shares/gingeraslab/www-data/
    # dobin/STAR/STAR.posix/doc/STARmanual.pdf
    #
    # The first 9 columns give information about the chimeric junction:
    #
    #   column 1: chromosome of the donor
    #   column 2: first base of the intron of the donor (1-based)
    #   column 3: strand of the donor
    #   column 4: chromosome of the acceptor
    #   column 5: first base of the intron of the acceptor (1-based)
    #   column 6: strand of the acceptor
    #   column 7: junction type: -1=encompassing junction
    #       (between the mates), 1=GT/AG, 2=CT/AC
    #   column 8: repeat length to the left of the junction
    #   column 9: repeat length to the right of the junction
    #
    # Columns 10-14 describe the alignments of the two chimeric segments,
    # it is SAM like. Alignments are given with respect to the (+) strand
    #
    #   column 10: read name
    #   column 11: first base of the first segment (on the + strand)
    #   column 12: CIGAR of the first segment
    #   column 13: first base of the second segment
    #   column 14: CIGAR of the second segment
    #

    names = [
        'seqname_a', 'location_a', 'strand_a', 'seqname_b', 'location_b',
        'strand_b', 'junction_type', 'repeat_length_left',
        'repeat_length_right', 'read_name', 'first_segment_base',
        'first_segment_cigar', 'second_segment_base', 'second_segment_cigar'
    ]

    junctions = pd.read_csv(
        str(chimeric_path),
        sep='\t',
        header=None,
        names=names,
        dtype={'seqname_a': str,
               'seqname_b': str})

    for strand_col in ['strand_a', 'strand_b']:
        junctions[strand_col] = junctions[strand_col].map({'+': 1, '-': -1})

    return junctions


def normalize_chimeric_junctions(chimeric_data, seqname=None):
    """Normalizes chimeric junction data so that seqnames are ordered.

    This function normalizes a chimeric junction DataFrame so that seqname_a
    and seqname_b are lexically ordered (i.e., seqname_a < seqname_b).
    Intra-chromosomal junctions (with seqname_a == seqname_b) are normalized
    using their location, so that location_a < location_b.

    Alternatively, if the seqname argument is given, the junctions are
    normalized so that seqname_a == seqname. Any junctions not including
    the seqname chromosome are dropped.

    Parameters
    ----------
    chimeric_data : pd.DataFrame
        DataFrame of chimeric junctions to normalize.
    seqname : str, optional
        Optional sequence name to normalize for. If given, junctions are
        normalized so that seqname_a == seqname. Any junctions not including
        the named sequence are dropped.

    Returns
    -------
    pd.DataFrame
        Normalized chimeric junction DataFrame.

    """

    if seqname is None:
        # Reads are 'normal' if seq_a < seq_b or, in the case that
        # seq_a == seq_b, if loc_a < loc_b.
        same_seq = chimeric_data.seqname_a == chimeric_data.seqname_b
        seq_ordered = (chimeric_data.seqname_a < chimeric_data.seqname_b)
        loc_ordered = (chimeric_data.location_a < chimeric_data.location_b)

        norm_mask = (seq_ordered | (same_seq & loc_ordered))
    else:
        # Subset to samples that have at least one end on seqname.
        chimeric_data = chimeric_data.ix[(chimeric_data.seqname_a == seqname) |
                                         (chimeric_data.seqname_b == seqname)]

        # Chimeric reads are 'normal' if first side is on seqname and
        # and if locations are ordered in the case that seq_a == seq_b.
        loc_ordered = (chimeric_data.location_a < chimeric_data.location_b)
        norm_mask = (chimeric_data.seqname_a == seqname) & loc_ordered

    normed = pd.concat(
        [
            chimeric_data.ix[norm_mask],
            _reverse_chimeric(chimeric_data.ix[~norm_mask])
        ],
        axis=0)

    return normed


def _reverse_chimeric(chimeric_data):
    return (chimeric_data.rename(columns={
        'seqname_a': 'seqname_b',
        'location_a': 'location_b',
        'strand_a': 'strand_b',
        'seqname_b': 'seqname_a',
        'location_b': 'location_a',
        'strand_b': 'strand_a',
        'flank_a': 'flank_b',
        'flank_b': 'flank_a'
    }).assign(
        strand_a=lambda df: df.strand_a * -1,
        strand_b=lambda df: df.strand_b * -1).get(chimeric_data.columns))


def extract_transposon_fusions(chimeric_data,
                               transposon_name,
                               merge_junction_dist=10,
                               max_spanning_dist=300,
                               max_junction_dist=10000):
    """Extracts transposon fusions from a STAR chimeric read dataframe."""

    # Subset and normalize chimeric reads for transposon.
    chimeric_data = chimeric_data.ix[
        (chimeric_data['seqname_a'] == transposon_name) ^
        (chimeric_data['seqname_b'] == transposon_name)] # yapf: disable

    chimeric_data = normalize_chimeric_junctions(
        chimeric_data, seqname=transposon_name)

    # Extract junction fusions and merge close junctions.
    junctions = list(
        extract_junction_fusions(
            chimeric_data, merge_dist=merge_junction_dist))

    # Assign spanning reads to junctions.
    junctions, unassigned = assign_spanning_reads(
        junctions,
        chimeric_data,
        max_dist_left=max_spanning_dist,
        max_dist_right=max_junction_dist)

    # Extract spanning fusions from unused reads.
    spanning = extract_spanning_fusions(unassigned, max_dist=max_spanning_dist)
    fusions = itertools.chain.from_iterable([junctions, spanning])

    # Convert to transposon fusions.
    for fusion in fusions:
        yield TransposonFusion.from_fusion(fusion, transposon_name)


def extract_junction_fusions(chimeric_data, merge_dist=None):
    """Extracts junction fusions from a STAR chimeric read dataframe."""

    # Ensure chimeric data only contains junction reads.
    chimeric_data = chimeric_data.query('junction_type >= 0')

    if len(chimeric_data) > 0:
        # Add flanking annotation.
        chimeric_data = chimeric_data.copy()

        flanks = (_flank_sizes(t) for t in chimeric_data.itertuples())
        chimeric_data['flank_a'], chimeric_data['flank_b'] = zip(*flanks)

        # Group by position and summarize.
        grouped = chimeric_data.groupby([
            'seqname_a', 'location_a', 'strand_a', 'seqname_b', 'location_b',
            'strand_b'
        ])

        summarized = (grouped.agg({
            'flank_a': 'max',
            'flank_b': 'max',
            'read_name': 'nunique'
        }).reset_index().assign(support_spanning=0)
                      .rename(columns={'read_name': 'support_junction'}))

        # Transform to Fusions.
        fusions = (Fusion(**toolz.keyfilter(lambda k: k not in {'Index'},
                                            row._asdict()))
                   for row in summarized.itertuples())

        # Merge fusions within dist.
        if merge_dist is not None:
            fusions = Fusion.merge(fusions, max_dist=merge_dist)

        for fusion in fusions:
            yield fusion


def _flank_sizes(row):
    return (_flank_size(row, donor=True), _flank_size(row, donor=False))


def _flank_size(row, donor=True):
    """Calculates flank size of a spanning read on one side of the fusion."""

    # TODO: account for deletions/insertions etc.

    if donor:
        mate_index = -1 if row.strand_a == 1 else 0
        cigar = row.first_segment_cigar
    else:
        mate_index = 0 if row.strand_b == 1 else -1
        cigar = row.second_segment_cigar

    cigar_split = cigar.split('p')
    size = _sum_matches(cigar_split[mate_index])

    return size


def _sum_matches(cigar_str):
    """Sums the number of matches in a CIGAR string."""
    matches = CIGAR_MATCH_REGEX.findall(cigar_str)
    return sum(map(int, matches))


def assign_spanning_reads(junctions, chimeric_data, max_dist_left,
                          max_dist_right):
    """Assigns spanning reads to the closest junction that they support."""

    # Ensure chimeric data only contains spanning reads.
    chimeric_data = chimeric_data.query('junction_type < 0')

    # Assign chimeric reads to junctions.
    trees_a = _build_lookup_trees(junctions, side='a')
    trees_b = _build_lookup_trees(junctions, side='b')

    assignments = collections.defaultdict(list)
    unassigned_idx = []

    for spanning_read in chimeric_data.itertuples():
        junc = _lookup_closest(spanning_read, trees_a, trees_b, max_dist_left,
                               max_dist_right)

        if junc is not None:
            assignments[junc].append(spanning_read)
        else:
            unassigned_idx.append(spanning_read.Index)

    # Augment junctions.
    new_juncs = [
        junc._replace(
            support_spanning=junc.support_spanning + len(assignments[junc]))
        for junc in junctions
    ]

    # Select unassigned reads.
    unassigned = chimeric_data.ix[unassigned_idx]

    return new_juncs, unassigned


def _build_lookup_trees(fusions, side='a'):
    # Define attr getters.
    get_seq = operator.attrgetter('seqname_' + side)
    get_loc = operator.attrgetter('location_' + side)
    get_strand = operator.attrgetter('strand_' + side)

    # Setup keyfunc.
    def _keyfunc(fusion):
        return (get_seq(fusion), get_strand(fusion))

    # Sort and group fusions.
    fusions = sorted(fusions, key=_keyfunc)
    grouped = itertools.groupby(fusions, key=_keyfunc)

    trees = {}
    for key, grp in grouped:
        trees[key] = IntervalTree.from_tuples(
            (get_loc(fus), get_loc(fus) + 1, fus) for fus in grp)

    return trees


def _lookup_closest(fusion,
                    trees_a,
                    trees_b,
                    max_dist_left,
                    max_dist_right,
                    slack=5):
    def _dist(fusion, hit):
        return (abs(fusion.location_a - hit.location_a) +
                abs(fusion.location_b - hit.location_b))

    # Identify candidate fusions.
    overlap_a = _lookup_tree(
        fusion, trees_a, side='a', max_dist=max_dist_left, slack=slack)
    overlap_b = _lookup_tree(
        fusion, trees_b, side='b', max_dist=max_dist_right, slack=slack)
    overlap = list(overlap_a & overlap_b)

    if len(overlap) > 0:
        distances = np.array(_dist(fusion, hit) for hit in overlap)
        return overlap[distances.argmin()]
    else:
        return None


def _lookup_tree(fusion, trees, side, max_dist=300, slack=5):

    # Define attr getters.
    get_seq = operator.attrgetter('seqname_' + side)
    get_loc = operator.attrgetter('location_' + side)
    get_strand = operator.attrgetter('strand_' + side)

    # Setup tree, query range.
    strand = get_strand(fusion)

    try:
        tree = trees[(get_seq(fusion), strand)]

        direction = (-1 if side == 'a' else 1) * strand

        loc = get_loc(fusion)
        if direction == 1:
            start, end = loc - max_dist, loc + slack
        else:
            start, end = loc - slack, loc + max_dist

        # Determine overlap and return fusions.
        overlap = set(tup[2] for tup in tree[start:end])

        return overlap
    except KeyError:
        return set()


def extract_spanning_fusions(chimeric_data, max_dist):
    """Extracts spanning fusions from a STAR chimeric read dataframe."""

    # Select spanning fusions.
    chimeric_data = chimeric_data.query('junction_type < 0')

    # Group by position and summarize.
    for grp in _groupby_position(chimeric_data, max_dist=max_dist):
        first = grp.iloc[0]

        agg_donor = max if first.strand_a == 1 else min
        agg_acc = min if first.strand_b == 1 else max

        yield Fusion(
            seqname_a=first.seqname_a,
            location_a=agg_donor(grp.location_a),
            strand_a=first.strand_a,
            seqname_b=first.seqname_b,
            location_b=agg_acc(grp.location_b),
            strand_b=first.strand_b,
            support_junction=0,
            support_spanning=len(grp),
            flank_a=0,
            flank_b=0)


def _groupby_position(spanning_data, max_dist):
    """Groups reads within given distance at both locations."""

    grp_cols = ['seqname_a', 'seqname_b', 'strand_a', 'strand_b']

    for _, grp in spanning_data.groupby(grp_cols):
        if len(grp) == 1:
            yield grp
        else:
            clust_a = _cluster_positions(
                grp, col='location_a', max_dist=max_dist)
            clust_b = _cluster_positions(
                grp, col='location_b', max_dist=max_dist)

            for _, loc_grp in grp.groupby([clust_a, clust_b]):
                yield loc_grp


def _cluster_positions(df, col, max_dist):
    locs = df[col].sort_values()

    diff = np.diff(np.insert(locs.values, 0, 0))
    groups = pd.Series(np.cumsum(diff > max_dist), index=locs.index)

    return groups.ix[df.index]
