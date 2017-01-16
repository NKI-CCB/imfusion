import argparse
import shutil

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import pytest

from frozendict import frozendict
import pandas as pd

from imfusion.insertions.aligners import star
from imfusion.model import Fusion, TransposonFusion, Insertion

# pylint: disable=no-self-use,redefined-outer-name

# TODO: add tests testing distances.


@pytest.fixture
def chimeric_junctions_path():
    """Returns path to example junction data."""
    return pytest.helpers.data_path(
        'Chimeric.out.junction', relative_to=__file__)


class TestReadChimericJunctions(object):
    """Tests for read_chimeric_junctions."""

    def test_example(self, chimeric_junctions_path):
        """Tests example case."""

        df = star.read_chimeric_junctions(chimeric_junctions_path)

        # Check shape.
        assert df.shape == (505, 14)

        # Check first row.
        first = df.iloc[0]
        assert first.seqname_a == 'T2onc'
        assert first.location_a == 1541
        assert first.strand_a == 1
        assert first.seqname_b == '16'
        assert first.location_b == 52141095
        assert first.strand_b == -1
        assert first.junction_type == 2
        assert first.repeat_length_left == 2
        assert first.repeat_length_right == 3
        assert first.read_name == 'ERR578333.12604936'
        assert first.first_segment_base == 1509
        assert first.first_segment_cigar == '32M68S'
        assert first.second_segment_base == 52139621
        assert first.second_segment_cigar == '8M1378N90M-70p68M32S'


@pytest.fixture
def chimeric_data():
    """Example containing spanning + junction reads from single fusion."""
    return _build_chimeric_data(
        [('1', 300, 1, 'T2onc', 420, 1, 2, '100M2208p38M62S', '62M38S', 'R1'),
         ('1', 300, 1, 'T2onc', 420, 1, 1, '100M2208p52M48S', '48M52S', 'R2'),
         ('1', 301, 1, 'T2onc', 420, 1, 1, '100M2208p52M48S', '48M52S', 'R3'),
         ('1', 300, 1, 'T2onc', 421, 1, 1, '100M2208p52M48S', '48M52S', 'R4'),
         ('1', 280, 1, 'T2onc', 435, 1, -1, '100M', '97M3S', 'S1'),
         ('1', 270, 1, 'T2onc', 445, 1, -1, '100M', '98M2S', 'S2'),
         ('1', 275, 1, 'T2onc', 435, 1, -1, '100M', '98M2S', 'S3')])


def _build_chimeric_data(rows):
    """Helper function to build a fusion dataframe from list of row data."""
    columns = ['seqname_a', 'location_a', 'strand_a',
               'seqname_b', 'location_b', 'strand_b',
               'junction_type', 'first_segment_cigar',
               'second_segment_cigar', 'read_name'] # yapf: disable
    return pd.DataFrame.from_records(rows, columns=columns)


class TestExtractTransposonFusions(object):
    """Tests for test_extract_transposon_fusions function."""

    def test_example(self, chimeric_data):
        """Test simple example."""

        fusions = list(star.extract_transposon_fusions(chimeric_data, 'T2onc'))

        assert fusions == [
            TransposonFusion(
                seqname='1',
                anchor_genome=300,
                anchor_transposon=420,
                strand_genome=-1,
                strand_transposon=-1,
                flank_genome=-62,
                flank_transposon=100,
                support_junction=4,
                support_spanning=3,
                metadata=frozendict({}))
        ]

    def test_example_file(self, chimeric_junctions_path):
        """Test data read from example file."""

        # Extract fusions.
        chimeric_data = star.read_chimeric_junctions(chimeric_junctions_path)
        fusions = set(star.extract_transposon_fusions(chimeric_data, 'T2onc'))

        # Check number of fusions.
        assert len(fusions) == 7

        # Check key example.
        expected = TransposonFusion(
            seqname='16',
            anchor_genome=52141095,
            anchor_transposon=1541,
            strand_genome=-1,
            strand_transposon=1,
            flank_genome=-78,
            flank_transposon=-76,
            support_junction=380,
            support_spanning=118,
            metadata=frozendict({}))
        assert expected in fusions


class TestExtractJunctionFusions(object):
    """Tests for extract_junction_fusions."""

    def test_example(self, chimeric_data):
        """Test simple example."""

        fusions = list(
            star.extract_junction_fusions(
                chimeric_data, merge_dist=10))

        assert fusions == [
            Fusion(
                seqname_a='1',
                location_a=300,
                strand_a=1,
                seqname_b='T2onc',
                location_b=420,
                strand_b=1,
                flank_a=52,
                flank_b=62,
                support_junction=4,
                support_spanning=0)
        ]


class TestExtractSpanningFusions(object):
    """Tests for extract_spanning_fusions."""

    def test_example(self, chimeric_data):
        """Test simple example."""

        fusions = list(
            star.extract_spanning_fusions(
                chimeric_data, max_dist=300))

        assert fusions == [
            Fusion(
                seqname_a='1',
                location_a=280,
                strand_a=1,
                seqname_b='T2onc',
                location_b=435,
                strand_b=1,
                flank_a=0,
                flank_b=0,
                support_junction=0,
                support_spanning=3)
        ]


@pytest.fixture
def star_align_kws(tmpdir):
    "Basic arguments for star_align." ""
    return {
        'fastq_path': Path('in.R1.fastq.gz'),
        'fastq2_path': Path('in.R2.fastq.gz'),
        'index_path': Path('/path/to/index'),
        'output_dir': Path(str(tmpdir / '_star')),
        'extra_args': None
    }


class TestStarAlign(object):
    """Tests for star_align function."""

    def test_basic(self, star_align_kws, mocker):
        """Tests basic call."""

        mock_align = mocker.patch.object(star.shell, 'run_command')
        star.star_align(**star_align_kws)

        # Check if output directory was created.
        assert star_align_kws['output_dir'].exists()

        # Check arguments to run command. Checks if we have an extra slash
        # added to output dir, if the fastqs were added properly, and if
        # readFilesCommand is added for the gzipped files.
        mock_align.assert_called_once_with(
            args=[
                'STAR', '--genomeDir', '/path/to/index', '--outFileNamePrefix',
                str(star_align_kws['output_dir']) + '/', '--readFilesIn',
                str(star_align_kws['fastq_path']),
                str(star_align_kws['fastq2_path']), '--readFilesCommand',
                'gunzip', '-c'
            ],
            stdout=None,
            stderr=None,
            logger=pytest.helpers.mock_any(object))

    def test_extra_args(self, star_align_kws, mocker):
        """Tests call with extra args."""

        star_align_kws['extra_args'] = {
            '--alignIntronMax': (200000, ),
            '--twopassMode': ('None', )
        }

        mock_align = mocker.patch.object(star.shell, 'run_command')
        star.star_align(**star_align_kws)

        args = pytest.helpers.mock_call(mock_align)[1]['args']

        # Check new argument.
        assert '--alignIntronMax' in args
        assert args[args.index('--alignIntronMax') + 1] == '200000'

        # Check overridden argument.
        assert '--twopassMode' in args
        assert args[args.index('--twopassMode') + 1] == 'None'

    def test_unpaired(self, star_align_kws, mocker):
        """Tests call for unpaired samples."""

        star_align_kws['fastq2_path'] = None

        mock_align = mocker.patch.object(star.shell, 'run_command')
        star.star_align(**star_align_kws)

        args = pytest.helpers.mock_call(mock_align)[1]['args']
        assert '--readFilesIn' in args
        assert args[args.index('--readFilesIn') + 1] == str(star_align_kws[
            'fastq_path'])

    def test_no_gzip(self, star_align_kws, mocker):
        """Tests call for non-gzipped file."""

        star_align_kws['fastq_path'] = Path('in1.fastq')
        star_align_kws['fastq2_path'] = Path('in2.fastq')

        mock_align = mocker.patch.object(star.shell, 'run_command')
        star.star_align(**star_align_kws)

        args = pytest.helpers.mock_call(mock_align)[1]['args']
        assert '--readFilesCommand' not in args


@pytest.fixture
def read_paths():
    """Example read paths."""
    return (Path('a.fastq.gz'), Path('b.fastq.gz'))


@pytest.fixture
def star_reference():
    """Returns example reference for STAR."""
    ref_path = pytest.helpers.data_path('reference', relative_to=__file__)
    return star.StarReference(ref_path)


@pytest.fixture
def star_output_dir(tmpdir, chimeric_junctions_path):
    """Simulated star output directory."""
    # Create directories.
    output_dir = Path(str(tmpdir / 'out'))

    star_dir = output_dir / '_star'
    star_dir.mkdir(parents=True)

    # Copy / simulate aligner output files.
    shutil.copy(str(chimeric_junctions_path),
                str(star_dir / 'Chimeric.out.junction'))  # yapf: disable

    pytest.helpers.touch(star_dir / 'Aligned.sortedByCoord.out.bam')

    return output_dir


@pytest.fixture
def cmdline_args(tmpdir):
    """Example command line arguments."""
    return [
        '--fastq', 'a.fastq.gz', '--fastq2', 'b.fastq.gz', '--reference',
        str(tmpdir), '--output_dir', '/path/to/out'
    ]


class TestStarAligner(object):
    """Tests for the StarAligner class."""

    def test_dependencies(self, star_reference, mocker):
        """Test dependencies for various configurations."""

        # Basic call.
        aligner = star.StarAligner(star_reference)
        assert aligner.dependencies == ['STAR']

        # Using external sort.
        aligner2 = star.StarAligner(star_reference, external_sort=True)
        assert aligner2.dependencies == ['STAR', 'samtools']

        # Using assembly.
        aligner3 = star.StarAligner(star_reference, assemble=True)
        assert aligner3.dependencies == ['STAR', 'stringtie']

        # Using external sort with sambamba installed.
        mocker.patch.object(
            star.shell, 'which', return_value='/usr/bin/sambamba')
        aligner4 = star.StarAligner(star_reference, external_sort=True)
        assert aligner4.dependencies == ['STAR', 'sambamba']

    def test_identify_insertions(self, read_paths, star_reference,
                                 star_output_dir, mocker):
        """Integration test using data from an example run."""

        # TODO: Check assembly?

        # Mock functions call.
        star_mock = mocker.patch.object(star, 'star_align')
        mocker.patch.object(star.util, 'count_lines', return_value=8e6)

        # Call identify insertions.
        fastq, fastq2 = read_paths

        aligner = star.StarAligner(star_reference)
        ins = list(
            aligner.identify_insertions(
                fastq, star_output_dir, fastq2_path=fastq2))

        # Check call to star_align.
        star_mock.assert_called_once_with(
            fastq_path=fastq,
            fastq2_path=fastq2,
            output_dir=star_output_dir / '_star',
            index_path=star_reference.index_path,
            extra_args={
                '--twopassMode': ('Basic', ),
                '--outReadsUnmapped': ('None', ),
                '--outSAMtype': ('BAM', 'SortedByCoordinate'),
                '--runThreadN': (1, ),
                '--chimSegmentMin': (12, ),
                '--outSAMstrandField': ('intronMotif', )
            },
            stdout=pytest.helpers.mock_any(object),
            logger=pytest.helpers.mock_any(object))

        # Check result, including specific Cblb insertion.
        assert len(ins) == 5

        assert ins[2].id == 'INS_3'
        assert ins[2].seqname == '16'
        assert ins[2].position == 52141095
        assert ins[2].strand == -1
        assert ins[2].support_junction == 380
        assert ins[2].support_spanning == 118
        assert ins[2].support == 498
        assert ins[2].metadata['gene_id'] == 'ENSMUSG00000022637'
        assert ins[2].metadata['transposon_anchor'] == 1541
        assert ins[2].metadata['feature_name'] == 'En2SA'
        assert ins[2].metadata['gene_name'] == 'Cblb'
        assert ins[2].metadata['feature_type'] == 'SA'
        assert ins[2].metadata['gene_strand'] == 1
        assert ins[2].metadata['orientation'] == 'antisense'
        assert ins[2].metadata['ffpm_junction'] == 190.0
        assert ins[2].metadata['ffpm_spanning'] == 59.0
        assert ins[2].metadata['ffpm'] == 249.0

    def test_from_args_basic(self, cmdline_args):
        """Tests creation with minimal arguments."""

        # Setup parser.
        parser = argparse.ArgumentParser()
        star.StarAligner.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        aligner = star.StarAligner.from_args(args)

        # Check args.
        assert args.fastq == Path('a.fastq.gz')
        assert args.fastq2 == Path('b.fastq.gz')
        assert args.output_dir == Path('/path/to/out')

        # Check aligner.
        assert not aligner._assemble
        assert aligner._assemble_args == {}
        assert aligner._min_flank == 12
        assert aligner._threads == 1
        assert not aligner._external_sort
        assert aligner._extra_args == {}
        assert aligner._merge_junction_dist == 10
        assert aligner._max_spanning_dist == 300
        assert aligner._max_junction_dist == 10000
        assert aligner._filter_features
        assert aligner._filter_orientation
        assert aligner._filter_blacklist is None

    def test_from_args_extra_args(self, cmdline_args):
        """Tests creation with extra options."""

        cmdline_args += [
            '--star_threads', '5', '--star_min_flank', '20',
            '--star_external_sort', '--star_args', "--limitBAMsortRAM 2000",
            '--merge_junction_dist', '20', '--max_spanning_dist', '600',
            '--max_junction_dist', '50000', '--assemble',
            '--no_filter_orientation', '--no_filter_feature',
            '--blacklisted_genes', 'En2'
        ] # yapf:disable

        # Setup parser.
        parser = argparse.ArgumentParser()
        star.StarAligner.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        aligner = star.StarAligner.from_args(args)

        # Check args.
        assert args.fastq == Path('a.fastq.gz')
        assert args.fastq2 == Path('b.fastq.gz')
        assert args.output_dir == Path('/path/to/out')

        # Check aligner.
        assert aligner._assemble
        assert aligner._assemble_args == {}
        assert aligner._min_flank == 20
        assert aligner._threads == 5
        assert aligner._external_sort
        assert aligner._extra_args == {'--limitBAMsortRAM': ('2000', )}
        assert aligner._merge_junction_dist == 20
        assert aligner._max_spanning_dist == 600
        assert aligner._max_junction_dist == 50000
        assert not aligner._filter_features
        assert not aligner._filter_orientation
        assert aligner._filter_blacklist == ['En2']
