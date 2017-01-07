import argparse
import shutil

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from frozendict import frozendict
import pytest

from imfusion.insertions.aligners import tophat
from imfusion.model import Insertion

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def read_paths():
    """Example read paths."""
    return [(Path('a.fastq.gz'), Path('b.fastq.gz'))]


@pytest.fixture
def tophat_path():
    """Return path to fusions.out file."""
    return pytest.helpers.data_path('fusions.out', relative_to=__file__)


@pytest.fixture
def tophat_reference():
    """Returns example reference for Tophat2."""
    ref_path = pytest.helpers.data_path('reference', relative_to=__file__)
    return tophat.TophatReference(ref_path)


@pytest.fixture
def tophat_output_dir(tmpdir, tophat_path):
    """Simulates Tophat2 output directory."""

    # Create directories.
    output_dir = Path(str(tmpdir / 'out'))

    tophat_dir = output_dir / '_tophat'
    tophat_dir.mkdir(parents=True)

    # Copy / simulate aligner output files.
    shutil.copy(str(tophat_path),
                str(tophat_dir / 'fusions.out'))  # yapf: disable

    pytest.helpers.touch(tophat_dir / 'Aligned.sortedByCoord.out.bam')

    return output_dir


@pytest.fixture
def cmdline_args(tmpdir):
    """Example command line arguments."""
    return [
        '--fastq', 'a.fastq.gz', '--fastq2', 'b.fastq.gz', '--reference',
        str(tmpdir), '--output_dir', '/path/to/out'
    ]


class TestTophatAligner(object):
    """Tests for the TophatAligner class."""

    def test_dependencies(self, tophat_reference):
        """Test dependencies for various configurations."""

        # Basic call.
        aligner = tophat.TophatAligner(tophat_reference)
        assert set(aligner.dependencies) == {'tophat2', 'bowtie'}

        # With assembly.
        aligner = tophat.TophatAligner(tophat_reference, assemble=True)
        assert set(aligner.dependencies) == {'tophat2', 'bowtie', 'stringtie'}

    def test_identify_insertions(self, read_paths, tophat_reference,
                                 tophat_output_dir, mocker):
        """Integration test using data from an example run."""

        # TODO: Check assembly?

        # Mock star_align call.
        tophat_mock = mocker.patch.object(tophat, 'tophat2_align')

        # Call identify insertions.
        aligner = tophat.TophatAligner(tophat_reference)
        ins = list(aligner.identify_insertions(read_paths, tophat_output_dir))

        # Check call to star_align.
        tophat_mock.assert_called_once_with(
            read_paths,
            output_dir=tophat_output_dir / '_tophat',
            index_path=tophat_reference.index_path,
            extra_args={
                '--fusion-search': (),
                '--transcriptome-index':
                (str(tophat_reference.transcriptome_path), ),
                '--num-threads': (1, ),
                '--bowtie1': (),
                '--fusion-anchor-length': (12, )
            },
            logger=pytest.helpers.mock_any(object))

        # Check result, including specific Cblb insertion.
        assert len(ins) == 7

        assert ins[2] == Insertion(
            id='INS_4',
            seqname='16',
            position=52141093,
            strand=-1,
            junction_support=462,
            spanning_support=103,
            support=565,
            metadata=frozendict({
                'transposon_anchor': 1539,
                'gene_name': 'Cblb',
                'feature_name': 'En2SA',
                'gene_strand': 1,
                'feature_type': 'SA',
                'orientation': 'antisense',
                'feature_strand': -1
            }))

    def test_from_args_basic(self, cmdline_args):
        """Tests creation with minimal arguments."""

        # Setup parser.
        parser = argparse.ArgumentParser()
        tophat.TophatAligner.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        aligner = tophat.TophatAligner.from_args(args)

        # Check args.
        assert args.fastq == [Path('a.fastq.gz')]
        assert args.fastq2 == [Path('b.fastq.gz')]
        assert args.output_dir == Path('/path/to/out')

        # Check aligner.
        # pylint: disable=w0212
        assert not aligner._assemble
        assert aligner._assemble_args == {}
        assert aligner._min_flank == 12
        assert aligner._threads == 1
        assert aligner._extra_args == {}
        assert aligner._filter_features
        assert aligner._filter_orientation
        assert aligner._filter_blacklist is None

    def test_from_args_extra_args(self, cmdline_args):
        """Tests creation with extra options."""

        cmdline_args += [
            '--tophat_threads', '5', '--tophat_min_flank', '20',
            '--tophat_args', "--limitBAMsortRAM 2000", '--assemble',
            '--no_filter_orientation', '--no_filter_feature',
            '--blacklisted_genes', 'En2'
        ] # yapf:disable

        # Setup parser.
        parser = argparse.ArgumentParser()
        tophat.TophatAligner.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        aligner = tophat.TophatAligner.from_args(args)

        # Check args.
        assert args.fastq == [Path('a.fastq.gz')]
        assert args.fastq2 == [Path('b.fastq.gz')]
        assert args.output_dir == Path('/path/to/out')

        # Check aligner.
        # pylint: disable=w0212
        assert aligner._assemble
        assert aligner._assemble_args == {}
        assert aligner._min_flank == 20
        assert aligner._threads == 5
        assert aligner._extra_args == {'--limitBAMsortRAM': ('2000', )}
        assert not aligner._filter_features
        assert not aligner._filter_orientation
        assert aligner._filter_blacklist == ['En2']


@pytest.fixture
def tophat_align_kws(tmpdir):
    """Example tophat2_align keyword arguments."""
    return {
        'fastqs': [('in.R1.fastq.gz', 'in.R2.fastq.gz')],
        'index_path': '/path/to/index',
        'output_dir': Path(str(tmpdir / '_tophat')),
        'extra_args': None
    }


class TestTophat2Align(object):
    """Tests for tophat2_align function."""

    def test_basic(self, tophat_align_kws, mocker):
        """Test basic invocation."""

        mock_align = mocker.patch.object(tophat.shell, 'run_command')
        tophat.tophat2_align(**tophat_align_kws)

        # Check if output directory was created.
        assert tophat_align_kws['output_dir'].exists()

        # Check arguments to run command. Checks if we have an extra slash
        # added to output dir, if the fastqs were added properly, and if
        # readFilesCommand is added for the gzipped files.
        mock_align.assert_called_once_with(
            args=[
                'tophat2', '--output-dir', str(tophat_align_kws['output_dir']),
                '/path/to/index', str(tophat_align_kws['fastqs'][0][0]),
                str(tophat_align_kws['fastqs'][0][1])
            ],
            stdout=None,
            stderr=None,
            logger=pytest.helpers.mock_any(object))

    def test_extra_args(self, tophat_align_kws, mocker):
        """Test adding extra arguments."""

        tophat_align_kws['extra_args'] = {
            '--mate-std-dev': (20, ),
            '--num-threads': (10, )
        }

        mock_align = mocker.patch.object(tophat.shell, 'run_command')
        tophat.tophat2_align(**tophat_align_kws)

        args = pytest.helpers.mock_call(mock_align)[1]['args']

        # Check new argument.
        assert '--mate-std-dev' in args
        assert args[args.index('--mate-std-dev') + 1] == '20'

        assert '--num-threads' in args
        assert args[args.index('--num-threads') + 1] == '10'

    def test_unpaired(self, tophat_align_kws, mocker):
        """Test single end case."""

        tophat_align_kws['fastqs'] = ['in.fastq.gz']

        mock_align = mocker.patch.object(tophat.shell, 'run_command')
        tophat.tophat2_align(**tophat_align_kws)

        args = pytest.helpers.mock_call(mock_align)[1]['args']
        assert args[-1] == 'in.fastq.gz'

    def test_multiple(self, tophat_align_kws, mocker):
        """Test case with multiple fastq files."""

        tophat_align_kws['fastqs'] = [('in1.R1.fastq.gz', 'in1.R2.fastq.gz'),
                                      ('in2.R1.fastq.gz', 'in2.R2.fastq.gz')]

        mock_align = mocker.patch.object(tophat.shell, 'run_command')
        tophat.tophat2_align(**tophat_align_kws)

        args = pytest.helpers.mock_call(mock_align)[1]['args']
        assert args[-2] == 'in1.R1.fastq.gz,in2.R1.fastq.gz'
        assert args[-1] == 'in1.R2.fastq.gz,in2.R2.fastq.gz'
