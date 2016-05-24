# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

import mock
import pytest

import pandas as pd

from imfusion.aligners.tophat2 import insertions as _mod
from imfusion.model import Fusion, Insertion

# pylint: disable=redefined-outer-name,no-self-use,too-few-public-methods


@pytest.fixture
def gtf_path():
    """Fixture for gtf file path."""
    return pytest.helpers.data_path('mm10.test.gtf.gz')


class TestIdentifyInsertions(object):
    """Tests for identify_insertions."""

    @pytest.fixture
    def fusion(self):
        return Fusion(seqname='1', anchor_genome=182436900,
                      anchor_transposon=420, strand_genome=1,
                      strand_transposon=1, flank_genome=35,
                      flank_transposon=35, spanning_reads=15,
                      supporting_mates=0, spanning_mates=0)


    @pytest.fixture
    def fusion_no_feat(self):
        return Fusion(seqname='1', anchor_genome=182436900,
                      anchor_transposon=10, strand_genome=1,
                      strand_transposon=1, flank_genome=35,
                      flank_transposon=35, spanning_reads=15,
                      supporting_mates=0, spanning_mates=0)


    @pytest.fixture
    def fusion_wrong_feat(self):
        return Fusion(seqname='1', anchor_genome=182436900,
                      anchor_transposon=110, strand_genome=1,
                      strand_transposon=1, flank_genome=35,
                      flank_transposon=35, spanning_reads=15,
                      supporting_mates=0, spanning_mates=0)

    @pytest.fixture
    def fusion_no_gene(self):
        return Fusion(seqname='1', anchor_genome=200, anchor_transposon=420,
                      strand_genome=1, strand_transposon=1, flank_genome=35,
                      flank_transposon=35, spanning_reads=15,
                      supporting_mates=0, spanning_mates=0)

    @pytest.fixture
    def transposon_features(self):
        """Fixture for example transposon features."""
        return pd.DataFrame({
            'name': ['En2SA', 'SA', 'SD', 'MSCV'],
            'start': [1479, 232, 1007, 100],
            'end': [1824, 428, 1184, 120],
            'strand': [-1, 1, 1, 1],
            'type': ['SA', 'SA', 'SD', 'Promotor']
        })


    @pytest.fixture
    def identify_kws(self, gtf_path, transposon_features, tmpdir):
        return {
            'fastqs': [Path('a.fastq'), Path('b.fastq')],
            'index_path': Path('/path/to/index'),
            'reference_gtf_path': gtf_path,
            'transposon_name': 'T2onc',
            'transposon_features': transposon_features,
            'work_dir': Path(str(tmpdir)),
            'min_flank': 20,
            'sample_id': 's1'
    }

    def test_example(self, identify_kws, fusion):
        """Integration test with an example fusion."""

        # Patch identify_fusions so that we don't actually
        # have to call TopHat2 on the command line.
        with mock.patch.object(_mod, 'identify_fusions',
                               return_value=[fusion]) as func:

            # Call insertions.
            insertions = list(_mod.identify_insertions(**identify_kws))

            # Check call.
            func.assert_called_once_with(
                identify_kws['fastqs'],
                index_path=identify_kws['index_path'],
                transposon_name=identify_kws['transposon_name'],
                work_dir=identify_kws['work_dir'],
                reference_gtf_path=identify_kws['reference_gtf_path'],
                min_flank=identify_kws['min_flank'],
                tophat_kws=None,
                transcriptome_index=None)

            # Check for our insertion.
            assert len(insertions) == 1

            # Check insertion attributes.
            ins = insertions[0]
            assert ins.seqname == fusion.seqname
            assert ins.position > fusion.anchor_genome
            assert ins.strand == 1
            assert ins.sample_id == identify_kws['sample_id']
            assert ins.gene_id == 'ENSMUSG00000026510'
            assert ins.gene_name == 'Trp53bp2'
            assert ins.orientation == 'sense'
            assert ins.feature_name == 'SA'
            assert ins.anchor_genome == fusion.anchor_genome
            assert ins.anchor_transposon == fusion.anchor_transposon
            assert ins.flank_genome == fusion.flank_genome
            assert ins.flank_transposon == fusion.flank_transposon
            assert ins.spanning_reads == fusion.spanning_reads
            assert ins.supporting_mates == fusion.supporting_mates
            assert ins.spanning_mates == fusion.spanning_mates

    def test_no_feature(self, identify_kws, fusion_no_feat):
        """Integration test with an example without feature."""

        # Patch identify_fusions so that we don't actually
        # have to call TopHat2 on the command line.
        with mock.patch.object(_mod, 'identify_fusions',
                               return_value=[fusion_no_feat]):
            insertions = list(_mod.identify_insertions(**identify_kws))
            assert len(insertions) == 0

    def test_wrong_feature(self, identify_kws, fusion_wrong_feat):
        """Integration test with an example with a non SA/SD feature."""

        # Patch identify_fusions so that we don't actually
        # have to call TopHat2 on the command line.
        with mock.patch.object(_mod, 'identify_fusions',
                               return_value=[fusion_wrong_feat]):
            insertions = list(_mod.identify_insertions(**identify_kws))
            assert len(insertions) == 0

    def test_no_gene(self, identify_kws, fusion_no_gene):
        """Integration test with an example without an annotated gene."""

        # Patch identify_fusions so that we don't actually
        # have to call TopHat2 on the command line.
        with mock.patch.object(_mod, 'identify_fusions',
                               return_value=[fusion_no_gene]):
            insertions = list(_mod.identify_insertions(**identify_kws))
            assert len(insertions) == 1
            assert insertions[0].gene_name is None
            assert insertions[0].orientation is None


@pytest.fixture
def fusion_path():
    """Fixture for example fusions."""
    return pytest.helpers.data_path('tophat/fusions.out',
                                    relative_to=__file__)

class TestReadFusions(object):
    """Test read_fusions."""

    def test_example(self, fusion_path):
        """Tests an example file."""

        fusions = _mod.read_fusions(fusion_path)

        # Check if we read four fusions and have 12 columns.
        assert fusions.shape == (4, 12)

        # Check one fusion for consistency.
        fusion = fusions.ix[2]
        assert fusion.seqname_a == '1'
        assert fusion.location_a == 130717807
        assert fusion.strand_a == 1
        assert fusion.seqname_b == 'T2onc'
        assert fusion.location_b == 1542
        assert fusion.strand_b == -1
        assert fusion.supp_reads == 15
        assert fusion.supp_mates == 0
        assert fusion.supp_spanning_mates == 0
        assert fusion.contradicting_reads == 20
        assert fusion.flank_a == 35
        assert fusion.flank_b == 35

    def test_read_fusions_empty(self, tmpdir):
        """Test with an empty file."""

        # Create empty dummy file.
        fusion_path = (tmpdir / 'fusions.out')
        with fusion_path.open('w'):
            pass

        # Check if this is correctly read.
        fusions = _mod.read_fusions(fusion_path)
        assert len(fusions) == 0


class TestExtractFusions(object):

    def test_example(self, fusion_path):
        """Tests extract_fusions function with an example file."""

        fusions = list(_mod.extract_fusions(fusion_path, 'T2onc'))

        # Check if we find two gene-transposon fusions.
        assert len(fusions) == 2

        # Check one fusion for consistency.
        fusion = fusions[0]
        assert fusion.seqname == '1'
        assert fusion.anchor_genome == 130717807
        assert fusion.anchor_transposon == 1542
        assert fusion.strand_genome == 1
        assert fusion.strand_transposon == -1
        assert fusion.spanning_reads == 15
        assert fusion.supporting_mates == 0
        assert fusion.spanning_mates == 0
        assert fusion.flank_genome == 35
        assert fusion.flank_transposon == 35

    def test_without_valid(self, fusion_path):
        """Tests extract_fusions without any valid fusions."""

        fusions = list(_mod.extract_fusions(fusion_path, 'INVALID'))
        assert len(fusions) == 0


class TestIdentifyFusions(object):

    TOPHAT_FILES = pytest.mark.datafiles(
        str(pytest.helpers.data_path('tophat', relative_to=__file__)),
        keep_top_dir=True)

    @TOPHAT_FILES
    @pytest.fixture()
    def identify_kwargs(self, gtf_path, datafiles):
        """Default kwargs for identify_fusions."""

        # Create mock index.
        index_path = Path(str(datafiles / 'test'))
        index_path.with_suffix('.1.ebwt').touch()

        return {
            'fastqs': [Path('a.fastq'), Path('b.fastq')],
            'index_path': index_path,
            'reference_gtf_path': gtf_path,
            'transposon_name': 'T2onc',
            'work_dir': Path(str(datafiles)),
            'use_existing': False
        }

    @pytest.fixture
    def tr_index_path(self, tmpdir):
        base_path = Path(str(tmpdir / 'tr_index'))
        Path(str(base_path) + '.1.ebwt').touch()
        return base_path

    @TOPHAT_FILES
    def test_example(self, identify_kwargs):
        """Checks basic call."""

        with mock.patch.object(_mod, 'tophat2_align') as mock_:
            # Call identify fusions with use_existing is False, to ensure
            # that we simulate the call to tophat2_align.
            fusions = list(_mod.identify_fusions(**identify_kwargs))

            # Check arguments to tophat_align.
            _, kwargs = pytest.helpers.mock_call(mock_)
            assert kwargs['fastqs'] == identify_kwargs['fastqs']
            assert kwargs['index_path'] == identify_kwargs['index_path']
            assert kwargs['output_dir'] == \
                identify_kwargs['work_dir'] / 'tophat'

            tophat_kws = kwargs['kwargs']
            assert '--fusion-search' in tophat_kws
            assert '--bowtie1' in tophat_kws
            assert tophat_kws['--fusion-anchor-length'] == 20

            # Check returned fusions.
            assert len(fusions) == 2

            # Check if alignment was symlinked.
            assert kwargs['output_dir'] / 'alignment.bam'

    @TOPHAT_FILES
    def test_identify_fusions_tr_index(self, identify_kwargs, tr_index_path):
        """Checks call with transcriptome-index."""

        with mock.patch.object(_mod, 'tophat2_align') as mock_:
            identify_kwargs['transcriptome_index'] = tr_index_path
            list(_mod.identify_fusions(**identify_kwargs))

            # Check arguments to tophat2.
            tophat_kws = pytest.helpers.mock_call(mock_)[1]['kwargs']
            assert '--transcriptome-index' in tophat_kws
            assert tophat_kws['--transcriptome-index'] == str(tr_index_path)

    def test_identify_fusions_existing(self, identify_kwargs, fusion_path):
        """Checks if tophat is called if use_existing is True but output
        does not exist and therefore should still be generated.
        """

        identify_kwargs['use_existing'] = True
        fusions_ = _mod.read_fusions(fusion_path)

        with mock.patch.object(_mod, 'tophat2_align') as mock_, \
             mock.patch.object(_mod, 'read_fusions', return_value=fusions_):
            list(_mod.identify_fusions(**identify_kwargs))
            assert mock_.call_count > 0


class TestTophat2Align(object):
    """Tests for tophat2_align."""

    @pytest.fixture
    def tophat2_kwargs(self, tmpdir):
        """Default tophat2 kwargs."""
        return {
            'fastqs': [Path('a.fastq'), Path('b.fastq')],
            'index_path': Path('/path/to/index'),
            'output_dir': tmpdir,
            'path': None,
            'check_python': False
        }

    def test_single(self, tophat2_kwargs):
        """Tests single-end example arguments."""

        with mock.patch('subprocess.check_call') as check_call:
            # Call tophat2 align.
            out_dir = _mod.tophat2_align(**tophat2_kwargs)
            assert out_dir == tophat2_kwargs['output_dir']

            # Extract call arguments.
            kwargs = check_call.call_args_list[0][1]
            args, stderr = kwargs['args'], kwargs['stderr']

            # Check executed call.
            assert args == ['tophat2', '--output-dir',
                            str(tophat2_kwargs['output_dir']),
                            str(tophat2_kwargs['index_path']),
                            ','.join(map(str, tophat2_kwargs['fastqs']))]
            assert stderr is not None

    def test_paired(self, tophat2_kwargs):
        """Tests if paired-end fastqs are handled correctly."""

        tophat2_kwargs['fastqs'] = [(Path('a.R1.fastq'), Path('a.R2.fastq')),
                                    (Path('b.R1.fastq'), Path('b.R2.fastq'))]

        with mock.patch('subprocess.check_call') as check_call:
            # Call tophat2 align.
            _mod.tophat2_align(**tophat2_kwargs)

            # Extract call arguments.
            args = check_call.call_args_list[0][1]['args']

            # Check if pairs were handled correctly.
            assert args[-2] == 'a.R1.fastq,b.R1.fastq'
            assert args[-1] == 'a.R2.fastq,b.R2.fastq'

    def test_extra_options(self, tophat2_kwargs):
        tophat2_kwargs['kwargs'] = {'--num-threads': 1,
                                    '--fusion-search': True,
                                    '--flag': False}

        with mock.patch('subprocess.check_call') as check_call:
            # Call tophat2 align.
            _mod.tophat2_align(**tophat2_kwargs)

            # Extract call arguments.
            args = check_call.call_args_list[0][1]['args']

            # Check if thread argument was inserted.
            thread_index = args.index('--num-threads')
            assert args[thread_index + 1] == '1'

            # Check if fusion flag was inserted without value.
            fusion_index = args.index('--fusion-search')
            assert args[fusion_index + 1] != 'True'

            # Check if negative flag was omitted.
            with pytest.raises(ValueError):
                args.index('--flag')

    def test_non_existing_output_dir(self, tophat2_kwargs, tmpdir):
        output_dir = Path(str(tmpdir / 'subdir'))
        tophat2_kwargs['output_dir'] = output_dir

        with mock.patch('subprocess.check_call'):
            _mod.tophat2_align(**tophat2_kwargs)
            assert output_dir.exists()
