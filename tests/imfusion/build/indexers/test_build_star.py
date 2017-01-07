import argparse

import pyfaidx
import pytest

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path

from imfusion.build.indexers import star

# pylint: disable=no-self-use,redefined-outer-name


class TestStarIndexer(object):
    """Tests for StarIndexer class."""

    def test_dependencies(self):
        """Check dependencies."""

        indexer = star.StarIndexer()
        assert indexer.dependencies == ['STAR']

    def test_build(self, build_kws, mocker):
        """Tests build using example files."""

        # Mock STAR call.
        mock = mocker.patch.object(star.shell, 'run_command')

        # Build reference.
        indexer = star.StarIndexer()
        indexer.build(**build_kws)

        # Check if reference files exist.
        ref = star.StarReference(build_kws['output_dir'])

        assert ref.base_path.exists()
        assert ref.fasta_path.exists()
        assert ref.gtf_path.exists()
        assert ref.indexed_gtf_path.exists()
        assert ref.index_path.exists()
        assert ref.transposon_name == 'T2onc'
        assert ref.transposon_path.exists()
        assert ref.features_path.exists()

        # Check presence of augmented reference sequences.
        refseq = pyfaidx.Fasta(str(ref.fasta_path))
        assert sorted(refseq.keys()) == ['1', '2', 'T2onc']

        # Check call to STAR for building the index.
        mock.assert_called_once_with(
            [
                'STAR', '--runMode', 'genomeGenerate', '--genomeDir',
                str(ref.index_path), '--genomeFastaFiles', str(ref.fasta_path),
                '--sjdbGTFfile', str(ref.gtf_path), '--sjdbOverhang', '100'
            ],
            stdout=pytest.helpers.mock_any(object))

    def test_from_args(self, cmdline_args):
        """Tests creation from command line."""

        # Setup parser.
        parser = argparse.ArgumentParser()
        star.StarIndexer.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        indexer = star.StarIndexer.from_args(args)

        # Check args.
        assert args.reference_seq == Path('/path/to/ref')
        assert args.reference_gtf == Path('/path/to/gtf')
        assert args.transposon_seq == Path('/path/to/tr')
        assert args.transposon_features == Path('/path/to/feat')

        # Check instance.
        assert indexer.overhang == 100

    def test_from_args_extra(self, cmdline_args):
        """Tests creation from command line using extra args."""

        cmdline_args += ['--overhang', '200']

        # Setup parser.
        parser = argparse.ArgumentParser()
        star.StarIndexer.configure_args(parser)

        # Construct aligner.
        args = parser.parse_args(cmdline_args)
        indexer = star.StarIndexer.from_args(args)

        # Check instance.
        assert indexer.overhang == 200


class TestStarReference(object):
    """Tests for StarReference class."""

    def test_non_existing(self, tmpdir):
        """Tests non-existing reference"""

        with pytest.raises(ValueError):
            star.StarReference(tmpdir / 'test')
