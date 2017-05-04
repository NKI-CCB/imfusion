import pytest

from imfusion.util import tabix


@pytest.fixture
def gtf_path():
    """Returns path to example gtf file containing Smn1."""
    return pytest.helpers.data_path(
        'Mus_musculus.GRCm38.76.Smn1.gtf', relative_to=__file__)


class TestReadGtfFrame(object):
    """Unit tests for the read_gtf_frame function."""

    def test_example(self, gtf_path):
        """Test example gtf containing Smn1."""

        gtf = tabix.read_gtf_frame(gtf_path)

        # Check shape.
        assert gtf.shape == (76, 9)
        assert list(gtf.columns) == tabix.GTF_COLUMNS

        # Check if start position = -1.
        assert gtf.iloc[0]['start'] == 100124851


class TestFlattenGtfFrame(object):
    """Unit tests for the flatten_gtf_frame function."""

    def test_example(self, gtf_path):
        """Test example gtf containing Smn1."""

        gtf = tabix.read_gtf_frame(gtf_path)
        gtf_flat = tabix.flatten_gtf_frame(gtf)

        # Check shape.
        assert gtf_flat.shape == (25, 9)
        assert list(gtf_flat.columns) == tabix.GTF_COLUMNS

        # Check positions.
        head = gtf_flat.head()
        assert list(head['start']) == [
            100124851, 100124857, 100124893, 100124923, 100125818]
        assert list(head['end']) == [
            100124857, 100124893, 100124923, 100124965, 100125931]

        # Check other values.
        first = gtf_flat.iloc[0]
        assert first.seqname == '13'
        assert first.source == 'imfusion_flatten'
        assert first.feature == 'exonic_part'
        assert first.score == '.'
        assert first.frame == '.'
        assert first.attribute == ('gene_id "ENSMUSG00000021645"; '
                                   'exonic_part_number "001"')

# class TestFlattenGtfFrame(object):
#     """Unit tests for flatten_intervals function."""

#     def test_smn1(self, test_intervals, flattened_test_intervals):
#         """Test Smn1 example."""
#         flattened = list(tabix.flatten_gtf_frame(test_intervals))
#         assert flattened == flattened_test_intervals

    # def test_empty(self):
    #     """Test empty example."""
    #     flattened = list(tabix._flatten_intervals([]))
    #     assert flattened == []

