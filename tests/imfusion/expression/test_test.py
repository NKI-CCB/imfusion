import pytest

from frozendict import frozendict

from imfusion.expression.counts import read_exon_counts
from imfusion.expression import test
from imfusion.model import Insertion

# pylint: disable=no-self-use,redefined-outer-name


@pytest.fixture
def insertions():
    """Example insertions in ENSMUSG00000051951 (5 exons)."""
    return [
        # After third exon.
        Insertion(
            id='1',
            seqname='1',
            position=3207327,
            strand=1,
            junction_support=1,
            spanning_support=1,
            support=2,
            metadata=frozendict({
                'gene_id': 'ENSMUSG00000051951',
                'sample': 'S1'})),

        # After first exon.
        Insertion(
            id='2',
            seqname='1',
            position=3214491,
            strand=-1,
            junction_support=1,
            spanning_support=1,
            support=2,
            metadata=frozendict({
                'gene_id': 'ENSMUSG00000051951',
                'sample': 'S2'}))
    ] # yapf: disable


@pytest.fixture
def exon_counts():
    """Example exon counts."""

    file_path = pytest.helpers.data_path(
        'exon_counts.txt', relative_to=__file__)
    return read_exon_counts(file_path)


class TestSplitCounts(object):
    """Tests for split_counts function."""

    def test_proper_example(self, insertions, exon_counts):
        """Tests example with two insertions, after exons 2 and 4."""

        counts = exon_counts.ix['ENSMUSG00000051951']
        before, after, dropped = test.split_counts(counts, insertions)

        assert len(before) == 1
        assert len(after) == 2
        assert dropped == set()

    def test_before_gene(self, insertions, exon_counts):
        """Tests dropping sample with an insertion before the gene"""

        insertions[1] = insertions[1]._replace(position=3215652)

        counts = exon_counts.ix['ENSMUSG00000051951']
        before, after, dropped = test.split_counts(counts, insertions)

        assert len(before) == 3
        assert len(after) == 2
        assert dropped == {'S2'}

    def test_after_gene(self, insertions, exon_counts):
        """Tests dropping sample with an insertion after the gene"""

        insertions[0] = insertions[0]._replace(position=3205801)

        counts = exon_counts.ix['ENSMUSG00000051951']
        before, after, dropped = test.split_counts(counts, insertions)

        assert len(before) == 1
        assert len(after) == 4
        assert dropped == {'S1'}

    def test_in_exon(self, insertions, exon_counts):
        """Tests insertion in fourth exon."""

        insertions[0] = insertions[0]._replace(position=3207217)

        counts = exon_counts.ix['ENSMUSG00000051951']
        before, after, dropped = test.split_counts(counts, insertions)

        assert len(before) == 1
        assert len(after) == 1
        assert dropped == set()


@pytest.fixture
def test_insertions():
    """Example insertions for Trp53bp2 and Nf1."""

    file_path = pytest.helpers.data_path(
        'insertions.txt', relative_to=__file__)
    return list(Insertion.from_csv(file_path, sep='\t'))


@pytest.fixture
def test_exon_counts():
    """Example exon counts for Trp53bp2 and Nf1."""

    file_path = pytest.helpers.data_path(
        'exon_counts_test.txt', relative_to=__file__)
    return read_exon_counts(file_path)


class TestDeExon(object):
    """Tests for de_exon function."""

    def test_pos_example_trp53bp2(self, test_insertions, test_exon_counts):
        """Tests positive example of DE in Trp53bp2 in the SB dataset."""

        result = test.de_exon(
            test_insertions, test_exon_counts, gene_id='ENSMUSG00000026510')

        assert result.p_value < 0.01
        assert result.direction == 1

    def test_neg_example_nf1(self, test_insertions, test_exon_counts):
        """Tests negative example of DE in Nf1 in the SB dataset."""

        result = test.de_exon(
            test_insertions, test_exon_counts, gene_id='ENSMUSG00000020716')

        assert result.p_value > 0.05
        assert result.direction == 1

    def test_plot_boxplot(self, test_insertions, test_exon_counts):
        """Exercises plotting of boxplot."""

        result = test.de_exon(
            test_insertions, test_exon_counts, gene_id='ENSMUSG00000026510')
        result.plot_boxplot()

    def test_plot_sums(self, test_insertions, test_exon_counts):
        """Exercises plotting of sums."""

        result = test.de_exon(
            test_insertions, test_exon_counts, gene_id='ENSMUSG00000026510')
        result.plot_sums()
