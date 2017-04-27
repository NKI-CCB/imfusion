from imfusion.util import genomic


class TestFlattenIntervals(object):
    """Unit tests for flatten_intervals function."""

    def test_example(self):
        """Tests example case."""

        ranges = [(10, 12), (12, 15), (15, 20), (20, 30), (30, 32), (45, 55)]
        flattened = list(genomic.flatten_intervals(ranges))

        assert flattened == [(10, 12), (12, 15), (15, 20), (20, 30), (30, 32),
                             (45, 55)]

    def test_empty(self):
        flattened = list(genomic.flatten_intervals([]))
        assert flattened == []
