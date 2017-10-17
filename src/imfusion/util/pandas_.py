from future.utils import native_str
import pandas as pd


class RecordSet(object):
    """Base class that provides functionality for serializing and
       deserializing namedtuple records into a DataFrame format.

    Subclasses should override the ``_tuple_class`` method to
    return the namedtuple class that should be used as a record.
    """

    def __init__(self, values):
        self._values = self._check_frame(values)

    @classmethod
    def _check_frame(cls, values):
        if values.shape[0] > 0:
            fields = cls._tuple_fields()

            for field in fields:
                if field not in values.columns:
                    raise ValueError(
                        'Missing required column {}'.format(field))

        return values.reindex(columns=fields)

    @classmethod
    def _tuple_class(cls):
        """Returns namedtuple class used to instantiate records."""
        raise NotImplementedError()

    @classmethod
    def _tuple_fields(cls):
        """Returns the fields in the named tuple class."""
        return cls._tuple_class()._fields

    @property
    def values(self):
        """Internal DataFrame representation of records."""
        return self._values

    def __getitem__(self, item):
        return self._values[item]

    def __setitem__(self, idx, value):
        self._values[idx] = value

    def __len__(self):
        return len(self._values)

    @property
    def loc(self):
        """Label-based indexer (similar to pandas .loc)."""
        return LocWrapper(self._values.loc, constructor=self._loc_constructor)

    @property
    def iloc(self):
        """Label-based indexer (similar to pandas .loc)."""
        return LocWrapper(self._values.iloc, constructor=self._loc_constructor)

    def _loc_constructor(self, values):
        if len(values.shape) != 2:
            return values
        return self.__class__(values)

    @classmethod
    def from_tuples(cls, tuples):
        """Builds a record set instance from the given tuples."""
        records = (tup._asdict() for tup in tuples)
        return cls(pd.DataFrame.from_records(records))

    def to_tuples(self):
        """Converts the record set into an iterable of tuples."""

        tuple_class = self._tuple_class()

        for row in self._values.itertuples():
            row_dict = row._asdict()
            row_dict.pop('Index', None)

            yield tuple_class(**row_dict)

    @classmethod
    def from_csv(cls, file_path, **kwargs):
        """Reads a record set from a csv file using pandas.read_csv."""
        values = pd.read_csv(native_str(file_path), **kwargs)
        return cls(values)

    def to_csv(self, file_path, **kwargs):
        """Writes the record set to a csv file using pandas' to_csv."""
        self._values.to_csv(native_str(file_path), **kwargs)

    def groupby(self, by, **kwargs):
        """Groups the set by values of the specified columns."""
        for key, group in self._values.groupby(by, **kwargs):
            yield key, self.__class__(group)

    def query(self, expr, **kwargs):
        """Queries the columns of the set with a boolean expression."""
        return self.__class__(self._values.query(expr, **kwargs))

    def sort_values(self, by, **kwargs):
        return self.__class__(self._values.sort_values(by=by, **kwargs))

    @classmethod
    def concat(cls, record_sets):
        """Concatenates multiple records sets into a single set."""
        return cls(pd.concat((rs.values for rs in record_sets), axis=0))


class LocWrapper(object):
    """Wrapper class that wraps an objects loc/iloc accessor."""

    def __init__(self, loc, constructor=None):
        if constructor is None:
            constructor = lambda x: x

        self._loc = loc
        self._constructor = constructor

    def __getitem__(self, item):
        result = self._loc[item]
        return self._constructor(result)


class MetadataRecordSet(RecordSet):
    """Base RecordSet that supports record metadata.

    Extension of the RecordSet class, which assumes that records contain
    a dict 'metadata' field which contains variable metadata. The
    MetadataRecordSet class ensures that this data is expanded from the
    original record when converted to the set's DataFrame format, and
    converted back again when transforming back to tuples.
    """

    METADATA_FIELD = 'metadata'

    @property
    def metadata_columns(self):
        """Available metadata columns."""
        return set(self._values.columns) - set(self._tuple_fields())

    @classmethod
    def _check_frame(cls, values):
        fields = [
            field for field in cls._tuple_fields()
            if field != cls.METADATA_FIELD
        ]

        if values.shape[0] > 0:
            for field in fields:
                if field not in values.columns:
                    raise ValueError(
                        'Missing required column {}'.format(field))

        extra_cols = set(values.columns) - set(fields)
        col_order = list(fields) + sorted(extra_cols)

        return values.reindex(columns=col_order)

    @classmethod
    def from_tuples(cls, tuples):
        """Builds a record set instance from the given tuples."""

        metadata_field = cls.METADATA_FIELD

        def _to_record(tup):
            record = tup._asdict()
            record.update(record.pop(metadata_field))
            return record

        records = (_to_record(tup) for tup in tuples)
        return cls(pd.DataFrame.from_records(records))

    def to_tuples(self):
        """Converts the record set into an iterable of tuples."""

        tuple_class = self._tuple_class()

        for row in self._values.itertuples():
            row_dict = row._asdict()
            row_dict.pop('Index', None)

            yield tuple_class(**row_dict)
