def reorder_columns(df, order):
    """Reorders dataframe columns, sorting any extra columns alphabetically."""
    extra_cols = set(df.columns) - set(order)
    return df[list(order) + sorted(extra_cols)]
