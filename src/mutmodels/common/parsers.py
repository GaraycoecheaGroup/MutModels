"""
functions to parse input data
"""

# parser for SigProfiler-Style count matrix tsv files
import pandas as pd
from .mutation_types import SBS96

# class counts for the standard SigProfilerMatrixGenerator matrices
MATRIX_CLASSES = {
    "SBS6": 6, "SBS24": 24, "SBS96": 96, "SBS384": 384, "SBS1536": 1536,
    "DBS78": 78, "DBS186": 186,
    "ID28": 28, "ID83": 83, "ID415": 415,
}

from .mutation_types import SBS96

def align(data, order=SBS96):
    """Reorder a Series (indexed by MutationType) or a DataFrame
    (row-indexed by MutationType) into canonical `order`."""
    order = list(order)                       # tuple -> list: avoids .loc tuple semantics
    missing = set(order) - set(data.index)
    extra   = set(data.index) - set(order)
    if missing or extra:
        raise ValueError(
            f"MutationType mismatch — missing: {sorted(missing)}, extra: {sorted(extra)}"
        )

    if data.index.has_duplicates:
            raise ValueError("Duplicate MutationType labels in index")

    return data.loc[order]                     # preserves Series/DataFrame type


def validate_datatype(data,data_type):
    if data_type is not None:
        if data_type not in MATRIX_CLASSES:
            raise ValueError(f"Unknown data_type {data_type!r}; "
                             f"known types: {', '.join(MATRIX_CLASSES)}")
        expected = MATRIX_CLASSES[data_type]
        if data.shape[0] != expected:
            raise ValueError(f"{data_type} expects {expected} classes, "
                             f"found {data.shape[0]}")

def read_matrix(path, data_type=None,counts=True):
    """Read any SigProfilerMatrixGenerator count matrix (SBS/DBS/ID).

    All these files share one layout: a tab-separated 'MutationType' column of
    mutation classes followed by one integer count column per sample. Only the
    number of classes and the label format differ between matrix types.

    path: path to the matrix file.
    data_type: optional name like 'SBS96' or 'DBS78'. If given, the row count
        is validated against the expected number of classes.

    Returns a DataFrame indexed by mutation class, one column per sample.
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "MutationType"

    validate_datatype(df,data_type)

    # if data_type is not None:
    #     if data_type not in MATRIX_CLASSES:
    #         raise ValueError(f"Unknown data_type {data_type!r}; "
    #                          f"known types: {', '.join(MATRIX_CLASSES)}")
    #     expected = MATRIX_CLASSES[data_type]
    #     if df.shape[0] != expected:
    #         raise ValueError(f"{data_type} expects {expected} classes, "
    #                          f"found {df.shape[0]}")

    # ensure alignment with our mut order
    if data_type == "SBS96":
        df = align(df,order=SBS96)

    if counts:
        df = df.astype(int)

    return df

def read_series(path,data_type=None,counts=False):
    """
    like matrix, but for a single column of values associated with a MutType index
    """
    series = align(pd.read_csv(path,sep='\t',index_col=0).squeeze("columns"))
    series.index.name = "MutationType"

    validate_datatype(series,data_type)

    if data_type == "SBS96":
        series = align(series,order=SBS96)
    
    if counts:
        series = series.astype(int)
    
    return series
