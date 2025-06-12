import pandas as pd
import pyranges as pr

def pyranges_column_name_change(df: pd.DataFrame) -> pd.DataFrame:
    """
    Changes the column names of a pandas DataFrame to match PyRanges expectations.
    """
    df_renamed = df.rename(columns={
        "sseqid": "Chromosome",
        "sstart": "Start",
        "send": "End"
    })

    return df_renamed


def get_interval_coincidence(df: pd.DataFrame, interval_df: pd.DataFrame) -> pd.DataFrame:
    """
    Uses PyRanges to find overlapping intervals between df and interval_df based on sseqid, start, end.
    Returns the rows in df that overlap with any row in interval_df.
    """

    # Rename columns to match PyRanges expectations
    df_renamed = pyranges_column_name_change(df)

    interval_df_renamed = pyranges_column_name_change(interval_df)

    # Convert to PyRanges
    pr_df = pr.PyRanges(df_renamed)
    pr_intervals = pr.PyRanges(interval_df_renamed)

    # Perform overlap (intersection)
    overlaps = pr_df.overlap(pr_intervals, strandedness=False)

    # Return as DataFrame and rename back
    result = overlaps.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result


def get_interval_not_coincidence(df: pd.DataFrame, interval_df: pd.DataFrame) -> pd.DataFrame:
    """
    Uses PyRanges to find overlapping intervals between df and interval_df based on sseqid, start, end.
    Returns the rows in df that do not overlap with any row in interval_df.
    """
    # Rename columns to match PyRanges expectations
    df_renamed = pyranges_column_name_change(df)

    interval_df_renamed = pyranges_column_name_change(interval_df)

    # Convert to PyRanges
    pr_df = pr.PyRanges(df_renamed)
    pr_intervals = pr.PyRanges(interval_df_renamed)

    # Perform overlap (intersection)
    non_overlapping = pr_df.overlap(pr_intervals, strandedness=False, invert=True)

    # Return as DataFrame and rename back
    result = non_overlapping.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result
