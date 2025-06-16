import pandas as pd
import pyranges as pr


def pyranges_column_name_change(df: pd.DataFrame) -> pd.DataFrame:
    """
    Renames specific columns in a Pandas DataFrame to standardized names related
    to genomic data representation. The function is expected to handle a DataFrame
    containing columns such as 'sseqid', 'sstart', and 'send', and convert these
    to 'Chromosome', 'Start', and 'End', respectively.

    Parameters:
    -----------
    df: pd.DataFrame
        A Pandas DataFrame containing the columns 'sseqid', 'sstart', and 'send'.
        These columns are expected to represent genomic data that will be renamed
        to standardized identifiers for further processing.

    Returns:
    --------
    pd.DataFrame
        A new DataFrame with renamed columns:
        - 'sseqid' renamed to 'Chromosome'
        - 'sstart' renamed to 'Start'
        - 'send' renamed to 'End'

    """
    df_renamed = df.rename(columns={
        "sseqid": "Chromosome",
        "sstart": "Start",
        "send": "End"
    })

    return df_renamed


def get_interval_coincidence(df: pd.DataFrame, interval_df: pd.DataFrame) -> pd.DataFrame:
    """
    Finds overlap intervals between two dataframes using PyRanges.

    This function takes two dataframes as input, converts them into PyRanges
    objects, calculates the overlapping intervals between them, and then returns
    the result as a dataframe. Input dataframes should contain genomic interval
    information. The column names of the dataframes will be preprocessed to match
    the requirements of PyRanges, and the result will have its columns renamed
    back before returning.

    Parameters
    ----------
    df : pd.DataFrame
        The main dataframe containing genomic intervals for which overlap
        will be calculated.
    interval_df : pd.DataFrame
        The dataframe containing intervals to check for overlap with the main
        dataframe.

    Returns
    -------
    pd.DataFrame
        A dataframe containing overlapping genomic intervals with standardized
        column names. The columns of the resulting dataframe will include
        'sseqid', 'sstart', and 'send'.
    """

    # Rename columns to match PyRanges expectations
    df_renamed = pyranges_column_name_change(df)
    interval_df_renamed = pyranges_column_name_change(interval_df)

    # Convert to PyRanges
    pr_df = pr.PyRanges(df_renamed)
    pr_intervals = pr.PyRanges(interval_df_renamed)

    # Perform overlap (intersection)
    overlaps = pr_df.overlap(pr_intervals, strandedness=False) # TODO: strandedness could be a parameter?

    # Return as DataFrame and rename back
    result = overlaps.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result


def get_interval_not_coincidence(df: pd.DataFrame, interval_df: pd.DataFrame) -> pd.DataFrame:
    """
    Gets intervals from the first DataFrame that do not coincide with intervals in
    the second DataFrame. This function converts input data into PyRanges objects to
    efficiently perform interval operations and identifies intervals from the first
    DataFrame that do not overlap with those in the second DataFrame.

    Parameters:
    ----------
        df: DataFrame
            The first DataFrame containing intervals to compare. Expected columns
            are compatible with PyRanges, such as Chromosome, Start, and End.
        interval_df: DataFrame
            The second DataFrame representing intervals for comparison. Similar
            structure to the first DataFrame.

    Returns:
    ----------
        DataFrame
            A DataFrame of intervals from the first DataFrame that do not overlap
            with any intervals in the second DataFrame. The resulting DataFrame has
            columns renamed back to match the original format, such as sseqid,
            sstart, and send.
    """
    # Rename columns to match PyRanges expectations
    df_renamed = pyranges_column_name_change(df)
    interval_df_renamed = pyranges_column_name_change(interval_df)

    # Convert to PyRanges
    pr_df = pr.PyRanges(df_renamed)
    pr_intervals = pr.PyRanges(interval_df_renamed)

    # Perform overlap (intersection)
    non_overlapping = pr_df.overlap(pr_intervals, strandedness=False, invert=True) # TODO: strandedness could be a parameter?

    # Return as DataFrame and rename back
    result = non_overlapping.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result

def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping or adjacent intervals within a DataFrame.

    This function takes a DataFrame containing genomic interval information, formats
    it to conform with PyRanges input requirements, merges overlapping or adjacent
    intervals, and returns the processed intervals in the original DataFrame structure.

    Parameters:
        df (pd.DataFrame): A DataFrame containing the input intervals.

    Returns:
        pd.DataFrame: A DataFrame with merged genomic intervals, formatted back to
        the original structure with renamed column labels.
    """

    # Rename to PyRanges expected format
    pr_df = pyranges_column_name_change(df)

    # Get `pr_df` as a PyRanges obejct
    pr_df = pr.PyRanges(pr_df)

    # Merge (collapse overlapping/adjacent intervals)
    merged = pr_df.merge()

    # Convert back to DataFrame and rename
    result = merged.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result

