"""
BLASTOISE Module: Genomic Range Operations
==========================================

This module provides functionality for handling genomic range operations in the BLASTOISE
pipeline. It leverages the PyRanges library to efficiently process and manipulate genomic
intervals, enabling operations such as finding overlaps, merging intervals, and comparing
datasets.

The module contains several key functions:
1. Format conversion functions (`to_pyranges_format`, `from_pyranges_format`) for
   translating between BLAST and PyRanges coordinate systems
2. Interval operation functions (`get_interval_overlap`, `merge_intervals`) for
   finding overlapping regions and merging adjacent intervals
3. Strand-specific processing functions (`get_merge_stranded`) for handling genomic
   data with strand orientation
4. Dataset comparison functions (`compare_genomic_datasets`) for identifying shared
   and unique regions between genomic datasets

These functions work together to provide a comprehensive toolkit for genomic interval
manipulation, which is essential for the identification and analysis of repetitive
elements in genomic data.

Author: R. Pacheco
"""
import pandas as pd
import pyranges as pr
from typing import Optional, List, Dict

# Update the mapping dictionaries to include strand information
BLAST_TO_PYRANGES = {
    'sseqid': 'Chromosome',
    'sstart': 'Start',
    'send': 'End',
    'sstrand': 'Strand'
}

PYRANGES_TO_BLAST = {
    'Chromosome': 'sseqid',
    'Start': 'sstart',
    'End': 'send',
    'Strand': 'sstrand'
}

def to_pyranges_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename specific columns in a DataFrame to standardized names for PyRanges.

    This function converts column names from BLAST format ('sseqid', 'sstart', 'send', 'sstrand')
    to PyRanges format ('Chromosome', 'Start', 'End', 'Strand') to ensure compatibility with
    PyRanges operations. The mapping is defined in the BLAST_TO_PYRANGES dictionary.
    Additionally, it converts strand notation from 'plus'/'minus' to '+'/'-' format.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing genomic interval data with columns 'sseqid', 'sstart',
        'send', and optionally 'sstrand'.

    Returns
    -------
    pd.DataFrame
        A new DataFrame with columns renamed according to PyRanges conventions
        and strand values converted to PyRanges format.

    Raises
    ------
    KeyError
        If any of the required columns ('sseqid', 'sstart', 'send') are missing
        from the input DataFrame.

    See Also
    --------
    from_pyranges_format : Convert column names from PyRanges format back to BLAST format.
    """
    # Create a copy to avoid modifying the original DataFrame
    df_renamed: pd.DataFrame = df.rename(columns=BLAST_TO_PYRANGES)
    
    # Convert strand notation if Strand column exists
    if 'Strand' in df_renamed.columns:
        # Map plus/minus to +/-
        strand_mapping: Dict[str, str] = {'plus': '+', 'minus': '-'}
        df_renamed['Strand'] = df_renamed['Strand'].map(strand_mapping).fillna(df_renamed['Strand'])
    
    return df_renamed


def from_pyranges_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename specific columns in a DataFrame from PyRanges format back to BLAST format.

    This function converts column names from PyRanges format ('Chromosome', 'Start', 'End', 'Strand')
    back to BLAST format ('sseqid', 'sstart', 'send', 'sstrand') and converts strand notation
    from '+'/'-' back to 'plus'/'minus' format.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with PyRanges-formatted column names.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns renamed back to BLAST format and strand values
        converted to BLAST format.

    See Also
    --------
    to_pyranges_format : Convert column names from BLAST format to PyRanges format.
    """
    # Create a copy to avoid modifying the original DataFrame
    df_renamed: pd.DataFrame = df.rename(columns=PYRANGES_TO_BLAST)
    
    # Convert strand notation back if sstrand column exists
    if 'sstrand' in df_renamed.columns:
        # Map +/- back to plus/minus
        strand_mapping: Dict[str, str] = {'+': 'plus', '-': 'minus'}
        df_renamed['sstrand'] = df_renamed['sstrand'].map(strand_mapping).fillna(df_renamed['sstrand'])
    
    return df_renamed


def fetch_overlapping_intervals(df: pd.DataFrame, interval_df: pd.DataFrame, invert: bool = False) -> pd.DataFrame:
    """
    Find intervals between two DataFrames that either overlap or don't overlap.

    This function takes two DataFrames containing genomic intervals, converts them
    into PyRanges objects, and identifies either the overlapping or non-overlapping
    intervals between them based on the 'invert' parameter. It's a key function for
    comparing genomic datasets and identifying regions of interest.

    Parameters
    ----------
    df : pd.DataFrame
        The main DataFrame containing genomic intervals for which overlap
        will be calculated. Must contain columns 'sseqid', 'sstart', and 'send'.
    interval_df : pd.DataFrame
        The DataFrame containing intervals to check for overlap with the main
        DataFrame. Must contain columns 'sseqid', 'sstart', and 'send'.
    invert : bool, default=False
        If True, returns intervals in df that do not overlap with interval_df.
        If False, returns intervals in df that do overlap with interval_df.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the requested genomic intervals with standardized
        column names ('sseqid', 'sstart', 'send'). If no intervals meet the criteria,
        an empty DataFrame with these columns is returned.

    Raises
    ------
    KeyError
        If any of the required columns are missing from either input DataFrame.
    ValueError
        If the input DataFrames cannot be converted to PyRanges objects.

    See Also
    --------
    merge_intervals : Merge overlapping or adjacent intervals within a DataFrame.
    """

    # Convert dataframes to PyRanges format
    pr_df: pr.PyRanges = pr.PyRanges(to_pyranges_format(df))
    pr_intervals: pr.PyRanges = pr.PyRanges(to_pyranges_format(interval_df))

    # Perform an overlap operation
    result_pr: pr.PyRanges = pr_df.overlap(pr_intervals, strandedness=False, invert=invert)

    # Get returning data frame
    out: pd.DataFrame = result_pr.df
    out = from_pyranges_format(out)
    return out


def merge_overlapping_intervals(
        df: pd.DataFrame,
        chr_col: str = 'sseqid',
        start_col: str = 'sstart',
        end_col: str = 'send',
        strand_col: str = 'sstrand',
        strand: bool = False,
) -> pd.DataFrame:
    """
    Merges overlapping intervals represented in a DataFrame and optionally considers the strand of the intervals.

    This function takes a pandas DataFrame containing interval data and merges overlapping intervals.
    The user can provide custom column names for chromosome, start, end, and strand. It also supports merging
    intervals by strand when the `strand` parameter is set to True. The result is a new DataFrame with merged intervals,
    sorted by chromosome and start positions.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing genomic interval data. It must include columns for chromosome, start, and end coordinates.
    chr_col : str, optional
        Column name in `df` representing the chromosome. Default is 'sseqid'.
    start_col : str, optional
        Column name in `df` representing the start coordinate of intervals. Default is 'sstart'.
    end_col : str, optional
        Column name in `df` representing the end coordinate of intervals. Default is 'send'.
    strand_col : str, optional
        Column name in `df` representing the strand information. Default is 'sstrand'.
    strand : bool, optional
        Boolean indicating whether to merge intervals while considering strand information. If True, the function will require
        `strand_col` to be specified and consider strand-specific merging. Default is False.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with merged intervals. The columns in the returned DataFrame correspond to the ones specified in the
        input (`chr_col`, `start_col`, `end_col`, and if applicable `strand_col`). The intervals are merged and sorted by
        chromosome and start coordinate.
    """
    # Return empty DataFrame with appropriate columns if input is empty
    if df.empty:
        columns: List[str] = [chr_col, start_col, end_col]
        if strand and strand_col:
            columns.append(strand_col)
        return pd.DataFrame(columns=columns)

    # Map to PyRanges format
    pr_input: pd.DataFrame
    if strand:
        # Gets PyRanges format with strand-specific column
        # Convert strand format from 'plus'/'minus' to '+'/'−'
        pr_input = df[[chr_col, start_col, end_col, strand_col]].copy()
        pr_input[strand_col] = pr_input[strand_col].map({'plus': '+', 'minus': '-'})
        pr_input = pr_input.rename(
            columns={chr_col: "Chromosome", start_col: "Start", end_col: "End", strand_col: "Strand"}
        )
    else:
        # Gets PyRanges format without strand-specific column
        pr_input = df[[chr_col, start_col, end_col]].rename(
            columns={chr_col: "Chromosome", start_col: "Start", end_col: "End"}
        )

    pr_df: pr.PyRanges = pr.PyRanges(pr_input)
    # If strand == True, will merge stranded
    # If strand == False, will ignore strands
    merged: pr.PyRanges = pr_df.merge(strand=strand)

    # Convert back to the original column names
    rename_map: Dict[str, str] = {"Chromosome": chr_col, "Start": start_col, "End": end_col}

    if strand:
        rename_map["Strand"] = strand_col

    out: pd.DataFrame = merged.df.rename(columns=rename_map)

    # Convert strand format back from '+'/'−' to 'plus'/'minus'
    if strand and strand_col:
        out[strand_col] = out[strand_col].map({'+': 'plus', '-': 'minus'})

    out = out.sort_values(by=[chr_col, start_col]).reset_index(drop=True)
    return out


def get_overlapping_info(df: pd.DataFrame, interval_df: Optional[pd.DataFrame] = None) -> dict[str, pd.DataFrame]:
    """
    Identifies and separates overlapping genomic intervals by strand relationship in the
    provided DataFrame(s) using PyRanges.

    This function takes a primary DataFrame containing genomic interval data and an optional
    interval DataFrame. It computes the overlaps between the primary DataFrame and the interval
    DataFrame (if provided) based on strand relationships. Overlaps are classified into two
    categories: overlaps on the "same strand" or on the "opposite strand".

    When interval_df is None (self-overlap), the function automatically excludes self-matches
    to avoid reporting every interval as overlapping with itself.

    The function returns a dictionary containing the overlapping intervals classified by these
    two strand relationships.

    Parameters
    ----------
    df : pd.DataFrame
        Primary DataFrame containing genomic interval data. This DataFrame must adhere
        to PyRanges-compatible format including required columns such as 'Chromosome',
        'Start', 'End', and 'Strand'.
    interval_df : pd.DataFrame or None, optional
        An optional DataFrame containing genomic intervals to compare against. If None,
        overlaps will be computed within the primary DataFrame itself.

    Returns
    -------
    dict of str, pd.DataFrame
        A dictionary containing two keys:
        - "same_strand": pd.DataFrame
            A DataFrame of intervals overlapping on the same strand.
        - "opposite_strand": pd.DataFrame
            A DataFrame of intervals overlapping on the opposite strand.
        Each DataFrame in the returned dictionary will have columns for sequences, start,
        and end positions (e.g., 'sseqid', 'sstart', 'send'), or will be empty if no
        overlaps are detected.
    """
    # Convert input DataFrame to PyRanges format
    pr_df = pr.PyRanges(to_pyranges_format(df))
    # Use input DataFrame if interval_df is None, otherwise convert interval_df to PyRanges
    pr_other = pr_df if interval_df is None else pr.PyRanges(to_pyranges_format(interval_df))
    
    # Flag to check if we're doing self-overlap (comparing df against itself)
    is_self_overlap = interval_df is None

    # Find overlapping intervals on the same strand
    same_pr = pr_df.overlap(pr_other, strandedness="same", invert=False)
    # Find overlapping intervals on opposite strands
    opp_pr = pr_df.overlap(pr_other, strandedness="opposite", invert=False)

    # Convert results back to DataFrames, use empty DataFrame if no overlaps found
    same_df = from_pyranges_format(same_pr.df) if hasattr(same_pr, "df") else \
        pd.DataFrame(columns=["sseqid", "sstart", "send"])
    opp_df = from_pyranges_format(opp_pr.df) if hasattr(opp_pr, "df") else \
        pd.DataFrame(columns=["sseqid", "sstart", "send"])

    # If this is a self-overlap, remove self-matches (intervals overlapping with themselves)
    if is_self_overlap:
        same_df = _remove_self_matches(same_df, df)
        opp_df = _remove_self_matches(opp_df, df)

    # Return dictionary with overlaps separated by strand relationship
    return {"same_strand": same_df, "opposite_strand": opp_df}


def _remove_self_matches(overlap_df: pd.DataFrame, original_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove self-matches from overlap results where an interval overlaps with itself.
    
    This function identifies exact matches between the overlap results and the original
    DataFrame to remove trivial self-overlaps that occur when comparing a DataFrame 
    against itself.
    
    Parameters
    ----------
    overlap_df : pd.DataFrame
        DataFrame containing overlap results from PyRanges.
    original_df : pd.DataFrame
        The original DataFrame that was used for self-overlap detection.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with self-matches removed.
    """
    if overlap_df.empty or original_df.empty:
        return overlap_df
    
    # Create a set of tuples representing the original intervals for fast lookup
    original_intervals = set()
    for _, row in original_df.iterrows():
        # Use the key columns that uniquely identify an interval
        key_cols = ['sseqid', 'sstart', 'send', 'sstrand']
        if all(col in row.index for col in key_cols):
            interval_key = tuple(row[col] for col in key_cols)
            original_intervals.add(interval_key)
    
    # Filter out self-matches from overlap results
    if not original_intervals:
        return overlap_df
    
    # Create mask to identify self-matches
    self_match_mask = pd.Series(False, index=overlap_df.index)
    
    for idx, row in overlap_df.iterrows():
        key_cols = ['sseqid', 'sstart', 'send', 'sstrand']
        if all(col in row.index for col in key_cols):
            interval_key = tuple(row[col] for col in key_cols)
            if interval_key in original_intervals:
                self_match_mask.loc[idx] = True
    
    # Return DataFrame with self-matches removed
    return overlap_df[~self_match_mask].copy()
