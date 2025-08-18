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
from idlelib.debugger_r import DictProxy

# noinspection PyPackageRequirements
import pandas as pd
# noinspection PyPackageRequirements
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
    df_renamed = df.rename(columns=BLAST_TO_PYRANGES)
    
    # Convert strand notation if Strand column exists
    if 'Strand' in df_renamed.columns:
        # Map plus/minus to +/-
        strand_mapping = {'plus': '+', 'minus': '-'}
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
    df_renamed = df.rename(columns=PYRANGES_TO_BLAST)
    
    # Convert strand notation back if sstrand column exists
    if 'sstrand' in df_renamed.columns:
        # Map +/- back to plus/minus
        strand_mapping = {'+': 'plus', '-': 'minus'}
        df_renamed['sstrand'] = df_renamed['sstrand'].map(strand_mapping).fillna(df_renamed['sstrand'])
    
    return df_renamed


def get_interval_overlap(df: pd.DataFrame, interval_df: pd.DataFrame, invert: bool = False) -> pd.DataFrame:
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
    pr_df = pr.PyRanges(to_pyranges_format(df))
    pr_intervals = pr.PyRanges(to_pyranges_format(interval_df))

    # Perform an overlap operation
    result_pr = pr_df.overlap(pr_intervals, strandedness=False, invert=invert)

    # Convert back to DataFrame with original column names
    if hasattr(result_pr, 'df'):
        return from_pyranges_format(result_pr.df)
    else:
        return pd.DataFrame(columns=['sseqid', 'sstart', 'send'])


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


def get_merge_stranded(data_input: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping genomic ranges while respecting strand orientation.

    This function processes genomic data to merge overlapping ranges based on strand
    information using PyRanges. It handles 'plus' and 'minus' strands separately,
    ensuring that only intervals on the same strand are merged. This is critical for
    maintaining biological relevance, as features on opposite strands represent
    different genomic elements even if they overlap in coordinates.

    Parameters
    ----------
    data_input : pd.DataFrame
        Input data containing genomic coordinates and strand information.
        Must contain columns 'sstrand', 'sseqid', 'sstart', and 'send'.
        The 'sstrand' column should contain values 'plus' or 'minus'.

    Returns
    -------
    pd.DataFrame
        A processed DataFrame containing merged genomic ranges with columns
        'sseqid', 'sstart', 'send', 'sstrand', and 'len'. The 'len' column
        contains the length of each interval. If the input DataFrame is empty,
        an empty DataFrame with these columns is returned.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the input DataFrame.
    ValueError
        If the 'sstrand' column contains values other than 'plus' or 'minus'.

    See Also
    --------
    merge_intervals : Merge overlapping intervals without considering strand.
    """

    result_dfs = []

    # Process each strand separately
    for strand in ['plus', 'minus']:
        # Filter and sort data for the current strand
        strand_df = data_input[data_input['sstrand'] == strand].copy()

        if not strand_df.empty:
            # Sort and merge intervals
            strand_df = strand_df.sort_values(by=['sseqid', 'sstart'])
            merged_df = merge_overlapping_intervals(strand_df)
            merged_df['sstrand'] = strand
            result_dfs.append(merged_df)

    # Combine results from both strands
    if result_dfs:
        all_data = pd.concat(result_dfs, ignore_index=True)
        # Add a length column
        all_data['len'] = all_data['send'] - all_data['sstart'] + 1
        return all_data
    else:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=['sseqid', 'sstart', 'send', 'sstrand', 'len'])


def _print_dataset_stats(_description: str, data: pd.DataFrame, total: int, label: str) -> None:
    """
    Print statistics about a dataset, including count and percentage.

    This helper function calculates and prints statistics about a dataset,
    including the count of rows and the percentage relative to a total.
    It's used to provide informative output during the comparison of genomic
    datasets, helping users understand the overlap between different datasets.

    Parameters
    ----------
    _description : str
        _description of the data category being reported (e.g., "Coincidence",
        "Non-coincidence").
    data : pd.DataFrame
        The dataset to report statistics for. The function will count the
        number of rows in this DataFrame.
    total : int
        The total number of items for percentage calculation. This is typically
        the total number of rows in the original dataset.
    label : str
        Label for the data being reported (e.g., "New data in Previous data",
        "Previous data NOT in New data").

    Returns
    -------
    None
        This function only prints to standard output and doesn't return a value.

    Notes
    -----
    The function handles the case where the total is zero, avoiding division by zero
    errors by not calculating a percentage in that case.
    """

    count = data.shape[0]
    if total > 0:
        percentage = count / total * 100
        print(f"\t\t\t\t- {label}: {count}/{total} - {percentage:.2f}%")
    else:
        print(f"\t\t\t\t- {label}: {count}/{total}")


def compare_genomic_datasets(
        main_data: pd.DataFrame,
        data_for_contrast: pd.DataFrame,
        strand: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compare two genomic datasets to identify overlapping and non-overlapping regions.

    This function performs a comprehensive comparison between two genomic datasets
    to identify regions that overlap between them and regions that are unique to each.
    It uses PyRanges tools to efficiently find overlaps, sorts the data by genomic
    coordinates, and merges overlapping regions within each category. The function
    also prints detailed statistics about the comparison results.

    Parameters
    ----------
    main_data : pd.DataFrame
        The primary dataset containing genomic data to be analyzed. Must have
        columns 'sseqid', 'sstart', 'send', and 'sstrand'. This is typically
        the newer or query dataset in the comparison.
    data_for_contrast : pd.DataFrame
        The secondary dataset used for comparison against the main dataset. Must
        have the same column structure as main_data. This is typically the older
        or reference dataset in the comparison.
    strand : str
        The genomic strand orientation to assign to the overlapping regions.
        Should be either 'plus' or 'minus'.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        A tuple containing three DataFrames:
        - coincidence_data : pd.DataFrame
            Regions that overlap between both datasets after merging. Contains
            columns 'sseqid', 'sstart', 'send', and 'sstrand'.
        - new_data : pd.DataFrame
            Regions unique to `main_data` (not found in `data_for_contrast`).
            May be empty if all regions in main_data overlap with data_for_contrast.
        - only_in_contrast_data : pd.DataFrame
            Regions unique to `data_for_contrast` (not found in `main_data`).
            May be empty if all regions in data_for_contrast overlap with main_data.

    Raises
    ------
    KeyError
        If any of the required columns are missing from either input DataFrame.
    ValueError
        If the strand parameter is not 'plus' or 'minus'.

    Notes
    -----
    The function prints detailed statistics about the comparison results to standard
    output, including counts and percentages of overlapping and non-overlapping regions.

    See Also
    --------
    get_interval_overlap : Find overlapping intervals between two DataFrames.
    merge_intervals : Merge overlapping intervals within a DataFrame.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Prepare input data by sorting by chromosome and start position
    # -----------------------------------------------------------------------------
    # Sort both datasets by chromosome (sseqid) and start position for consistent processing
    main_data = main_data.sort_values(by=['sseqid', 'sstart'])
    data_for_contrast = data_for_contrast.sort_values(by=['sseqid', 'sstart'])

    # Store the original row counts for statistical reporting
    main_data_len = main_data.shape[0]
    data_contrast_len = data_for_contrast.shape[0]

    # -----------------------------------------------------------------------------
    # STEP 2: Identify overlapping regions between the two datasets
    # -----------------------------------------------------------------------------
    print("\n\t\t\t- Coincidence data:")

    # Find regions in main_data that overlap with data_for_contrast
    main_exists_in_contrast = get_interval_overlap(main_data, data_for_contrast, invert=False)
    _print_dataset_stats("Coincidence", main_exists_in_contrast, main_data_len, 
                        "New data in Previous data")

    # Find regions in data_for_contrast that overlap with main_data
    contrast_exists_in_main = get_interval_overlap(data_for_contrast, main_data, invert=False)
    _print_dataset_stats("Coincidence", contrast_exists_in_main, data_contrast_len, 
                        "Previous data in New data")

    # -----------------------------------------------------------------------------
    # STEP 3: Merge all overlapping regions into a single coherent set
    # -----------------------------------------------------------------------------
    # Combine overlapping regions from both directions (main→contrast and contrast→main)
    overlapping_data = pd.concat([main_exists_in_contrast, contrast_exists_in_main], ignore_index=True)

    # Merge any overlapping intervals within the combined set
    merged_data = merge_overlapping_intervals(overlapping_data)
    print(f"\t\t\t\t- Merged data: {merged_data.shape[0]}")

    # Create the final coincidence dataset with strand information
    coincidence_data = merged_data.copy()
    coincidence_data['sstrand'] = strand  # Assign the specified strand to all overlapping regions

    # -----------------------------------------------------------------------------
    # STEP 4: Identify non-overlapping (unique) regions in each dataset
    # -----------------------------------------------------------------------------
    print("\n\t\t\t- NOT coincidence data:")

    # Find regions in main_data that do NOT overlap with data_for_contrast (unique to main)
    main_not_in_contrast = get_interval_overlap(main_data, data_for_contrast, invert=True)
    _print_dataset_stats("Non-coincidence", main_not_in_contrast, main_data_len, 
                        "New data NOT in Previous data")

    # Find regions in data_for_contrast that do NOT overlap with main_data (unique to contrast)
    contrast_not_in_main = get_interval_overlap(data_for_contrast, main_data, invert=True)
    _print_dataset_stats("Non-coincidence", contrast_not_in_main, data_contrast_len, 
                        "Previous data NOT in New data")

    # -----------------------------------------------------------------------------
    # STEP 5: Prepare final return values, handling empty datasets gracefully
    # -----------------------------------------------------------------------------
    # For new_data and only_in_contrast_data, return empty DataFrame if no unique regions found
    new_data = main_not_in_contrast if not main_not_in_contrast.empty else pd.DataFrame()
    only_in_contrast_data = contrast_not_in_main if not contrast_not_in_main.empty else pd.DataFrame()

    # Return the three datasets: overlapping regions, unique to main, unique to contrast
    return coincidence_data, new_data, only_in_contrast_data
