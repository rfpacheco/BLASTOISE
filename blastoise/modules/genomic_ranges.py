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

# noinspection PyPackageRequirements
import pandas as pd
# noinspection PyPackageRequirements
import pyranges as pr
from typing import Optional


# Column mapping constants
BLAST_TO_PYRANGES = {
    "sseqid": "Chromosome",
    "sstart": "Start",
    "send": "End"
}

PYRANGES_TO_BLAST = {
    "Chromosome": "sseqid",
    "Start": "sstart",
    "End": "send"
}


def to_pyranges_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename specific columns in a DataFrame to standardized names for PyRanges.

    This function converts column names from BLAST format ('sseqid', 'sstart', 'send')
    to PyRanges format ('Chromosome', 'Start', 'End') to ensure compatibility with
    PyRanges operations. The mapping is defined in the BLAST_TO_PYRANGES dictionary.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing genomic interval data with columns 'sseqid', 'sstart',
        and 'send'.

    Returns
    -------
    pd.DataFrame
        A new DataFrame with columns renamed according to PyRanges conventions.

    Raises
    ------
    KeyError
        If any of the required columns ('sseqid', 'sstart', 'send') are missing
        from the input DataFrame.

    See Also
    --------
    from_pyranges_format : Convert column names from PyRanges format back to BLAST format.
    """

    return df.rename(columns=BLAST_TO_PYRANGES)


def from_pyranges_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename PyRanges columns back to the original BLAST format.

    This function converts column names from PyRanges format ('Chromosome', 'Start', 'End')
    back to BLAST format ('sseqid', 'sstart', 'send') after PyRanges operations have been
    completed. The mapping is defined in the PYRANGES_TO_BLAST dictionary.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with PyRanges column names ('Chromosome', 'Start', 'End').

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns renamed back to BLAST format ('sseqid', 'sstart', 'send').

    Raises
    ------
    KeyError
        If any of the required columns ('Chromosome', 'Start', 'End') are missing
        from the input DataFrame.

    See Also
    --------
    to_pyranges_format : Convert column names from BLAST format to PyRanges format.
    """

    return df.rename(columns=PYRANGES_TO_BLAST)


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


def merge_intervals(
    df: pd.DataFrame,
    chr_col: str = "sseqid",
    start_col: str = "sstart",
    end_col: str = "send",
) -> pd.DataFrame:
    """
    Merge overlapping or adjacent genomic intervals within a DataFrame.

    This function takes a DataFrame containing genomic interval information, maps the
    provided columns to PyRanges format, merges overlapping or adjacent intervals, and
    returns the processed intervals using the same column names provided.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing genomic intervals to be merged. Must contain columns
        specified by `chr_col`, `start_col`, and `end_col`.
    chr_col : str, default="sseqid"
        Column name representing the chromosome/sequence identifier.
    start_col : str, default="sstart"
        Column name representing the start coordinate.
    end_col : str, default="send"
        Column name representing the end coordinate.

    Returns
    -------
    pd.DataFrame
        A DataFrame with merged genomic intervals, formatted back to the original
        column names given by `chr_col`, `start_col`, and `end_col`. If the input
        DataFrame is empty or no intervals can be merged, an empty DataFrame with
        these columns is returned.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the input DataFrame.
    ValueError
        If the input DataFrame cannot be converted to a PyRanges object.

    See Also
    --------
    get_interval_overlap : Find overlapping intervals between two DataFrames.
    get_merge_stranded : Merge intervals while respecting strand information.
    """

    required = {chr_col, start_col, end_col}  # Needed columns to be present
    missing = required - set(df.columns)  # Check which necessary columns are not present in the data frame
    if missing:
        raise KeyError(f"Missing required columns: {', '.join(sorted(missing))}")

    if df.empty:
        return pd.DataFrame(columns=[chr_col, start_col, end_col])

    # Map to PyRanges format
    pr_input = df[[chr_col, start_col, end_col]].rename(
        columns={chr_col: "Chromosome", start_col: "Start", end_col: "End"}
    )

    pr_df = pr.PyRanges(pr_input)
    merged = pr_df.merge()

    # Convert back to the original column names
    if hasattr(merged, "df"):
        out = merged.df.rename(
            columns={"Chromosome": chr_col, "Start": start_col, "End": end_col}
        )[[chr_col, start_col, end_col]]

        # Ensure integer type for coordinates, if possible
        for col in (start_col, end_col):
            if pd.api.types.is_numeric_dtype(out[col]):
                out[col] = out[col].astype(int)

        # Optional: sort for consistency
        out = out.sort_values(by=[chr_col, start_col], kind="mergesort").reset_index(drop=True)
        return out
    else:
        return pd.DataFrame(columns=[chr_col, start_col, end_col])


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
            merged_df = merge_intervals(strand_df)
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
    merged_data = merge_intervals(overlapping_data)
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
