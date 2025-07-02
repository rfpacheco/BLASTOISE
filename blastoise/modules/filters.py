"""
BLASTOISE Module: Genomic Data Filtering and Masking
===================================================

This module provides functionality for filtering and masking genomic data in the BLASTOISE
pipeline. It handles the application of length filters, coordinate normalization, and
removal of regions that overlap with predefined masked zones.

The module contains two main functions:
1. `global_filters_main`: Applies length filtering and coordinate normalization to genomic data
2. `remove_masking_zone`: Removes regions from genomic data that overlap with specified masked zones

These functions work together to ensure that the genomic data used in the BLASTOISE pipeline
meets quality criteria and excludes regions that should not be considered in the analysis,
such as known repetitive elements or low-complexity regions.

Author: R. Pacheco
"""

import time
import pandas as pd

from .genomic_ranges import get_merge_stranded, get_interval_overlap
from .strand_location import match_data_and_remove


def global_filters_main(data_input: pd.DataFrame, min_length: int) -> pd.DataFrame:
    """
    Apply filtering and processing operations to genomic data.

    This function processes genomic data through a series of filtering and modification
    steps to ensure it meets quality criteria. It calculates sequence lengths, filters
    out sequences shorter than the specified minimum length, removes dashes from
    coordinates, and optionally merges overlapping intervals using PyRanges.

    The processing follows these main steps:
    1. Calculate sequence lengths based on start and end coordinates
    2. Filter out sequences shorter than the minimum length
    3. Remove dashes from coordinates for consistency
    4. Merge overlapping intervals using PyRanges (if data remains after filtering)

    Parameters
    ----------
    data_input : pd.DataFrame
        The input DataFrame containing genomic coordinates. Must include columns
        'send' and 'sstart' for calculating sequence lengths.
    min_length : int 
        The minimum acceptable sequence length. Sequences shorter than this
        value will be filtered out.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the filtered and processed genomic data. If the
        input data is empty or all sequences are filtered out, an empty DataFrame
        is returned.

    Raises
    ------
    KeyError
        If the required columns 'send' and 'sstart' are missing from the input DataFrame.
    TypeError
        If the column values cannot be converted to the required data types.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Calculate sequence lengths based on start and end coordinates
    # -----------------------------------------------------------------------------
    print("\t\t\t- Setting length")
    # Calculate the length of each sequence (inclusive of start and end positions)
    data_input["length"] = abs(data_input["send"] - data_input["sstart"]) + 1

    # -----------------------------------------------------------------------------
    # STEP 2: Filter out sequences shorter than the minimum length
    # -----------------------------------------------------------------------------
    print("\t\t\t- Filtering by length")
    # Keep only sequences that meet or exceed the minimum length requirement
    data_filtered = data_input[data_input["length"].astype(int) >= min_length]

    # -----------------------------------------------------------------------------
    # STEP 3: Remove dashes from coordinates for consistency
    # -----------------------------------------------------------------------------
    print("\t\t\t- Removing dashes from coordinates")
    # Replace any dash characters in the data with empty strings
    # This ensures consistent formatting of coordinate values
    data_filtered = data_filtered.apply(lambda x: x.replace("-", ""))

    # Create a copy of the filtered data to avoid modifying the original
    final_data = data_filtered.copy()

    # -----------------------------------------------------------------------------
    # STEP 4: Merge overlapping intervals using PyRanges (if data remains)
    # -----------------------------------------------------------------------------
    if not data_filtered.empty:  # Only process if there's data after filtering
        tic = time.perf_counter()
        print("\t\t\t- Processing with PyRange")
        # Use PyRanges to merge overlapping genomic intervals
        data_pyranges = get_merge_stranded(data_filtered)
        final_data = data_pyranges.copy()
        # Report execution time for performance monitoring
        toc = time.perf_counter()
        print(f"\t\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    # Return the filtered and processed data
    return final_data


def remove_masking_zone(data_input: pd.DataFrame, masking_data: pd.DataFrame) -> pd.DataFrame:
    """
    Remove genomic regions that overlap with specified masked zones.

    This function filters genomic data by removing regions that overlap with predefined
    masked zones. It uses PyRanges tools to efficiently identify overlapping regions,
    then removes these regions from the input data. This is particularly useful for
    excluding known repetitive elements, low-complexity regions, or other areas that
    should not be considered in downstream analysis.

    The function performs two main steps:
    1. Identify regions in the input data that overlap with masking data
    2. Remove these overlapping regions from the input data

    Parameters
    ----------
    data_input : pd.DataFrame
        The DataFrame containing genomic intervals to be filtered. Must contain
        columns that define genomic coordinates (e.g., 'sseqid', 'sstart', 'send').
    masking_data : pd.DataFrame
        The DataFrame containing masked zones to exclude. Must have the same
        coordinate column structure as data_input.

    Returns
    -------
    pd.DataFrame
        A filtered DataFrame with overlapping regions removed. If no regions
        overlap with the masking data, the original input data is returned
        (as a copy to avoid modifying the original).

    Raises
    ------
    KeyError
        If the required coordinate columns are missing from either DataFrame.

    Notes
    -----
    The function uses the `get_interval_overlap` function to identify overlapping
    regions and the `match_data_and_remove` function to remove them from the input data.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Identify regions that overlap with the masking data
    # -----------------------------------------------------------------------------
    # Find all regions in data_input that overlap with any region in masking_data
    overlaps_with_masking = get_interval_overlap(data_input, masking_data)

    # -----------------------------------------------------------------------------
    # STEP 2: Remove overlapping regions from the input data
    # -----------------------------------------------------------------------------
    # Create a copy of the input data to avoid modifying the original
    data_masked = data_input.copy()

    # Only perform removal if overlapping regions were found
    if not overlaps_with_masking.empty:
        # Remove all overlapping regions from the input data
        data_masked = match_data_and_remove(data_input, overlaps_with_masking)

    # Return the filtered data with masked zones removed
    return data_masked
