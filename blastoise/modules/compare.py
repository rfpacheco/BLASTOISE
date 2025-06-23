"""
BLASTOISE Module: Genomic Dataset Comparison
===========================================

This module provides functionality for comparing genomic datasets in the BLASTOISE
pipeline. It handles the comparison of genomic features between different datasets,
identifying common and unique regions across both plus and minus strands.

The module contains one main function:
1. `compare_main`: Performs comparative analysis between two genomic datasets,
   processing plus and minus strands separately and combining the results.

This function serves as a critical component in the BLASTOISE pipeline for identifying
similarities and differences between genomic datasets, which is essential for tracking
the evolution of repetitive element discovery across iterations.

Author: R. Pacheco
Version: 0.4.2
License: MIT
"""

import pandas as pd

from modules.genomic_ranges import compare_genomic_datasets
from modules.files_manager import columns_to_numeric


def compare_main(
        newest_data: pd.DataFrame,
        contrast_data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Perform comparative analysis between two genomic datasets across both DNA strands.

    This function analyzes two genomic datasets to identify common and unique features
    between them. It processes plus and minus strands separately using PyRanges tools
    for efficient overlap analysis, then combines the results. The function handles
    empty datasets gracefully and ensures proper numeric type conversion for all
    coordinate columns.

    The analysis follows these main steps:
    1. Separate data by strand (plus/minus)
    2. Process each strand independently using compare_genomic_datasets
    3. Convert coordinate columns to numeric types
    4. Combine results from both strands

    Parameters
    ----------
    newest_data : pd.DataFrame
        A DataFrame of the newest data to analyze, containing genomic feature information 
        such as strand and coordinates. Must contain columns 'sseqid', 'sstart', 'send',
        and 'sstrand'.
    contrast_data : pd.DataFrame
        A DataFrame of contrast data to compare against the newest data, with similar 
        genomic feature information. Must have the same column structure as newest_data.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        A tuple containing three DataFrames:
        - coincidence_data: Genomic features common between `newest_data` and
          `contrast_data`.
        - new_data: Genomic features exclusive to the `newest_data`.
        - old_data_exclusive: Genomic features exclusive to the `contrast_data`.

    Notes
    -----
    Each returned DataFrame may be empty if no features meet the corresponding criteria.
    All coordinate columns are converted to numeric types to ensure consistency.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Separate data by strand (plus/minus)
    # -----------------------------------------------------------------------------
    # Filter the newest data by strand
    newest_data_plus = newest_data[newest_data["sstrand"] == "plus"].copy()
    newest_data_minus = newest_data[newest_data["sstrand"] == "minus"].copy()

    # Filter the contrast data by strand
    contrast_data_plus = contrast_data[contrast_data["sstrand"] == "plus"].copy()
    contrast_data_minus = contrast_data[contrast_data["sstrand"] == "minus"].copy()

    # -----------------------------------------------------------------------------
    # STEP 2: Process plus strand data
    # -----------------------------------------------------------------------------
    print("")
    print("\t\t- '+' strand:")

    # Compare genomic datasets for plus strand
    coincidence_plus, new_data_plus, old_data_exclusive_plus = compare_genomic_datasets(
        newest_data_plus, contrast_data_plus, "plus"
    )

    # Convert coordinate columns to numeric type for each result dataset
    if not new_data_plus.empty:  # If the data frame has lines
        new_data_plus = columns_to_numeric(new_data_plus)  # Convert columns to numeric
    if not coincidence_plus.empty:  # If the data frame has lines
        coincidence_plus = columns_to_numeric(coincidence_plus)  # Convert columns to numeric
    if not old_data_exclusive_plus.empty:  # If the data frame has lines
        old_data_exclusive_plus = columns_to_numeric(old_data_exclusive_plus)  # Convert columns to numeric
    # -----------------------------------------------------------------------------
    # STEP 3: Process minus strand data
    # -----------------------------------------------------------------------------
    print("")
    print("\t\t- '-' strand:")

    # Compare genomic datasets for minus strand
    coincidence_minus, new_data_minus, old_data_exclusive_minus = compare_genomic_datasets(
        newest_data_minus, contrast_data_minus, "minus"
    )

    # Convert coordinate columns to numeric type for each result dataset
    if not new_data_minus.empty:  # If the data frame is not empty
        new_data_minus = columns_to_numeric(new_data_minus)
    if not coincidence_minus.empty:  # If the data frame is not empty
        coincidence_minus = columns_to_numeric(coincidence_minus)
    if not old_data_exclusive_minus.empty:  # If the data frame is not empty
        old_data_exclusive_minus = columns_to_numeric(old_data_exclusive_minus)
    # -----------------------------------------------------------------------------
    # STEP 4: Combine results from both strands
    # -----------------------------------------------------------------------------
    # Combine new data from plus and minus strands
    if not new_data_plus.empty or not new_data_minus.empty:
        # Concatenate non-empty DataFrames from both strands
        new_data = pd.concat([new_data_plus, new_data_minus], ignore_index=True)
    else:  # If both DataFrames are empty
        # Create an empty DataFrame with the expected structure
        new_data = pd.DataFrame()

    # Combine coincidence data from plus and minus strands
    if not coincidence_plus.empty or not coincidence_minus.empty:
        # Concatenate non-empty DataFrames from both strands
        coincidence_data = pd.concat([coincidence_plus, coincidence_minus], ignore_index=True)
    else:  # If both DataFrames are empty
        # Create an empty DataFrame with the expected structure
        coincidence_data = pd.DataFrame()

    # Combine old data exclusive from plus and minus strands
    if not old_data_exclusive_plus.empty or not old_data_exclusive_minus.empty:
        # Concatenate non-empty DataFrames from both strands
        old_data_exclusive = pd.concat([old_data_exclusive_plus, old_data_exclusive_minus], ignore_index=True)
    else:  # If both DataFrames are empty
        # Create an empty DataFrame with the expected structure
        old_data_exclusive = pd.DataFrame()
    # -----------------------------------------------------------------------------
    return coincidence_data, new_data, old_data_exclusive
