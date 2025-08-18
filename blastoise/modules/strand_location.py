"""
BLASTOISE Module: Strand Orientation and Genomic Interval Processing
===================================================================

This module provides functionality for handling strand orientation and overlapping
genomic intervals in the BLASTOISE pipeline. It processes genomic coordinates to
ensure proper strand assignment and resolves overlapping sequences across different
strands.

The module contains several key functions:
1. Data matching and filtering functions (`match_data`, `match_data_and_remove`, 
   `match_data_and_set_false`) for comparing and processing genomic intervals
2. Overlapping data processing functions (`process_overlapping_data`, 
   `smart_merge_across_flips`) for handling complex overlapping scenarios
3. Strand orientation functions (`set_overlapping_status`, `set_strand_direction`) 
   for determining the correct strand assignment for sequences
4. Cleanup functions (`del_last_overlapping_elem`) for removing redundant overlapping
   elements after processing

These functions work together to ensure that genomic intervals are properly oriented
and that overlapping sequences are handled consistently, which is critical for
accurate identification of repetitive elements in the genome.

Author: R. Pacheco
"""

# Import needed modules
import pandas as pd
import os
import time
from joblib import Parallel, delayed

from typing import Hashable

from .genomic_ranges import get_merge_stranded
from .genomic_ranges import get_interval_overlap, merge_intervals
from blastoise.extra.utils.csv_to_gff import csv_to_gff


# noinspection DuplicatedCode
def match_data(data_input: pd.DataFrame, to_discard: pd.DataFrame) -> pd.DataFrame:
    """
    Match data from two DataFrames based on specific genomic coordinate columns.

    This function performs an inner join operation between the input DataFrame
    `data_input` and a subset of columns from the DataFrame `to_discard`. The
    matching is conducted over the genomic coordinate columns 'sseqid', 'sstart', 
    'send', and 'sstrand'. The resulting DataFrame contains rows where the values 
    in these columns align in both DataFrames.

    Parameters
    ----------
    data_input : pd.DataFrame
        The primary DataFrame to be matched against. Must contain the columns
        'sseqid', 'sstart', 'send', and 'sstrand'.
    to_discard : pd.DataFrame
        The DataFrame from which specific columns are used to find matches in 
        `data_input`. Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing rows from `data_input` that share the same values 
        in the specified columns with the corresponding rows in `to_discard`.

    Raises
    ------
    KeyError
        If any of the required columns are missing from either DataFrame.
    """

    matches = data_input.merge(
        to_discard[['sseqid', 'sstart', 'send', 'sstrand']],
        on=['sseqid', 'sstart', 'send', 'sstrand'],
        how='inner'
    )

    return matches


def match_data_and_remove(data_input: pd.DataFrame, to_discard: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows from a DataFrame based on matching genomic coordinates from another DataFrame.

    This function identifies rows from the given `data_input` DataFrame that matches rows in the 
    `to_discard` DataFrame based on genomic coordinate columns and removes these matched rows. 
    The matching is performed using the columns 'sseqid', 'sstart', 'send', and 'sstrand' 
    as the indices for comparison.

    Parameters
    ----------
    data_input : pd.DataFrame
        The input DataFrame from which rows should be evaluated and potentially removed.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.
    to_discard : pd.DataFrame
        A DataFrame containing rows to be matched and removed from the `data_input` DataFrame.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.

    Returns
    -------
    pd.DataFrame
        A copy of the `data_input` DataFrame with rows matching the criteria removed.

    Raises
    ------
    KeyError
        If any of the required columns are missing from either DataFrame.
    """

    matches = match_data(data_input, to_discard)

    # Remove only the matching data
    data_input = data_input.loc[
        ~data_input.set_index(['sseqid', 'sstart', 'send', 'sstrand']).index.isin(
            matches.set_index(['sseqid', 'sstart', 'send', 'sstrand']).index
        )
    ].copy()

    return data_input


def match_data_and_set_false(data_input: pd.DataFrame, to_discard: pd.DataFrame) -> None:
    """
    Mark matching rows in a DataFrame as False in the 'analyze' column.

    This function takes a primary dataset and a dataset of elements to be discarded, 
    then matches elements based on genomic coordinate columns. It uses these matches 
    to update the 'analyze' column in the primary dataset, marking the matched 
    elements as `False`.

    Parameters
    ----------
    data_input : pd.DataFrame
        The primary dataset to be analyzed and updated. Must contain the columns:
        'sseqid', 'sstart', 'send', 'sstrand', and 'analyze'.
    to_discard : pd.DataFrame
        The dataset containing elements that should be matched and marked as not 
        for analysis in the primary dataset. Must contain the columns: 'sseqid', 
        'sstart', 'send', and 'sstrand'.

    Returns
    -------
    None
        The function operates directly on the `data_input` DataFrame, updating its 
        'analyze' column in-place. No new object is returned.

    Raises
    ------
    KeyError
        If any of the required columns are missing from either DataFrame.
    """

    # Get any element in `data_input` with the same start and end coordinates as `to_discard` dataset
    matches = match_data(data_input, to_discard)

    # Use the index of the matches to update the 'analyze' column in data_input
    data_input.loc[
        data_input.set_index(['sseqid', 'sstart', 'send', 'sstrand']).index.isin(
            matches.set_index(['sseqid', 'sstart', 'send', 'sstrand']).index
        ),
        'analyze'
    ] = False


def process_overlapping_data(
        new_df: pd.DataFrame,
        row_df: pd.DataFrame,  # will have only 1 sequence
        idx: Hashable,
        all_og_inrange: pd.DataFrame,  # will have only 1 sequence
        all_elems_inrange: pd.DataFrame  # will have multiple sequences
) -> None:
    """
    Process overlapping genomic intervals based on strand alignment.

    This function modifies the provided DataFrame (`new_df`) by analyzing and adjusting 
    rows based on the comparison of overlapping sequences and their strands. It handles 
    two main cases:
    1. When the current row is on the same strand as the original sequence, it merges
       all overlapping sequences on that strand and updates the coordinates.
    2. When the current row is on a different strand, it removes overlapping sequences
       on the same strand as the current row.

    Parameters
    ----------
    new_df : pd.DataFrame
        The main DataFrame containing elements to be modified based on overlapping data.
        Must contain the columns 'sseqid', 'sstart', 'send', 'sstrand', and 'analyze'.
    row_df : pd.DataFrame
        A single-row DataFrame representing the sequence currently under analysis.
        Must contain the column 'sstrand'.
    idx : Hashable
        The index of the row in `new_df` that corresponds to `row_df`.
    all_og_inrange : pd.DataFrame
        A single-row DataFrame representing the original sequence for comparison.
        Must contain the column 'sstrand'.
    all_elems_inrange : pd.DataFrame
        A DataFrame containing multiple sequences that overlap with `row_df` or
        `all_og_inrange`. Must contain the column 'sstrand'.

    Returns
    -------
    None
        This function directly modifies the `new_df` without returning a new DataFrame.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the DataFrames.
    IndexError
        If `all_og_inrange` or `row_df` do not contain exactly one row.
    """

    same_strand = row_df.iloc[0]['sstrand'] == all_og_inrange.iloc[0]['sstrand']
    if same_strand:  # If row is in the same strand as the original data
        # Take from `all_elems_inrange` the ones in the same strand as `all_og_inrange`
        all_elems_inrange_same_strand = all_elems_inrange[
            all_elems_inrange['sstrand'] == all_og_inrange.iloc[0]['sstrand']
        ]
        # Merge 'all_elems_inrange_same_strand'. 'row_df' should be there as well, since it is in range as well
        merged_elem = merge_intervals(all_elems_inrange_same_strand)

        # Let's set false all elems in `all_elems_inrange` in the original `new_df`
        match_data_and_set_false(new_df, all_elems_inrange)

        # Change the 'row_df' to True in the original data
        new_df.loc[idx, 'analyze'] = True

        # Change the values for 'sstart' and 'send' coordinates
        new_df.loc[idx, ['sstart', 'send']] = merged_elem[['sstart', 'send']].values[0]
    else: # `row_df` and the og sequence are in different strands
        # From `all_elems_inrange` take the elements in the same strand as `row_df`
        all_elems_inrange_same_strand = all_elems_inrange[
            all_elems_inrange['sstrand'] == row_df.iloc[0]['sstrand']
            ]
        # And remove these elements from the `new_df`
        match_data_and_set_false(new_df, all_elems_inrange_same_strand)
