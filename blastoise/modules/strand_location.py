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
Version: 0.4.2
License: MIT
"""

# Import needed modules
import pandas as pd
import os
import time
from joblib import Parallel, delayed

from typing import Hashable

from modules.genomic_ranges import get_merge_stranded
from modules.genomic_ranges import get_interval_overlap, merge_intervals
from extra.utils.csv_to_gff import csv_to_gff


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


def smart_merge_across_flips(
    all_og_inrange: pd.DataFrame,
    all_elems_inrange: pd.DataFrame,
    strand_col: str = "sstrand",
    start_col: str = "sstart",
    end_col: str = "send",
) -> pd.DataFrame:
    """
    Intelligently merge overlapping genomic intervals across strand flips.

    This function processes genomic intervals by grouping them into "strand blocks" 
    determined by consecutive strands and applying a sliding window approach to 
    determine merging criteria. It handles three main cases:

    1. Outer blocks with matching strands: When two blocks with the same strand 
       orientation are separated by a block with a different orientation, and the 
       outer blocks overlap, they are merged while the middle block is removed.

    2. Adjacent blocks with overlap: When adjacent blocks overlap but don't meet 
       the criteria for case 1, the overlapping elements are removed to maintain 
       strand specificity.

    3. Default case: Blocks that don't match any special criteria are kept intact.

    Parameters
    ----------
    all_og_inrange : pd.DataFrame
        A DataFrame containing the original intervals with their corresponding strand 
        and positions. Must contain the columns specified by `strand_col`, `start_col`, 
        and `end_col`.
    all_elems_inrange : pd.DataFrame
        A DataFrame containing elements that may potentially overlap or coincide with 
        intervals in `all_og_inrange`. Must have the same column structure as 
        `all_og_inrange`.
    strand_col : str, optional
        Column name representing the strand direction ("plus" or "minus"), 
        defaults to "sstrand".
    start_col : str, optional
        Column name representing the start position of intervals, 
        defaults to "sstart".
    end_col : str, optional
        Column name representing the end position of intervals, 
        defaults to "send".

    Returns
    -------
    pd.DataFrame
        A DataFrame that retains strand-specific intervals, intelligently merges 
        overlapping blocks across strand flips, and progressively processes all 
        interval blocks.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the DataFrames.
    """

    # ==================================================================
    # 0) Build the strand-block structure
    # ==================================================================
    # Detect points where strand changes between consecutive rows
    change_points = (
            all_og_inrange[strand_col] != all_og_inrange[strand_col].shift()
    )
    # Create blocks of consecutive rows with the same strand by cumulative sum of change points
    # noinspection PyUnresolvedReferences
    block_id = change_points.cumsum()
    all_og_inrange = all_og_inrange.assign(_block_id=block_id)

    # Using the original `all_in_range` data, get the minimum and maximum coordinates for each block
    block_summary = (
        all_og_inrange.groupby("_block_id")
        .agg(block_start=(start_col, "min"),  # Get min start coord
             block_end=(end_col, "max"),  # Get max end coord
             strand=(strand_col, "first"))  # Get strand of block
    )

    # ==================================================================
    # 1) Sliding window on blocks: b0, b1, b2
    # ==================================================================
    # Initialize list to store processed block ids and other control variables
    ids = block_summary.index.to_list()  # Get a list of block IDs
    n_blocks = len(ids)  # Total number of blocks
    final_chunks = []  # Store final processed chunks
    skip_next_blocks = set()  # Track blocks to skip

    for idx, b0 in enumerate(ids):
        # Skip if the block was already processed 
        if b0 in skip_next_blocks:
            continue

        # Get next blocks (b1,b2) if they exist, otherwise None
        b1 = ids[idx + 1] if idx + 1 < n_blocks else None
        b2 = ids[idx + 2] if idx + 2 < n_blocks else None

        # =================================================================
        # Case: We have 3 consecutive blocks (b0,b1,b2)
        # =================================================================
        # If b1 and b2 exists
        if b1 is not None and b2 is not None:
            # Get strand orientation for each block
            s0 = block_summary.loc[b0, "strand"]
            s1 = block_summary.loc[b1, "strand"]
            s2 = block_summary.loc[b2, "strand"]

            # Check if outer blocks (b0,b2) overlap
            ## Get overlapping elements for b0 and b2 blocks 
            elems_of_b0 = get_interval_overlap(all_elems_inrange,
                                                   all_og_inrange.loc[all_og_inrange["_block_id"] == b0])
            elems_of_b2 = get_interval_overlap(all_elems_inrange,
                                                   all_og_inrange.loc[all_og_inrange["_block_id"] == b2])
            ## Check if they share overlapping elements
            elems_of_b0_vs_b2 = get_interval_overlap(elems_of_b0, elems_of_b2, invert=False)
            outer_overlap = False
            if not elems_of_b0_vs_b2.empty:  # If it has rows
                outer_overlap = True  # Then the b1 and b1 share overlapping elements

            # Check if adjacent blocks (b0,b1) overlap  
            ## Get overlapping elements for b1 block
            elems_of_b1 = get_interval_overlap(all_elems_inrange,
                                                   all_og_inrange.loc[all_og_inrange["_block_id"] == b1])
            ## Check if b0,b1 share overlapping elements
            elems_of_b0_vs_b1 = get_interval_overlap(elems_of_b0, elems_of_b1, invert=False)
            inner_overlap = False
            if not elems_of_b0_vs_b1.empty:  # If it has rows
                inner_overlap = True  # Then the adyacent b0 and b1 share overlapping elements

            # Rule 1: Outer blocks have same strand & overlap
            if s0 == s2 and s0 != s1 and outer_overlap:
                # When the outer blocks (b0, b2) have the same strand orientation 
                # and different from the middle block (b1), and they overlap

                # Combine overlapping elements from blocks b0 & b2
                outer_rows = pd.concat([elems_of_b0, elems_of_b2])

                # Get original elements from b0 and b2 blocks
                og_b1_and_b2 = pd.concat([
                    all_og_inrange.loc[all_og_inrange["_block_id"] == b0],
                    all_og_inrange.loc[all_og_inrange["_block_id"] == b2]
                ])

                # Merge all elements and remove duplicates 
                all_elems = pd.concat([outer_rows, og_b1_and_b2]).drop_duplicates()
                all_elems.sort_values(start_col)

                # Merge overlapping intervals into single ranges
                collapsed_outer = merge_intervals(all_elems)
                collapsed_outer["sstrand"] = s0  # Add strand value

                # Add merged blocks to the final output
                final_chunks.append(collapsed_outer)

                # Mark the middle and last blocks to skip in future iterations
                skip_next_blocks.update({b1, b2})

                continue

            # Rule 2: Only adjacent blocks (b0,b1) overlap
            # Rule 2: Only adjacent blocks (b0,b1) overlap
            if inner_overlap:
                # Remove overlapping elements between b0,b1 that connect the two blocks

                # Get elements from b1 that overlap with b0
                elems_of_b1_vs_b0 = get_interval_overlap(elems_of_b1, elems_of_b0, invert=False)

                # Combine and deduplicate overlapping elements
                elems_to_remove = pd.concat(
                    [elems_of_b1_vs_b0, elems_of_b0_vs_b1]
                ).drop_duplicates()
                elems_to_remove.sort_values(start_col)

                # Remove overlapping elements
                cleaned = match_data_and_remove(all_elems_inrange, elems_to_remove)

                # Add cleaned data to the final output
                final_chunks.append(cleaned)

                # Skip processing b0 again since it's handled
                skip_next_blocks.add(b0)

                continue

        # Default: Keep the block untouched if no rule applies
        untouched = all_og_inrange[all_og_inrange["_block_id"] == b0]
        final_chunks.append(untouched)

    # =================================================================
    # Build final output - merge chunks & cleanup
    # =================================================================
    result = (pd.concat(final_chunks, ignore_index=True, copy=False)
              .drop(columns="_block_id", errors="ignore")  # Safe drop temp column
              .sort_values(start_col)
              .reset_index(drop=True))

    return result

# noinspection PyUnusedLocal
def _set_overlapping_status_single(
        chrom: str,
        new_df_chr: pd.DataFrame,
        og_df_chr: pd.DataFrame,
        run_phase: int
) -> pd.DataFrame:
    """
    Process overlapping genomic intervals for a single chromosome.

    This function analyzes and processes overlapping intervals within two DataFrames 
    that represent genomic data for a single chromosome. It iterates through each row 
    in the new DataFrame, identifies overlapping sequences in the original DataFrame, 
    and applies appropriate processing based on strand orientation and overlap patterns.

    The function handles several cases:
    1. Single strand in original data: Merges overlapping sequences on the same strand
    2. Two strands with simple overlap: Removes overlapping connections between strands
    3. Complex strand patterns: Applies specialized processing for interlaced strands

    Parameters
    ----------
    chrom : str
        The name of the chromosome being processed.
    new_df_chr : pd.DataFrame
        A DataFrame representing new genomic intervals for the chromosome.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.
        An 'analyze' column will be added if not present.
    og_df_chr : pd.DataFrame
        A DataFrame representing the original genomic intervals for the chromosome.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.
    run_phase : int
        The current phase/iteration of the analysis.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the resolved intervals after processing overlaps.
        Only rows with 'analyze' set to True are included in the result.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the DataFrames.
    """
    # To avoid taking in data already analyzed, insert a boolean True column. NOTE: important

    new_df_chr['analyze'] = True
    for idx, row in new_df_chr.iterrows():
        # Skip already processed elements to avoid processing them again
        if not new_df_chr.loc[idx, 'analyze']:
            continue

        # Get row as a pd.DataFrame and not a pd.Series
        row = new_df_chr.loc[idx:idx, :].copy()

        # Get elements from `og_df_chr` that overlap with `row`.
        # NOTO: 'vs' will be used instead of 'overlap'
        og_vs_row = get_interval_overlap(og_df_chr, row, invert=False)

        # And now, all the elements that overlap with `og_vs_row`. # NOTE: take only True values in `new_df_chr`
        new_elems_in_og_vs_row = get_interval_overlap(
            new_df_chr[new_df_chr['analyze'] == True], og_vs_row, invert=False
        )

        # Now get all new elements that overlap with 'new_elems_in_og_vs_row'.
        # NOTE: take only True values in `new_df_chr`
        new_elems_inrange_of_og = get_interval_overlap(
            new_df_chr[new_df_chr['analyze'] == True], new_elems_in_og_vs_row, invert=False
        )

        # And get all original elements in the whole range of `all_elems_in_range`
        og_inrange = get_interval_overlap(og_df_chr, new_elems_inrange_of_og, invert=False)

        # Now, get all elems that overlap all this `og_inrange`
        all_elems_vs_og_inrange = get_interval_overlap(new_df_chr[new_df_chr['analyze'] == True], og_inrange, invert=False)

        # If these new elems in `all_elems_vs_og_inrange` overlap some other new element coming from a little far
        # away original element, we need to detect them.
        all_og_inrange = get_interval_overlap(og_df_chr, all_elems_vs_og_inrange, invert=False)
        all_elems_inrange = get_interval_overlap(new_df_chr[new_df_chr['analyze'] == True], all_og_inrange, invert=False)

        # Let's count how many original elements are in "minus" and "plus" strand. There could be 4 cases
        ## 1) There's only 1 original element in range.
        ## 2) There are 2 or more original elements in the range. Same strand.
        ## 3) There are 2 or more original elements in the range. Different strands
        how_many_og = {}
        for og_strand in all_og_inrange['sstrand'].unique():
            how_many_og[og_strand] = all_og_inrange[all_og_inrange['sstrand'] == og_strand].shape[0]

        # Case 1)
        # If there's only 1 element, be it 'minus' or 'plus'
        if len(how_many_og) == 1: # Only one 'strand' is present in `og_inrange`
            # In this case, it doesn't matter if its 1 sequences in `og_inrange` or > 2 sequences. The merged will
            # be implemented the same way
            process_overlapping_data(
                new_df_chr,
                row,
                idx,
                all_og_inrange,
                all_elems_inrange
            )
        else: # When there are hits in different strand sequences from `original_og_inrange`
            if sum(how_many_og.values()) == 2:  # The normal case is when 1 sequence is in 'plus' and the other in 'minus'
                # In this case the first step to avoid overlaps is to remove all the sequences in `all_elems_inrange`
                # of the original element A that overlaps all the elements in `all_elems_inrange` of element B
                og_a = og_inrange.loc[0:0, :]
                og_b = og_inrange.loc[1:1, :]
                elems_of_a = get_interval_overlap(all_elems_inrange, og_a, invert=False)
                elems_of_b = get_interval_overlap(all_elems_inrange, og_b, invert=False)
                elems_to_remove_a = get_interval_overlap(elems_of_a, elems_of_b, invert=False)
                elems_to_remove_b = get_interval_overlap(elems_of_b, elems_of_a, invert=False)
                elems_to_remove = pd.concat(
                    [elems_to_remove_a, elems_to_remove_b]
                ).drop_duplicates().sort_values(['sstart', 'send'])
                match_data_and_set_false(new_df_chr, elems_to_remove)
                is_row_removed = match_data(row, elems_to_remove)
                if not is_row_removed.empty:
                    continue
                else:
                    # Now that the connection between the 2 `original_og_inrange` is removed. The rest will behave like
                    # as if it were only one `original_og_inrange`
                    og_vs_row = get_interval_overlap(og_inrange, row, invert=False) # Selects original data that overlaps with row
                    all_elems_vs_og = get_interval_overlap(new_df_chr[new_df_chr['analyze'] == True], og_vs_row, invert=False)
                    process_overlapping_data(
                        new_df_chr,
                        row,
                        idx,
                        og_vs_row,
                        all_elems_vs_og
                    )
            else:

                # These are not normal cases. For example, it is when in the original data is
                # 'minus' -- 'plus' -- 'minus', or 'plus' -- 'minus' -- 'plus', or 'plus' -- 'plus' -- 'minus'.
                ## First, in the `all_og_inrange` detect the strand flip
                change_points = (all_og_inrange['sstrand'] != all_og_inrange['sstrand'].shift())

                # Now check the blocks by number
                # noinspection PyUnresolvedReferences
                block_id = change_points.cumsum()

                if block_id.unique().shape[0] == 2:
                    # if there are only 2 blocks, it means a structure without interlaps like:
                    # 'plus' -- 'plus' -- 'minus' or 'plus' -- 'minus' -- 'minus'
                    ## Take from `all_elems_inrange` the elems with each strand
                    elems_in_plus = all_elems_inrange[all_elems_inrange['sstrand'] == 'plus']
                    elems_in_minus = all_elems_inrange[all_elems_inrange['sstrand'] == 'minus']

                    # Take the elems that overlap each one
                    elems_in_plus_vs_minus = get_interval_overlap(elems_in_plus, elems_in_minus, invert=False)
                    elems_in_minus_vs_plus = get_interval_overlap(elems_in_minus, elems_in_plus, invert=False)
                    elems_to_remove = pd.concat(
                        [elems_in_plus_vs_minus, elems_in_minus_vs_plus]
                    ).sort_values(['sstart', 'send'])
                    match_data_and_set_false(new_df_chr, elems_to_remove)
                    is_row_removed = match_data(row, elems_to_remove)
                    if not is_row_removed.empty:
                        continue
                    else:
                        elems_in_row = get_interval_overlap(new_df_chr[new_df_chr['analyze'] == True], row, invert=False)
                        process_overlapping_data(
                            new_df_chr,
                            row,
                            idx,
                            og_vs_row,
                            elems_in_row
                        )
                else:
                    # This one is harder, since there are interlaps for sure, something like:
                    # 'plus' -- 'minus' -- 'plus'
                    # Create one-row summary per block to make the window logic easier
                    ## First make all these `all_elems_inrange` in the original `df` as False
                    match_data_and_set_false(new_df_chr, all_elems_inrange)
                    all_elems_inrange_resolved = smart_merge_across_flips(all_og_inrange, all_elems_inrange)
                    idx_locator = idx
                    for _, elem in all_elems_inrange_resolved.iterrows():
                        # Replace each idx in `new_df_chr` for each element in `all_elems_inrange_resolved`
                        new_df_chr.loc[idx_locator, ['sstart', 'send']] = elem[['sstart', 'send']]
                        new_df_chr.loc[idx_locator, 'sstrand'] = elem['sstrand']
                        new_df_chr.loc[idx_locator, 'analyze'] = True
                        idx_locator += 1

    final_data = new_df_chr[new_df_chr['analyze'] == True].copy()

    return final_data


def set_overlapping_status(
        new_df: pd.DataFrame,
        og_df: pd.DataFrame,
        run_phase: int,
        n_jobs: int = -1
) -> pd.DataFrame:
    """
    Process overlapping genomic intervals in parallel across chromosomes.

    This function is a parallel wrapper around the `_set_overlapping_status_single` 
    function. It splits the input DataFrames by chromosome, processes each chromosome 
    in parallel using the specified number of jobs, and then combines the results.

    The parallelization strategy follows these steps:
    1. Split both DataFrames by chromosome
    2. Process each chromosome independently in parallel
    3. Combine the results and maintain the original order

    Parameters
    ----------
    new_df : pd.DataFrame
        A DataFrame containing new genomic intervals to analyze.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.
    og_df : pd.DataFrame
        A DataFrame containing original/reference genomic intervals.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.
    run_phase : int
        The current phase/iteration of the analysis.
    n_jobs : int, default -1
        Number of parallel jobs to run. -1 means using all available processors,
        1 means fallback to the original single-process execution.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame containing the processed intervals from all chromosomes,
        sorted by 'sseqid' and 'sstart'.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the DataFrames.
    """

    # Fast exit: keep the exact behavior if the caller explicitly disables parallelism.
    if n_jobs == 1:
        return _set_overlapping_status_single("ALL", new_df, og_df, run_phase)

    # ==================================================================
    # 1. Split the two dataframes by chromosome
    # ==================================================================
    new_groups = {c: df for c, df in new_df.groupby("sseqid", sort=False)}
    og_groups  = {c: df for c, df in og_df.groupby("sseqid", sort=False)}

    # Use only chromosomes that are present in *new* – everything else would produce empty output anyway.
    chromosomes = list(new_groups)

    # ==================================================================
    # 2. Run one worker per chromosome
    # ==================================================================
    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_set_overlapping_status_single)(
            chrom,
            new_groups[chrom],
            og_groups.get(chrom, og_df.iloc[0:0],),  # empty DF if missing # TODO: check slicing
            run_phase
        )
        for chrom in chromosomes
    )

    # ==================================================================
    # 3. Concatenate the per-chromosome outputs and keep original order
    # ==================================================================
    combined_df = pd.concat(results, ignore_index=True)
    combined_df.sort_values(by=['sseqid', 'sstart'], inplace=True)
    return combined_df

def set_strand_direction(
        data_input: pd.DataFrame,
        run_phase: int,
        folder_path: str,
        n_jobs: int = -1
) -> pd.DataFrame:
    """
    Determine the correct strand orientation for genomic sequences.

    This function analyzes and processes genomic sequence data to determine the correct 
    strand orientation for each sequence. It handles two main categories of sequences:

    1. New sequences: Sequences that don't overlap with any sequence from previous iterations
    2. Overlapping sequences: Sequences that overlap with sequences from previous iterations

    For each category, the function processes plus and minus strands separately, resolves
    overlaps between strands, and ensures data consistency across iterations.

    The processing follows these main steps:
    1. Separate original coordinates and new coordinates
    2. Identify new elements and overlapping elements
    3. Process overlapping elements using parallel processing
    4. Handle strand-specific processing for new elements
    5. Combine the results into a single consistent DataFrame

    Parameters
    ----------
    data_input : pd.DataFrame
        Input DataFrame containing sequence data with columns:
        'og_sseqid', 'og_sstart', 'og_send', 'og_sstrand', 'sseqid', 'sstart', 
        'send', and 'sstrand'.
    run_phase : int
        The current phase/iteration of the analysis.
    folder_path : str
        Path to the folder where temporary files will be stored.
    n_jobs : int, default -1
        Number of parallel jobs to run. -1 means using all available processors.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing processed sequences with correct strand orientation,
        sorted by 'sseqid' and 'sstart'.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the DataFrame.
    ValueError
        If the folder_path is invalid or cannot be created.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Prepare original and new coordinate data
    # -----------------------------------------------------------------------------
    # Extract and format original coordinates before extension
    og_data = data_input[['og_sseqid', 'og_sstart', 'og_send', 'og_sstrand']].copy()
    og_data.columns = ['sseqid', 'sstart', 'send', 'sstrand']  # Rename columns for consistency
    og_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    og_data.drop_duplicates(inplace=True)

    # Extract and format new coordinates
    new_data = data_input[['sseqid', 'sstart', 'send', 'sstrand']].copy()
    new_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    new_data.drop_duplicates(inplace=True)

    # -----------------------------------------------------------------------------
    # STEP 2: Separate new elements from overlapping elements
    # -----------------------------------------------------------------------------
    # Identify elements that don't overlap with original data (truly new elements)
    new_elems = get_interval_overlap(new_data, og_data, invert=True)
    print(f"\t\t\t- New elements: {new_elems.shape[0]}")

    # Identify elements that overlap with original data
    overlapping_elems = pd.merge(
        new_data, new_elems,
        on=['sseqid', 'sstart', 'send', 'sstrand'],
        how='left', indicator=True
    ).query('_merge == "left_only"').drop('_merge', axis=1)
    print(f"\t\t\t- Overlapping elements: {overlapping_elems.shape[0]}")

    # -----------------------------------------------------------------------------
    # STEP 3: Process overlapping elements
    # -----------------------------------------------------------------------------
    # Initialize strand-specific DataFrames
    og_and_overlap_elems_plus = pd.DataFrame()  # Will hold plus strand overlapping elements
    og_and_overlap_elems_minus = pd.DataFrame()  # Will hold minus strand overlapping elements

    if not overlapping_elems.empty:
        tic = time.perf_counter()
        print(f"\t\t\t- Checking strand orientation in overlapping elements:")

        # Save intermediate files for debugging and visualization
        save_folder = os.path.join(folder_path, 'RUNS')
        os.makedirs(save_folder, exist_ok=True)
        overlapping_elems.to_csv(
            os.path.join(save_folder, f"run_{run_phase - 1}_new_df.csv"), index=False
        )
        csv_to_gff(
            os.path.join(save_folder, f"run_{run_phase - 1}_new_df.csv")
        )
        og_data.to_csv(
            os.path.join(save_folder, f"run_{run_phase - 1}_og_df.csv"), index=False
        )
        csv_to_gff(
            os.path.join(save_folder, f"run_{run_phase - 1}_og_df.csv")
        )

        # Process overlapping elements using parallel processing
        overlapping_elems = set_overlapping_status(overlapping_elems, og_data, run_phase, n_jobs=n_jobs)
        overlapping_elems = get_merge_stranded(overlapping_elems)  # Merge overlapping intervals
        toc = time.perf_counter()
        print(f"\t\t\t\t- Execution time: {toc - tic:0.2f} seconds")

        # Separate overlapping elements by strand
        og_and_overlap_elems_plus = overlapping_elems[overlapping_elems['sstrand'] == 'plus'].copy()
        og_and_overlap_elems_minus = overlapping_elems[overlapping_elems['sstrand'] == 'minus'].copy()

    # -----------------------------------------------------------------------------
    # STEP 4: Process new elements
    # -----------------------------------------------------------------------------
    if not new_elems.empty:
        # Separate new elements by strand
        new_elems_plus = new_elems[new_elems['sstrand'] == 'plus'].copy()
        new_elems_minus = new_elems[new_elems['sstrand'] == 'minus'].copy()

        # -----------------------------------------------------------------------------
        # STEP 4.1: Remove overlaps between strands within new elements
        # -----------------------------------------------------------------------------
        # Identify elements that overlap between strands
        new_elems_plus_in_minus = get_interval_overlap(
            new_elems_plus, new_elems_minus, invert=False
        )
        new_elems_minus_in_plus = get_interval_overlap(
            new_elems_minus, new_elems_plus, invert=False
        )

        # Remove overlapping elements from each strand
        if not new_elems_plus_in_minus.empty:
            new_elems_plus = match_data_and_remove(new_elems_plus, new_elems_plus_in_minus)

        if not new_elems_minus_in_plus.empty:
            new_elems_minus = match_data_and_remove(new_elems_minus, new_elems_minus_in_plus)

        # Merge overlapping intervals within each strand
        new_elems_plus = get_merge_stranded(new_elems_plus)
        new_elems_minus = get_merge_stranded(new_elems_minus)

        # -----------------------------------------------------------------------------
        # STEP 4.2: Remove overlaps between new elements and overlapping elements
        # -----------------------------------------------------------------------------
        # Handle plus strand new elements that overlap with minus strand overlapping elements
        if not og_and_overlap_elems_minus.empty:
            new_elems_plus_vs_og_and_overlap_minus = get_interval_overlap(
                new_elems_plus, og_and_overlap_elems_minus, invert=False
            )
            if not new_elems_plus_vs_og_and_overlap_minus.empty:
                new_elems_plus = match_data_and_remove(new_elems_plus, new_elems_plus_vs_og_and_overlap_minus)

        # Handle minus strand new elements that overlap with plus strand overlapping elements
        if not og_and_overlap_elems_plus.empty:
            new_elems_minus_vs_og_and_overlap_plus = get_interval_overlap(
                new_elems_minus, og_and_overlap_elems_plus, invert=False
            )
            if not new_elems_minus_vs_og_and_overlap_plus.empty:
                new_elems_minus = match_data_and_remove(new_elems_minus, new_elems_minus_vs_og_and_overlap_plus)

        # Combine plus and minus strand new elements
        new_elems = pd.concat([new_elems_plus, new_elems_minus])
        new_elems.sort_values(by=['sseqid', 'sstart'], inplace=True)

    # -----------------------------------------------------------------------------
    # STEP 5: Combine results and return
    # -----------------------------------------------------------------------------
    # Combine new elements and processed overlapping elements
    result = pd.concat([new_elems, overlapping_elems])
    result.sort_values(by=['sseqid', 'sstart'], inplace=True)

    return result


def del_last_overlapping_elem(last_run_elems: pd.DataFrame) -> pd.DataFrame:
    """
    Remove redundant overlapping elements between different strands.

    This function performs a final cleanup step after strand processing to ensure
    that no redundant overlapping elements remain between different strands. When
    sequences on different strands (plus and minus) overlap, this function determines
    which sequence to keep based on length - the longer sequence is retained while
    the shorter one is removed.

    This is particularly important after operations like `smart_merge_across_flips`,
    which may generate sequences that unify blocks b0 and b2 while leaving block b1
    on another strand. In such cases, b1 should be removed to avoid redundancy.

    Parameters
    ----------
    last_run_elems : pd.DataFrame
        ADataFrame containing genomic elements discovered in the current iteration.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.

    Returns
    -------
    pd.DataFrame
        A cleaned DataFrame where, in cases of overlapping elements between strands,
        the smaller element has been removed. If no overlapping elements are found,
        the original DataFrame is returned unchanged.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the DataFrame.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Separate elements by strand
    # -----------------------------------------------------------------------------
    # Divide input elements into plus and minus strand DataFrames
    last_run_elems_plus = last_run_elems[last_run_elems['sstrand'] == 'plus'].copy()
    last_run_elems_minus = last_run_elems[last_run_elems['sstrand'] == 'minus'].copy()

    # -----------------------------------------------------------------------------
    # STEP 2: Identify overlapping elements between strands
    # -----------------------------------------------------------------------------
    # Find plus strand elements that overlap with minus strand elements
    last_run_elems_plus_vs_minus = get_interval_overlap(last_run_elems_plus, last_run_elems_minus, invert=False)
    # Find minus strand elements that overlap with plus strand elements
    last_run_elems_minus_vs_plus = get_interval_overlap(last_run_elems_minus, last_run_elems_plus, invert=False)

    # -----------------------------------------------------------------------------
    # STEP 3: Process overlapping elements if any exist
    # -----------------------------------------------------------------------------
    if not last_run_elems_plus_vs_minus.empty:
        # Initialize empty DataFrame to collect elements that should be removed
        elems_to_remove = pd.DataFrame()
        # Use plus strand overlapping elements as the basis for analysis
        data_to_analyze = last_run_elems_plus_vs_minus.copy()

        # -----------------------------------------------------------------------------
        # STEP 3.1: Compare each overlapping element and determine which to keep
        # -----------------------------------------------------------------------------
        for idx, row in data_to_analyze.iterrows():
            # Calculate the length of the current plus strand element
            row_len = row.send - row.sstart + 1

            # Find minus strand elements that overlap with this specific plus strand element
            overlapping_with_row = get_interval_overlap(
                last_run_elems_minus_vs_plus,
                # Convert the row to a DataFrame for compatibility with get_interval_overlap
                last_run_elems_plus_vs_minus.iloc[idx:idx+1, :],
                invert=False
            )

            # -----------------------------------------------------------------------------
            # STEP 3.2: Compare lengths and mark shorter elements for removal
            # -----------------------------------------------------------------------------
            for idx2, elem in overlapping_with_row.iterrows():
                # Calculate the length of the current minus strand element
                elem_len = elem.send - elem.sstart + 1

                # Keep the longer element and remove the shorter one
                # Note: It's assumed that elem_len and row_len cannot be equal
                if elem_len < row_len:
                    # If a minus strand element is shorter, mark it for removal
                    elems_to_remove = pd.concat(
                        [elems_to_remove, overlapping_with_row.iloc[idx2:idx2+1, :]],
                    )
                else:
                    # If a plus strand element is shorter, mark it for removal
                    elems_to_remove = pd.concat(
                        [elems_to_remove, data_to_analyze.iloc[idx:idx+1, :]]
                    )

        # -----------------------------------------------------------------------------
        # STEP 4: Remove the identified elements from the original DataFrame
        # -----------------------------------------------------------------------------
        last_run_elems = match_data_and_remove(last_run_elems, elems_to_remove)

    # Return the cleaned DataFrame with redundant overlapping elements removed
    return last_run_elems
