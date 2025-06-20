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
    Matches data from two DataFrame objects based on specific columns and returns the result.

    This function performs an inner join operation between the input DataFrame
    `data_input` and a subset of columns from the DataFrame `to_discard`. The
    matching is conducted over the columns 'sseqid', 'sstart', 'send', and
    'sstrand'. The resulting DataFrame contains rows where the values in these
    columns align in both DataFrames.

    Parameters:
    ----------
    data_input : pandas.DataFrame
        The primary DataFrame to be matched against.
    to_discard : pandas.DataFrame
        The DataFrame from which specific columns are used to find matches in `data_input`.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing rows from `data_input` that share the same values in the specified columns with the
        corresponding rows in `to_discard`.
    """
    matches = data_input.merge(
        to_discard[['sseqid', 'sstart', 'send', 'sstrand']],
        on=['sseqid', 'sstart', 'send', 'sstrand'],
        how='inner'
    )

    return matches


def match_data_and_remove(data_input: pd.DataFrame, to_discard: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows from a DataFrame based on matching criteria from another DataFrame.

    This function identifies rows from the given `data_input` DataFrame that matches the rows in the `to_discard`
    DataFrame based on specific indexing criteria and removes the matched rows. The matching is performed using the
    columns 'sseqid', 'sstart', 'send' and 'sstrand' as the indices for comparison.

    Parameters:
    -----------
    data_input: pandas.DataFrame
        The input DataFrame from which rows should be evaluated and potentially removed.
    to_discard: pandas.DataFrame
        A DataFrame containing rows to be matched and removed from the `data_input` DataFrame based on specified
        criteria.

    Returns:
    --------
    pandas.DataFrame
        A copy of the `data_input` DataFrame with rows matching the criteria removed.
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
    Matches data based on specific conditions and updates a column in the input dataset.

    This function takes a primary dataset and a dataset of elements to be discarded, then matches elements based on
    specific coordinate-based keys. It uses these matches to update a boolean column in the primary dataset, marking
    the matched elements as `False` in the 'analyze' column.

    Parameters:
    -----------
    data_input : pandas.DataFrame
        The primary dataset to be analyzed and updated. It must contain the columns: 'sseqid', 'sstart', 'send', and
        'sstrand'.

    to_discard : pandas.DataFrame
        The dataset containing elements that should be matched and marked as not for analysis in the primary dataset.
        Similar to `data_input`, it must also have the columns: 'sseqid', 'sstart', 'send', and 'sstrand'.

    Returns:
    --------
    None
        The function operates directly on the `data_input` DataFrame, updating its 'analyze' column in-place. No new
        object is returned.
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
        all_elems_inrange: pd.DataFrame  # will have multipel sequences
) -> None:
    """
    Processes overlapping data between sequences in a DataFrame to handle intervals based on strand alignment.
    This function modifies the provided DataFrame (`new_df`) by analyzing and adjusting certain rows based on 
    the comparison of overlapping sequences and their strands, as well as merging intervals when applicable.

    Parameters:
    -----------
    new_df: pandas.DataFrame
        The main DataFrame containing elements that need to be modified based on overlapping data. Modifications
        include updating rows to mark intervals merged or analyzed.
    row_df: pandas.DataFrame
        A single-row DataFrame representing the sequence currently under analysis. This DataFrame has one sequence that
        will be compared against overlapping sequences.
    idx (Hashable):
        The index of the row in the `new_df` that corresponds to the `row_df`.
    all_og_inrange: pandas.DataFrame
        A single-row DataFrame representing the original sequence against which comparisons of overlapping sequences
        will be made. It contains one sequence.
    all_elems_inrange: pandas.DataFrame
        A DataFrame containing multiple sequences/rows that are in the overlapping range with the `row_df` or
        `all_og_inrange`. This serves as the main pool of sequences to determine strand alignment and merging.

    Returns:
    -------
    None: This function directly modifies the `new_df` without returning a new DataFrame.
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
    Merges overlapping intervals across strand flips intelligently while retaining strand-specific intervals that do
    not overlap under specific conditions.

    This function processes genomic intervals, represented in the form of data frames, by grouping them into
    "strand blocks" determined by consecutive strands and sliding a window over these blocks to decide merging
    criteria. It handles three main cases: merging outer blocks with matching strands, handling overlapping
    blocks, and keeping untouched blocks intact. Results in a cleaned or merged data frame of intervals.

    Parameters:
    -----------
    all_og_inrange: pd.DataFrame
        Data frame containing the original intervals in range with their corresponding strand and positions.
    all_elems_inrange: pd.DataFrame
        Data frame containing elements that may potentially overlap or coincide with intervals in `all_og_inrange`.
    strand_col: str, optional
        Column name in `all_og_inrange` representing the strand direction ("+" or "-") of intervals, defaults to
        "sstrand".
    start_col: str, optional
        Column name in `all_og_inrange` representing the start position of intervals, defaults to "sstart".
    end_col: str, optional
        Column name in `all_og_inrange` representing the end position of intervals, defaults to "send".

    Returns:
    --------
    pd.DataFrame:
        A data frame that retains strand-specific intervals, intelligently merges overlapping blocks across strand
        flips, and progressively processes all interval blocks.
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
    Set the overlapping status for single chromosome data.

    This function analyzes and processes overlapping intervals within two dataframes that represent genomic data for a
    single chromosome. It ensures that overlapping entries are handled correctly based on specific rules and conditions,
    updating flags and making adjustments to resolve overlaps.

    Parameters:
    ----------
    chrom: str
        The name of the chromosome being processed.
    new_df_chr: pd.DataFrame
        A dataframe representing new genomic intervals for the chromosome. It must have the column 'analyze'.
    og_df_chr: pd.DataFrame
        A dataframe representing the original genomic intervals for the chromosome.
    run_phase: int
        An integer representing the phase of the analysis.

    Returns:
    --------
    pd.DataFrame
        A dataframe containing the resolved intervals after processing overlaps. The returned dataframe includes
        new intervals with updated 'analyze' flags and resolved overlaps.

    Raises:
    None
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
    Parallel wrapper around the original algorithm.

    Parameters
    ----------
    new_df : DataFrame
        Data to analyze.  Must contain the column ``'sseqid'``.
    og_df  : DataFrame
        Reference data.  Must contain the column ``'sseqid'``.
    run_phase : int
        Iteration counter
    n_jobs : int, default ``-1``
        Number of processes Joblib should spawn.  (``-1`` ⇒ use all cores, ``1`` ⇒ fallback to the original
        single-process execution.)
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
        folder_path: str
) -> pd.DataFrame:
    """
     Analyzes and processes genomic sequence data to determine the correct strand orientation for each sequence.
     It handles both completely new sequences and sequences that overlap with the 'n-1' iteration result, merging and 
     extending where appropriate, ensuring data consistency between strands.

    Parameters:
    -----------
    data_input: pandas.DataFrame
        Input DataFrame containing sequence data with columns such as 'og_sseqid', 'og_sstart', 'og_send',
        'og_sstrand', 'sseqid', 'sstart', 'send', and 'sstrand'.
    run_phase: int
        Indicates the run step iteratoin.
    folder_path: int
        Path to the folder where temporary files are stored.

    Returns:
    --------
    pd.DataFrame:
        A DataFrame with processed data that includes new elements and overlapping elements with orientation
        properly set.
    """

    # Save the original coordinates before the extension in `og_data`
    og_data = data_input[['og_sseqid', 'og_sstart', 'og_send', 'og_sstrand']].copy()
    og_data.columns = ['sseqid', 'sstart', 'send', 'sstrand'] # Change column names
    og_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    og_data.drop_duplicates(inplace=True)

    # Save all the new elements discovered in `new_data`
    new_data = data_input[['sseqid', 'sstart', 'send', 'sstrand']].copy()
    new_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    new_data.drop_duplicates(inplace=True)

    # Get the `new_data` that doesn't overlap with `og_data`
    # `new_elems` are considered new elements. Let's save them correctly. From 'new_data', extract the elements
    # that have the same 'sseqid, 'sstart', 'send' and 'sstrand' as in `new_elems`
    new_elems = get_interval_overlap(new_data, og_data, invert=True)
    print(f"\t\t\t- New elements: {new_elems.shape[0]}")

    # And the rest elements, which for sure, overlap.
    ## Keep the data from `new_data` that don't connect with the `new_elems`
    overlapping_elems = pd.merge(
        new_data, new_elems,
        on=['sseqid', 'sstart', 'send', 'sstrand'],
        how='left', indicator=True
    ).query('_merge == "left_only"').drop('_merge', axis=1)
    print(f"\t\t\t- Overlapping elements: {overlapping_elems.shape[0]}")

    # Now we need to check if the overlapping is in the same strand as the original data
    og_and_overlap_elems_plus = pd.DataFrame() # Initialization
    og_and_overlap_elems_minus = pd.DataFrame() # Initialization
    if not overlapping_elems.empty: # If it has rows
        tic = time.perf_counter()
        print(f"\t\t\t- Checking strand orientation in overlapping elements:")
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
        # Set the number of parallel jobs for multiprocessing
        multiprocessing_jobs = -1  # Use all available CPU cores
        # Process overlapping elements using parallel processing
        overlapping_elems = set_overlapping_status(overlapping_elems, og_data, run_phase, n_jobs=multiprocessing_jobs)
        overlapping_elems = get_merge_stranded(overlapping_elems) # Merge the data
        toc = time.perf_counter()
        print(f"\t\t\t\t- Execution time: {toc - tic:0.2f} seconds")

        # Now dive the data in 'minus' and plus' strand if there are rows with 'plus' or 'minus'
        og_and_overlap_elems_plus = overlapping_elems[overlapping_elems['sstrand'] == 'plus'].copy()
        og_and_overlap_elems_minus = overlapping_elems[overlapping_elems['sstrand'] == 'minus'].copy()

    if not new_elems.empty: # If it has rows
        # Now dive the data in 'minus' and plus' strand
        new_elems_plus = new_elems[new_elems['sstrand'] == 'plus'].copy()
        new_elems_minus = new_elems[new_elems['sstrand'] == 'minus'].copy()

        # Select elements that overlap between each strand in the 'new_elems'
        new_elems_plus_in_minus = get_interval_overlap(
            new_elems_plus, new_elems_minus, invert=False
        )
        new_elems_minus_in_plus = get_interval_overlap(
            new_elems_minus, new_elems_plus, invert=False
        )

        # And remove them
        if not new_elems_plus_in_minus.empty:
            new_elems_plus = match_data_and_remove(new_elems_plus, new_elems_plus_in_minus)

        if not new_elems_minus_in_plus.empty:
            new_elems_minus = match_data_and_remove(new_elems_minus, new_elems_minus_in_plus)

        # Not let's merge the results
        new_elems_plus = get_merge_stranded(new_elems_plus)
        new_elems_minus = get_merge_stranded(new_elems_minus)

        # Remove the elements that overlap with the original+overlapping sequences in the opposite strand
        if not og_and_overlap_elems_minus.empty:  # If it has rows
            # Get `new_elems_plus` that overlap in `og_and_overlap_elems_minus`
            new_elems_plus_vs_og_and_overlap_minus = get_interval_overlap(
                new_elems_plus, og_and_overlap_elems_minus, invert=False
            )
            if not new_elems_plus_vs_og_and_overlap_minus.empty: # If it has rows
                # Remove from `new_elems_plus` the elements from `new_elems_plus_vs_og_and_overlap_minus`
                new_elems_plus = match_data_and_remove(new_elems_plus, new_elems_plus_vs_og_and_overlap_minus)

        if not og_and_overlap_elems_plus.empty:  # If it has rows
            # Get `new_elems_minus` that overlap in `og_and_overlap_elems_plus`
            new_elems_minus_vs_og_and_overlap_plus = get_interval_overlap(
                new_elems_minus, og_and_overlap_elems_plus, invert=False
            )
            if not new_elems_minus_vs_og_and_overlap_plus.empty: # If it has rows
                # Remove from `new_elems_minus` the elements from `new_elems_minus_vs_og_and_overlap_plus`
                new_elems_minus = match_data_and_remove(new_elems_minus, new_elems_minus_vs_og_and_overlap_plus)

        # Join again in new_elems and remove tmp files
        new_elems = pd.concat([new_elems_plus, new_elems_minus])
        new_elems.sort_values(by=['sseqid', 'sstart'], inplace=True)

    # Combine new elements and processed leftover elements
    result = pd.concat([new_elems, overlapping_elems])
    result.sort_values(by=['sseqid', 'sstart'], inplace=True)

    return result


def del_last_overlapping_elem(last_run_elems: pd.DataFrame) -> pd.DataFrame:
    """
    After completing `smart_merge_across_flips`, sometimes the output will generate a sequence that unifies the blocks
    b0 and b2, leaving the bloc b1 in another strand. In this case, b1 should be removed

    Parameters:
    -----------
    last_run_elems: pd.DataFrame
        A dataframe with new elements discovered in the current iteration.

    Returns:
    --------
    pd.DataFrame:
        In case there are overlapping elements between the strands, the small one will be removed. If there are not,
        the original table will be returned
    """
    # Divide `last_run_elems` into its strands
    last_run_elems_plus = last_run_elems[last_run_elems['sstrand'] == 'plus'].copy()
    last_run_elems_minus = last_run_elems[last_run_elems['sstrand'] == 'minus'].copy()

    # Check the overlaps between each other
    last_run_elems_plus_vs_minus = get_interval_overlap(last_run_elems_plus, last_run_elems_minus, invert=False)
    last_run_elems_minus_vs_plus = get_interval_overlap(last_run_elems_minus, last_run_elems_plus, invert=False)

    # If there's data, iterate from one of the data frames
    if not last_run_elems_plus_vs_minus.empty:
        elems_to_remove = pd.DataFrame()
        data_to_analyze = last_run_elems_plus_vs_minus.copy()
        for idx, row in data_to_analyze.iterrows():
            # Check where it overlaps against `last_run_elems_minus_vs_plus`
            row_len = row.send - row.sstart + 1

            # Check for where it overlaps
            overlapping_with_row = get_interval_overlap(
                last_run_elems_minus_vs_plus,
                last_run_elems_plus_vs_minus.iloc[idx:idx+1, :],  # Get the pd.DataFrame version
                invert=False
            )

            # Check which element is bigger between row and the elems in `overlapping_with_row`
            for _, elem in overlapping_with_row.iterrows():
                elem_len = elem.send - elem.sstart + 1
                # It's not possible for elem_len == row_len
                if elem_len < row_len:
                    elems_to_remove = pd.concat([elems_to_remove, elem_len])
                else:
                    elems_to_remove = pd.concat([elems_to_remove, row])

        # Remove the elements from `elems_to_remove` from `last_run_elems`
        last_run_elems = match_data_and_remove(last_run_elems, elems_to_remove)

    return last_run_elems