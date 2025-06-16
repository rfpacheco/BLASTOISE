# Import needed modules
import pandas as pd
import os
import time
from joblib import Parallel, delayed

from typing import Hashable

from modules.bedops import get_bedops_bash_file, bedops_contrast
from modules.genomic_ranges import get_interval_coincidence, get_interval_not_coincidence, merge_intervals, get_merge_stranded
from extra.csv_to_gff import csv_to_gff


# noinspection DuplicatedCode
def match_data(data_input, to_discard):
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
        The DataFrame from which specific columns are used to find matches
        in `data_input`.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing rows from `data_input` that share the same
        values in the specified columns with the corresponding rows in
        `to_discard`.
    """
    matches = data_input.merge(
        to_discard[['sseqid', 'sstart', 'send', 'sstrand']],
        on=['sseqid', 'sstart', 'send', 'sstrand'],
        how='inner'
    )

    return matches


def match_data_and_remove(data_input, to_discard):
    """
    Remove rows from a DataFrame based on matching criteria from another DataFrame.

    This function identifies rows from the given `data_input` DataFrame that match the rows
    in the `to_discard` DataFrame based on specific indexing criteria, and removes the
    matched rows. The matching is performed using the columns 'sseqid', 'sstart', 'send',
    and 'sstrand' as the indices for comparison.

    Parameters:
    -----------
        data_input: pandas.DataFrame
            The input DataFrame from which rows should be evaluated and potentially removed.
        to_discard: pandas.DataFrame
            A DataFrame containing rows to be matched and removed from the `data_input`
            DataFrame based on specified criteria.

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


def match_data_and_set_false(data_input, to_discard):
    """
    Matches data based on specific conditions and updates a column in the input dataset.

    This function takes a primary dataset and a dataset of elements to be discarded, then matches
    elements based on specific coordinate-based keys. It uses these matches to update a boolean
    column in the primary dataset, marking the matched elements as `False` in the 'analyze' column.

    Parameters:
    -----------
    data_input : pandas.DataFrame
        The primary dataset to be analyzed and updated. It must
        contain the columns: 'sseqid', 'sstart', 'send', and 'sstrand'.

    to_discard : pandas.DataFrame
        The dataset containing elements that should be matched
        and marked as not for analysis in the primary dataset.
        Similar to `data_input`, it must also have the columns:
        'sseqid', 'sstart', 'send', and 'sstrand'.

    Returns:
    --------
    None
        The function operates directly on the `data_input` DataFrame, updating
        its 'analyze' column in-place. No new object is returned.
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
    # ──────────────────────────────────────────────────────────────────
    # 0) Build the strand-block structure
    # ──────────────────────────────────────────────────────────────────
    change_points = (all_og_inrange[strand_col] != all_og_inrange[strand_col].shift())
    block_id = change_points.cumsum()
    all_og_inrange = all_og_inrange.assign(_block_id=block_id)

    block_summary = (
        all_og_inrange.groupby("_block_id")
          .agg(block_start=(start_col, "min"),
               block_end=(end_col,   "max"),
               strand=(strand_col,   "first"))
    )

    # ──────────────────────────────────────────────────────────────────
    # Containers
    # ──────────────────────────────────────────────────────────────────
    final_chunks = []
    skip_next_blocks = set()

    # ──────────────────────────────────────────────────────────────────
    # 1) Sliding window on blocks: b0, b1, b2
    # ──────────────────────────────────────────────────────────────────
    ids = block_summary.index.to_list()
    n_blocks = len(ids)

    for idx, b0 in enumerate(ids):
        if b0 in skip_next_blocks:
            continue

        b1 = ids[idx + 1] if idx + 1 < n_blocks else None
        b2 = ids[idx + 2] if idx + 2 < n_blocks else None

        # ──────────────────────────────────────────────────────────────
        # 2) We have at least three consecutive blocks
        # ──────────────────────────────────────────────────────────────
        if b1 is not None and b2 is not None:
            s0 = block_summary.loc[b0, "strand"]
            s1 = block_summary.loc[b1, "strand"]
            s2 = block_summary.loc[b2, "strand"]

            # do outer blocks overlap?
            ## Get all the overlapping elems for block b0 and b2 using `get_interval_coincidence`
            elems_of_b0 = get_interval_coincidence(all_elems_inrange,
                                                   all_og_inrange.loc[all_og_inrange["_block_id"] == b0])
            elems_of_b2 = get_interval_coincidence(all_elems_inrange,
                                                   all_og_inrange.loc[all_og_inrange["_block_id"] == b2])
            ## Check if they have elements overlapping in common
            elems_of_b0_vs_b2 = get_interval_coincidence(elems_of_b0, elems_of_b2)
            outer_overlap = False
            if not elems_of_b0_vs_b2.empty:
                outer_overlap = True

            # do b0 and b1 overlap?
            ## Get all the overlapping elems for block b1 using `get_interval_coincidence`
            elems_of_b1 = get_interval_coincidence(all_elems_inrange,
                                                   all_og_inrange.loc[all_og_inrange["_block_id"] == b1])
            ## Check if they have elements overlapping in common
            elems_of_b0_vs_b1 = get_interval_coincidence(elems_of_b0, elems_of_b1)
            inner_overlap = False
            if not elems_of_b0_vs_b1.empty:
                inner_overlap = True

            # ─ Rule 1: outer blocks share strand & overlap ───────────
            if s0 == s2 and s0 != s1 and outer_overlap:
                outer_rows = pd.concat([elems_of_b0, elems_of_b2])
                og_b1_and_b2 = pd.concat([
                    all_og_inrange.loc[all_og_inrange["_block_id"] == b0],
                    all_og_inrange.loc[all_og_inrange["_block_id"] == b2]
                ])
                all_elems = pd.concat([outer_rows, og_b1_and_b2]).drop_duplicates()
                all_elems.sort_values(start_col)
                collapsed_outer = merge_intervals(all_elems)
                collapsed_outer["sstrand"] = s0
                final_chunks.append(collapsed_outer)
                skip_next_blocks.update({b1, b2})

                continue

            # ─ Rule 2: only b0 ↔ b1 overlap ─────────────────────────
            if inner_overlap:
                elems_of_b1_vs_b0 = get_interval_coincidence(elems_of_b1, elems_of_b0)
                elems_to_remove = pd.concat([elems_of_b1_vs_b0, elems_of_b0_vs_b1]).drop_duplicates()
                elems_to_remove.sort_values(start_col)
                cleaned = match_data_and_remove(all_elems_inrange, elems_to_remove)
                final_chunks.append(cleaned)
                skip_next_blocks.add(b0)
                continue

        # ─ Default: keep block untouched ─────────────────────────────
        untouched = all_og_inrange[all_og_inrange["_block_id"] == b0]
        final_chunks.append(untouched)

    # ──────────────────────────────────────────────────────────────────
    # 3) Build final output
    # ──────────────────────────────────────────────────────────────────
    result = (pd.concat(final_chunks, ignore_index=True, copy=False)
                .drop(columns="_block_id", errors="ignore")   # ← safe drop
                .sort_values(start_col)
                .reset_index(drop=True))
    return result


def _set_overlapping_status_single(chrom: str,
                                   new_df_chr: pd.DataFrame,
                                   og_df_chr: pd.DataFrame,
                                   run_phase: int) -> pd.DataFrame:
    """
    **INTERNAL** helper that contains the ORIGINAL sequential algorithm.

    It is executed by a single worker and receives data restricted to one
    chromosome.  Code below the ‘BEGIN ORIGINAL CODE’ comment is *exactly*
    what used to be inside the old `set_overlapping_status` body – so
    nothing changes logically, we only moved it into a private function.
    """
    # ------------------------------------------------------------------
    # BEGIN ORIGINAL CODE (cut-and-paste of the old implementation)
    # ------------------------------------------------------------------
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
        og_vs_row = get_interval_coincidence(og_df_chr, row)

        # And now, all the elements that overlap with `og_vs_row`. # NOTE: take only True values in `new_df_chr`
        new_elems_in_og_vs_row = get_interval_coincidence(new_df_chr[new_df_chr['analyze'] == True], og_vs_row)

        # Now get all new elements that overlap with 'new_elems_in_og_vs_row'. # NOTE: take only True values in `new_df_chr`
        new_elems_inrange_of_og = get_interval_coincidence(new_df_chr[new_df_chr['analyze'] == True], new_elems_in_og_vs_row)

        # And get all original elements in the whole range of `all_elems_in_range`
        og_inrange = get_interval_coincidence(og_df_chr, new_elems_inrange_of_og)

        # Now, get all elems that overlap all this `og_inrange`
        all_elems_vs_og_inrange = get_interval_coincidence(new_df_chr[new_df_chr['analyze'] == True], og_inrange)

        # If these new elems in `all_elems_vs_og_inrange` overlap some other new element comming from a little far
        # away original element, we need to detect them.
        all_og_inrange = get_interval_coincidence(og_df_chr, all_elems_vs_og_inrange)
        all_elems_inrange = get_interval_coincidence(new_df_chr[new_df_chr['analyze'] == True], all_og_inrange)

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
                elems_of_a =get_interval_coincidence(all_elems_inrange, og_a)
                elems_of_b =get_interval_coincidence(all_elems_inrange, og_b)
                elems_to_remove_a = get_interval_coincidence(elems_of_a, elems_of_b)
                elems_to_remove_b = get_interval_coincidence(elems_of_b, elems_of_a)
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
                    og_vs_row = get_interval_coincidence(og_inrange, row) # Selects original data that overlaps with row
                    all_elems_vs_og = get_interval_coincidence(new_df_chr[new_df_chr['analyze'] == True], og_vs_row)
                    process_overlapping_data(
                        new_df_chr,
                        row,
                        idx,
                        og_vs_row,
                        all_elems_vs_og
                    )
            else:

                # These are not normal cases. For example, is when in the original data is 'minus' -- 'plus' -- 'minus',
                # or 'plus' -- 'minus' -- 'plus', or 'plus' -- 'plus' -- 'minus'.
                ## First, in the `all_og_inrange` detect the strand flip
                change_points = (all_og_inrange['sstrand'] != all_og_inrange['sstrand'].shift())

                # Now check the blocks by number
                block_id = change_points.cumsum()

                if block_id.unique().shape[0] == 2:
                    # if there are only 2 blocks, it means a structure without interlaps like:
                    # 'plus' -- 'plus' -- 'minus' or 'plus' -- 'minus' -- 'minus'
                    ## Take from `all_elems_inrange` the elems with each strand
                    elems_in_plus = all_elems_inrange[all_elems_inrange['sstrand'] == 'plus']
                    elems_in_minus = all_elems_inrange[all_elems_inrange['sstrand'] == 'minus']

                    # Take the elems that overlap each one
                    elems_in_plus_vs_minus = get_interval_coincidence(elems_in_plus, elems_in_minus)
                    elems_in_minus_vs_plus = get_interval_coincidence(elems_in_minus, elems_in_plus)
                    elems_to_remove = pd.concat(
                        [elems_in_plus_vs_minus, elems_in_minus_vs_plus]
                    ).sort_values(['sstart', 'send'])
                    match_data_and_set_false(new_df_chr, elems_to_remove)
                    is_row_removed = match_data(row, elems_to_remove)
                    if not is_row_removed.empty:
                        continue
                    else:
                        elems_in_row = get_interval_coincidence(new_df_chr[new_df_chr['analyze'] == True], row)
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


def set_overlapping_status(new_df: pd.DataFrame,
                           og_df: pd.DataFrame,
                           run_phase: int,
                           n_jobs: int = -1) -> pd.DataFrame:
    """
    Parallel wrapper around the original algorithm.

    Parameters
    ----------
    new_df : DataFrame
        Data to analyse.  Must contain column ``'sseqid'``.
    og_df  : DataFrame
        Reference data.  Must contain column ``'sseqid'``.
    run_phase : int
        Iteration counter
    n_jobs : int, default ``-1``
        Number of processes Joblib should spawn.  (``-1`` ⇒ use all cores,
        ``1`` ⇒ fallback to the original single-process execution.)
    """

    # Fast exit: keep the exact behaviour if the caller explicitly disables
    # parallelism.
    if n_jobs == 1:
        return _set_overlapping_status_single("ALL", new_df, og_df, run_phase)

    # ------------------------------------------------------------------
    # 1. Split the two dataframes by chromosome
    # ------------------------------------------------------------------
    new_groups = {c: df for c, df in new_df.groupby("sseqid", sort=False)}
    og_groups  = {c: df for c, df in og_df.groupby("sseqid", sort=False)}

    # Use only chromosomes that are present in *new* – everything else
    # would produce empty output anyway.
    chromosomes = list(new_groups)

    # ------------------------------------------------------------------
    # 2. Run one worker per chromosome
    # ------------------------------------------------------------------
    results = Parallel(n_jobs=n_jobs, backend="loky")(
        delayed(_set_overlapping_status_single)(
            chrom,
            new_groups[chrom],
            og_groups.get(chrom, og_df.iloc[0:0],),  # empty DF if missing
            run_phase
        )
        for chrom in chromosomes
    )

    # ------------------------------------------------------------------
    # 3. Concatenate the per-chromosome outputs and keep original order
    # ------------------------------------------------------------------
    combined_df = pd.concat(results, ignore_index=True)
    return combined_df

# noinspection DuplicatedCode
def set_strand_direction(data_input: pd.DataFrame,
                         run_phase: int,
                         folder_path: str
 ) -> pd.DataFrame:
    """
     Analyzes and processes genomic sequence data to determine the correct strand orientation for each sequence.
     It handles both new sequences and overlapping sequences, merging and extending where appropriate, ensuring
     data consistency between strands.

     Workflow:
     1. **Prepare the Input Data**:
        - Extracts the original dataset and reconfigures its columns for consistency and comparison.
        - Sorts and deduplicates data by sequence identifiers (`sseqid`) and start positions (`sstart`) to simplify
            operations.

     2. **Split Data into Original and New Sequences**:
        - Separates the input data into `og_data` (already detected sequences) and `new_data` (new detected elements).

     3. **Identify Overlaps and Non-Overlapping Sequences**:
        - Determines where the `new_data` sequences overlap or do not overlap with the `og_data` using BEDOPS tools.
        - Non-overlapping sequences are treated as *new elements*, while overlapping sequences move to further
            processing.

     4. **Process Overlapping Elements**:
        - For overlapping elements, evaluates the strand orientation:
            a. If an overlap occurs with the same strand as the existing sequences, they are **merged and extended**
                if there's room for extension.
            b. If the strands differ, conflicting sequences are removed (sequences from inverse strands cannot coexist).
            c. For multiple overlapping matches, the new elements are merged with the existing sequences only with the
                ones they share the same orientation.

     5. **Output Combined Results**:
        - Combines the processed overlapping sequences with the new elements from `new_data` that had no conflicts.
        - Produces a clean, consistent DataFrame with updated strand orientations and merged sequences.

     6. **Clean up Temporary Files**:
        - Filters and removes temporary files generated during the BEDOPS processing.

     Key Features:
     - **Strand Orientation Evaluation**:
       Ensures that sequences are merged correctly if overlaps are found and validates strand (plus/minus direction)
        alignment.
     - **Handling of Multiple Overlaps**:
       Manages situations where overlaps occur across multiple existing sequences, ensuring data integrity and avoiding
        strand conflicts.
     - **Performance Tracking**:
       Tracks execution time of strand computations for performance insights.
     - **Integration with BEDOPS**:
       Relies on BEDOPS to perform overlap, contrast, and merging operations efficiently.

    Parameters:
    -----------
        data_input (pd.DataFrame): Input DataFrame containing sequence data with columns such as
                                   'og_sseqid', 'og_sstart', 'og_send', 'og_sstrand', 'sseqid',
                                   'sstart', 'send', and 'sstrand'.
        run_phase (int): Indicates the phase of the analysis. Used to determine whether to run the
        folder_path (str): Path to the folder where temporary files are stored.

    Returns:
    --------
        pd.DataFrame: A DataFrame with processed data that includes new elements and overlapping
                      elements with orientation properly set.
    """

    # Save the original coordinates before the extension in `og_data`
    og_data = data_input[['og_sseqid', 'og_sstart', 'og_send', 'og_sstrand']].copy()
    og_data.columns = ['sseqid', 'sstart', 'send', 'sstrand'] # Change column names
    og_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    og_data.drop_duplicates(inplace=True)
    print(f"\t\t\t- Already existing elements row length: {og_data.shape[0]}")

    # Save all the new elements discovered in `new_data`
    new_data = data_input[['sseqid', 'sstart', 'send', 'sstrand']].copy()
    new_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    new_data.drop_duplicates(inplace=True)
    print(f"\t\t\t- Potential newly discovered elements row length: {new_data.shape[0]}")

    # Get the `new_data` that doesn't overlap with `og_data`
    new_no_overlaps_og = get_interval_not_coincidence(new_data, og_data)

    # `new_no_overlaps_og` are considered new elements. Let's save them correctly. From 'new_data', extract the elements
    # that have the same 'sseqid, 'sstart', 'send' and 'sstrand' as in `new_no_overlaps_og`
    new_elems = pd.merge(new_data, new_no_overlaps_og,
                         on=['sseqid', 'sstart', 'send', 'sstrand'],
                         how='inner')
    print(f"\t\t\t- New elements: {new_elems.shape[0]}")

    # And the rest elements, which for sure, overlap
    overlapping_elems = pd.merge(new_data, new_no_overlaps_og,
                                 on=['sseqid', 'sstart', 'send', 'sstrand'],
                                 how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
    print(f"\t\t\t- Overlapping elements: {overlapping_elems.shape[0]}")

    # Now we need to check if the overlapping is in the same strand as the original data
    original_elems_plus = pd.DataFrame() # Initialization
    original_elems_minus = pd.DataFrame() # Initialization
    if not overlapping_elems.empty:
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
        multiprocessing_jobs = -1
        overlapping_elems = set_overlapping_status(overlapping_elems, og_data, run_phase, n_jobs=multiprocessing_jobs)
        overlapping_elems = get_merge_stranded(overlapping_elems) # Merge the data
        toc = time.perf_counter()
        print(f"\t\t\t\t- Execution time: {toc - tic:0.2f} seconds")

        # Now dive the data in 'minus' and plus' strand if there are rows with 'plus' or 'minus'
        original_elems_plus = overlapping_elems[overlapping_elems['sstrand'] == 'plus'].copy()
        original_elems_minus = overlapping_elems[overlapping_elems['sstrand'] == 'minus'].copy()

    if not new_elems.empty:
        # Now dive the data in 'minus' and plus' strand
        new_elems_plus = new_elems[new_elems['sstrand'] == 'plus'].copy()
        new_elems_minus = new_elems[new_elems['sstrand'] == 'minus'].copy()

        # First, for each new element before being merged. Let's remove elements in the contrary strands that overlap with each other
        new_elems_plus_overlap_with_new_elems_minus = get_interval_coincidence(
            new_elems_plus, new_elems_minus
        )
        new_elems_minus_overlap_with_new_elems_plus = get_interval_coincidence(
            new_elems_minus, new_elems_plus
        )

        # And remove them
        if not new_elems_plus_overlap_with_new_elems_minus.empty:
            new_elems_plus = match_data_and_remove(new_elems_plus, new_elems_plus_overlap_with_new_elems_minus)

        if not new_elems_minus_overlap_with_new_elems_plus.empty:
            new_elems_minus = match_data_and_remove(new_elems_minus, new_elems_minus_overlap_with_new_elems_plus)

        # Not let's merge the results
        new_elems_plus = merge_intervals(new_elems_plus)
        new_elems_minus = merge_intervals(new_elems_minus)

        # And make again the tmp files
        new_elems_plus_bedops = get_bedops_bash_file(new_elems_plus)
        new_elems_minus_bedops = get_bedops_bash_file(new_elems_minus)

        # Remove the elements that overlap with the original sequences
        if not original_elems_minus.empty:  # Removing new elements in minus strand that overlaps with original strands
            new_elems_plus_overlaps_original_minus = get_interval_coincidence(
                new_elems_plus, original_elems_minus
            )
            if not new_elems_plus_overlaps_original_minus.empty: # MAIN TAKE
                new_elems_plus = match_data_and_remove(new_elems_plus, new_elems_plus_overlaps_original_minus)

        if not original_elems_plus.empty:  # Removing new elements in plus strand that overlaps with original strands
            new_elems_minus_overlaps_original_plus = bedops_contrast(new_elems_minus_bedops, original_elems_plus,
                                                                     'coincidence')
            if not new_elems_minus_overlaps_original_plus.empty: # MAIN TAKE
                new_elems_minus = match_data_and_remove(new_elems_minus, new_elems_minus_overlaps_original_plus)

        # Join again in new_elems and remove tmp files
        os.remove(new_elems_plus_bedops)  # Removing tmp file with no more use
        os.remove(new_elems_minus_bedops)  # Removing tmp file with no more use
        new_elems = pd.concat([new_elems_plus, new_elems_minus])
        new_elems.sort_values(by=['sseqid', 'sstart'], inplace=True)

    # Combine new elements and processed leftover elements
    result = pd.concat([new_elems, overlapping_elems])

    return result
