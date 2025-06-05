# Import needed modules
import pandas as pd
import os
import time

from modules.bedops import get_bedops_bash_file, bedops_contrast, bedops_main


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


def filter_redundant_seqs(df_to_filter, df_to_contrast, df_to_discard):
    """
    From a set of sequences selects those that share the same strand as `df_to_contrast`, this will be later merged
     with `df_to_contrast` coordinates to generate one merged sequence. The elements that do not share the same
     coordinates and strand orientation as the 'merged one' will be selected and discarded later.

    Parameters
    ----------
    df_to_filter : pd.DataFrame
        Contains sequences that need to be filtered.
    df_to_contrast : pd.DataFrame
        Contains a sequence used to get the correct strand for `df_to_filter` as well as more coordinates to be merged
        with `df_to_filter`.
    df_to_discard :
        Empty data frame that will be filled with sequences to be discarded later.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        same_elems_merged: pd.DataFrame
            Contains ONE sequences in the same strand as `df_to_contrast`. It comes from the merged sequences from
            the data in `df_to_filter` that has the same strand as `df_to_contrast`.
        df_to_discard: pd.DataFrame
            Contains sequences that shall be discarded later on. Sequences that don't share the same strand as
            `df_to_contrast` or don't have the same coordinates as in `same_elems_merged`.
    """


    # Get the elements from `df_to_filter` that has the same strand as `df_to_contrast`
    same_elems = df_to_filter[df_to_filter['sstrand'] == df_to_contrast['sstrand'].iloc[0]].copy()
    # Same but in opposite
    not_same_elems = df_to_filter[df_to_filter['sstrand'] != df_to_contrast['sstrand'].iloc[0]].copy()

    if not not_same_elems.empty: # If it has lines
        df_to_discard = pd.concat([df_to_discard, not_same_elems]) # Fil the discard dataset

    same_elems_bedops = get_bedops_bash_file(same_elems) # tmp bedops files

    # Merge only the elems from `same_elems` using its bedops file
    same_elems_merged = bedops_contrast(
        same_elems_bedops,
        "",  # To force `same_elems_bedops` to merge by itself, the second argument is an empty string
        'merge'
    )

    # Add the 'strand' column to the merged pd.DataFrame
    same_elems_merged['sstrand'] = df_to_contrast['sstrand'].iloc[0]

    # Get elems in `same_elems` that don't have the same coordinates as in `same_elems_merged`
    merged_start = same_elems_merged.iloc[0]['sstart'] # Get start coord
    merged_end = same_elems_merged.iloc[0]['send'] # Get end coord

    # Filter rows that don't match the start and end coordinates from the merged data
    filtered = same_elems[
        (same_elems['sstart'] != merged_start) | (same_elems['send'] != merged_end)
    ]

    if not filtered.empty: # If it has lines
        df_to_discard = pd.concat([df_to_discard, filtered])

    # Remove tmp files
    os.remove(same_elems_bedops)

    return same_elems_merged, df_to_discard


def edit_og_data_and_get_merged_seq(original_overlapping_with_row_bedops, 
                                    new_overlapping_data_that_overlaps_with_selected_original, 
                                    row_df, to_discard, new_overlapping_data, index):
    """
    Edits and merges data based on overlap information between datasets, updating the primary data
    and filtering redundant sequences.

    This function processes overlapping genomic or sequence data. It works by merging overlapping
    sequences from an input dataset with corresponding data from a contrast dataset. The function also
    handles redundant sequence filtering, updates the 'analyze' column of the input dataset to mark
    matches for exclusion from analysis, and removes temporary files used in processing.

    Parameters:
    ----------
    original_overlapping_with_row_bedops : str
        Path to the file containing overlapping data rows that belong to the contrast dataset.
    new_overlapping_data_that_overlaps_with_selected_original : pandas.DataFrame
        Subset of the input dataset containing rows that overlap with the contrast dataset.
    row_df : pandas.Series
        A single row dataframe or series representing the sequence to be compared and updated.
    to_discard : pandas.DataFrame
        Containing sequences that should be excluded or filtered out in the process.
    new_overlapping_data : pandas.DataFrame
        The primary dataset containing genomic or sequence data, which will be updated during the
        function execution.
    index : int
        Index of the row in `new_overlapping_data` that corresponds to the row being updated.

    Returns:
    ----------
    pandas.DataFrame
        The updated `to_discard` DataFrame after processing and filtering redundant sequences.
    """

    # Get the merged coordinates data from `new_overlapping_data_that_overlaps_with_selected_original` and the data set to discard
    merged_coordinates, to_discard = filter_redundant_seqs(
        new_overlapping_data_that_overlaps_with_selected_original,
        row_df,
        to_discard
    )

    merged_coordinates_bedops = get_bedops_bash_file(merged_coordinates) # tmp bedops file

    # Merge our sequence from `merged_coordinates_bedops` with the contrast one in `original_overlapping_with_row_bedops`
    extended_sequence = bedops_contrast(
        merged_coordinates_bedops,
        original_overlapping_with_row_bedops,
        'merge'
    )

    # Replace in the original index in `new_overlapping_data` the new merged sequence
    new_overlapping_data.loc[index, ['sstart', 'send']] = extended_sequence.iloc[0][['sstart', 'send']]

    # Data input will be an 'in/out' parameter. Discard is an 'in' parameter
    match_data_and_set_false(new_overlapping_data, to_discard)

    # Remove tmp files
    os.remove(merged_coordinates_bedops)


def set_overlapping_status(new_overlapping_data, original_overlapping_data):
    """
    Analyzes and updates overlapping genomic regions based on strand compatibility. For each 
    genomic region in the input data, determines if it overlaps with regions in the contrast data. 
    If overlaps exist and strands match, the regions are extended and merged. Regions with 
    mismatched strands are marked for removal.

    Parameters
    ----------
    new_overlapping_data : pd.DataFrame
        Set of sequences that overlap with one or more sequences in `original_overlapping_data`.

    original_overlapping_data : str  
        Path to a BEDOPS file containing the original sequences that `new_overlapping_data` overlaps with.

    Returns
    -------
    pd.DataFrame
        Processed dataframe with:
        - Updated start/end positions for merged overlapping regions
        - Removed rows for non-overlapping or strand-mismatched regions
        - Original rows converted to integer coordinates
    """

    # Set the 'analyze' column to True
    new_overlapping_data['analyze'] = True  # Only True values will be analyzed, the rest will be skipped.
    new_overlapping_data_bedops = get_bedops_bash_file(new_overlapping_data)

    for index, row in new_overlapping_data.iterrows():
        if not new_overlapping_data.at[index, 'analyze']:
            # IMPORTANT: with `at` we access the original DataFrame, not the copy. Because in `.iterrows()`, we access
            # the copy and not the original one
            continue # If the column is false, skip it

        # Set row as a dataframe of one line
        row_df = pd.DataFrame(row).T
        row_df_bedops = get_bedops_bash_file(row_df) # tmp bedops file

        # Get the elements from `original_overlapping_data` that overlap with `row_df`
        original_overlaps_with_row_df = bedops_contrast(
            original_overlapping_data,
            row_df_bedops,
            'coincidence'
        )

        # TODO: create a "fast" method to filter from 'new_overlapping_data_bedops' only the 'True' elements from 'new_overlapping_data'
        # Get elements from the original `new_overlapping_data` that overlap with `original_overlaps_with_row_df`
        original_overlaps_with_row_df_bedops = get_bedops_bash_file(original_overlaps_with_row_df) # tmp bedops file
        new_overlapping_data_that_overlaps_with_selected_original = bedops_contrast(
            new_overlapping_data_bedops,
            original_overlaps_with_row_df_bedops,
            'coincidence'
        )

        # Create the "discard" dataframe
        discard_df = pd.DataFrame(columns=['sseqid', 'sstart', 'send', 'sstrand'])

        # If they are in the same strand, should extend.
        if original_overlaps_with_row_df.shape[0] == 1: # means row_df overlapped with only 1 elem
            if row_df['sstrand'].iloc[0] == original_overlaps_with_row_df['sstrand'].iloc[0]: # Both in same strand
                edit_og_data_and_get_merged_seq(
                    original_overlaps_with_row_df_bedops,
                    new_overlapping_data_that_overlaps_with_selected_original,
                    row_df,
                    discard_df,
                    new_overlapping_data,
                    index
                )
            else:  # They are not in the same strand
                # Take from `new_overlapping_data_that_overlaps_with_selected_original` the elements that are in the same strand as `row_df`
                elem_in_same_strand = new_overlapping_data_that_overlaps_with_selected_original[
                    new_overlapping_data_that_overlaps_with_selected_original['sstrand'] == row_df['sstrand'].iloc[0]
                ]
                match_data_and_set_false(new_overlapping_data, elem_in_same_strand)
        else: # means `row_df` overlapped with > 1 elem
            # ---------------------------------------------------------------
            dict_counter = {}
            for _, elem in original_overlaps_with_row_df.iterrows():
                if elem['sstrand'] in dict_counter:
                    dict_counter[elem['sstrand']] += 1
                else:
                    dict_counter[elem['sstrand']] = 1
            dict_len = len(dict_counter)
            # ---------------------------------------------------------------
            if dict_len == 1: # The overlapping elems are in the same strand
                if row_df['sstrand'].iloc[0] == list(dict_counter.keys())[0]:
                    edit_og_data_and_get_merged_seq(
                        original_overlaps_with_row_df_bedops,
                        new_overlapping_data_that_overlaps_with_selected_original,
                        row_df,
                        discard_df,
                        new_overlapping_data,
                        index
                    )
                else:  # They are not in the same strand
                    # In this case, `row_df` doesn't match in 'strand' with any `original_overlaps_with_row_df`
                    # So the every element in `new_overlapping_data_that_overlaps_with_selected_original` with the same 'strand as 'row_df'
                    # shall be removed
                    elem_in_same_strand = new_overlapping_data_that_overlaps_with_selected_original[
                        new_overlapping_data_that_overlaps_with_selected_original['sstrand'] == row_df['sstrand'].iloc[0]
                    ]
                    match_data_and_set_false(new_overlapping_data, elem_in_same_strand)
            elif dict_len > 1:  # The overlapping elems are in different strands
                # We will add the data to the strand it matches.
                # The one with inverse match, will be ignored
                for _, elem in original_overlaps_with_row_df.iterrows():
                    if elem['sstrand'] == row_df['sstrand'].iloc[0]: # Only when they are in the same strand
                        elem = pd.DataFrame(elem).T
                        elem_bedops = get_bedops_bash_file(elem)
                        edit_og_data_and_get_merged_seq(elem_bedops,
                                                        new_overlapping_data_that_overlaps_with_selected_original,
                                                        row_df,
                                                        discard_df,
                                                        new_overlapping_data,
                                                        index
                        )
                        os.remove(elem_bedops) # Remove tmp file
                    else: # if the element doesn't match the strand
                        pass
            else: # if dict_len == 0
                pass

        os.remove(row_df_bedops)
        os.remove(original_overlaps_with_row_df_bedops)

    new_overlapping_data.dropna(inplace=True)
    new_overlapping_data['sstart'] = new_overlapping_data['sstart'].astype(int)
    new_overlapping_data['send'] = new_overlapping_data['send'].astype(int)

    final_data = new_overlapping_data[new_overlapping_data['analyze'] == True].copy()

    # Remove tmp files
    os.remove(new_overlapping_data_bedops)

    return final_data


def set_strand_direction(data_input):
    """
     Analyzes and processes genomic sequence data to determine the correct strand orientation for each sequence.
     It handles both new sequences and overlapping sequences, merging and extending where appropriate, ensuring
     data consistency between strands.

     Workflow:
     1. **Prepare the Input Data**:
        - Extracts the original dataset and reconfigures its columns for consistency and comparison.
        - Sorts and deduplicates data by sequence identifiers (`sseqid`) and start positions (`sstart`) to simplify operations.

     2. **Split Data into Original and New Sequences**:
        - Separates the input data into `original_contrast_data` (already detected sequences) and `new_data` (new detected elements).

     3. **Identify Overlaps and Non-Overlapping Sequences**:
        - Determines where the `new_data` sequences overlap or do not overlap with the `original_contrast_data` using BEDOPS tools.
        - Non-overlapping sequences are treated as *new elements*, while overlapping sequences move to further processing.

     4. **Process Overlapping Elements**:
        - For overlapping elements, evaluates the strand orientation:
            a. If an overlap occurs with the same strand as the existing sequences, they are **merged and extended** if there's room for extension.
            b. If the strands differ, conflicting sequences are removed (sequences from inverse strands cannot coexist).
            c. For multiple overlapping matches, the new elements are merged with the existing sequences only with the ones they share the same orientation.

     5. **Output Combined Results**:
        - Combines the processed overlapping sequences with the new elements from `new_data` that had no conflicts.
        - Produces a clean, consistent DataFrame with updated strand orientations and merged sequences.

     6. **Clean up Temporary Files**:
        - Filters and removes temporary files generated during the BEDOPS processing.

     Key Features:
     - **Strand Orientation Evaluation**:
       Ensures that sequences are merged correctly if overlaps are found and validates strand (plus/minus direction) alignment.
     - **Handling of Multiple Overlaps**:
       Manages situations where overlaps occur across multiple existing sequences, ensuring data integrity and avoiding strand conflicts.
     - **Performance Tracking**:
       Tracks execution time of strand computations for performance insights.
     - **Integration with BEDOPS**:
       Relies on BEDOPS to perform overlap, contrast, and merging operations efficiently.

    Arguments:
        data_input (pd.DataFrame): Input DataFrame containing sequence data with columns such as
                                   'og_sseqid', 'og_sstart', 'og_send', 'og_sstrand', 'sseqid',
                                   'sstart', 'send', and 'sstrand'.

    Returns:
        pd.DataFrame: A DataFrame with processed data that includes new elements and overlapping
                      elements with orientation properly set.
    """

    # Save the original coordinates before the extension in `original_contrast_data`
    original_contrast_data = data_input[['og_sseqid', 'og_sstart', 'og_send', 'og_sstrand']].copy()
    original_contrast_data.columns = ['sseqid', 'sstart', 'send', 'sstrand'] # Change column names
    original_contrast_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    original_contrast_data.drop_duplicates(inplace=True)
    print(f"\t\t\t- Already existing elements row length: {original_contrast_data.shape[0]}")

    # Save all the new elements discovered in `new_data`
    new_data = data_input[['sseqid', 'sstart', 'send', 'sstrand']].copy()
    new_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    new_data.drop_duplicates(inplace=True)
    print(f"\t\t\t- Potential newly discovered elements row length: {new_data.shape[0]}")

    # From both pd.DataFrames, create a tmp .bed file using BEDOPS tools.
    original_bedops = get_bedops_bash_file(original_contrast_data)
    new_bedops = get_bedops_bash_file(new_data)

    # Get the newly discovered elements that do not match with the already existing sequences
    new_data_no_overlaps_original = bedops_contrast(new_bedops, original_bedops, 'opposite') # TODO: there's a special case, new data. But that will overlap other data after `set_overlapping_status` takes place

    # `new_data_no_overlaps_original` are considered new elements. Let's save them correctly.
    # From 'new_data', extract the elements that have the same 'sseqid, 'sstart', 'send' and 'sstrand' as in
    # `new_data_no_overlaps_original`
    new_elems = pd.merge(new_data, new_data_no_overlaps_original,
                         on=['sseqid', 'sstart', 'send', 'sstrand'],
                         how='inner')
    print(f"\t\t\t- New elements: {new_elems.shape[0]}")

    # And the rest elements, which for sure, overlap
    overlapping_elems = pd.merge(new_data, new_data_no_overlaps_original,
                                 on=['sseqid', 'sstart', 'send', 'sstrand'],
                                 how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
    print(f"\t\t\t- Overlapping elements: {overlapping_elems.shape[0]}")

    # Now we need to check if the overlapping is in the same strand as the original data
    original_elems_plus_bedops = '' # Initialization
    original_elems_minus_bedops = '' # Initialization
    if not overlapping_elems.empty:
        tic = time.perf_counter()
        print(f"\t\t\t- Checking strand orientation in overlapping elements:")
        overlapping_elems = set_overlapping_status(overlapping_elems, original_bedops)
        overlapping_elems = bedops_main(overlapping_elems) # Merge the data
        toc = time.perf_counter()
        print(f"\t\t\t\t- Execution time: {toc - tic:0.2f} seconds")

        # Now dive the data in 'minus' and plus' strand if there are rows with 'plus' or 'minus'
        original_elems_plus = overlapping_elems[overlapping_elems['sstrand'] == 'plus'].copy()
        original_elems_minus = overlapping_elems[overlapping_elems['sstrand'] == 'minus'].copy()

        # Transform to bedops tmp files
        original_elems_plus_bedops = get_bedops_bash_file(original_elems_plus)
        original_elems_minus_bedops = get_bedops_bash_file(original_elems_minus)

    # Let's merge the `new_elems`
    if not new_elems.empty:
        new_elems = bedops_main(new_elems)

        # Now dive the data in 'minus' and plus' strand
        new_elems_plus = new_elems[new_elems['sstrand'] == 'plus'].copy()
        new_elems_minus = new_elems[new_elems['sstrand'] == 'minus'].copy()

        # Transform to bedops tmp files
        new_elems_plus_bedops = get_bedops_bash_file(new_elems_plus)
        new_elems_minus_bedops = get_bedops_bash_file(new_elems_minus)

        # Remove the elements that overlap with the original sequences
        if original_elems_minus_bedops != '':  # Removing new elements in minus strand that overlaps with original strands
            new_elems_plus_overlaps_original_minus = bedops_contrast(new_elems_plus_bedops, original_elems_minus_bedops,
                                                                     'coincidence')
            if not new_elems_plus_overlaps_original_minus.empty: # MAIN TAKE
                new_elems_plus = match_data_and_remove(new_elems_plus, new_elems_plus_overlaps_original_minus)
            os.remove(original_elems_minus_bedops)  # Removing tmp file with no more use

        if original_elems_plus_bedops != '':  # Removing new elements in plus strand that overlaps with original strands
            new_elems_minus_overlaps_original_plus = bedops_contrast(new_elems_minus_bedops, original_elems_plus_bedops,
                                                                     'coincidence')
            if not new_elems_minus_overlaps_original_plus.empty: # MAIN TAKE
                new_elems_minus = match_data_and_remove(new_elems_minus, new_elems_minus_overlaps_original_plus)
            os.remove(original_elems_plus_bedops)  # Removing tmp file with no more use

        # Join again in new_elems and remove tmp files
        os.remove(new_elems_plus_bedops)  # Removing tmp file with no more use
        os.remove(new_elems_minus_bedops)  # Removing tmp file with no more use
        new_elems = pd.concat([new_elems_plus, new_elems_minus])
        new_elems.sort_values(by=['sseqid', 'sstart'], inplace=True)

    # Combine new elements and processed leftover elements
    result = pd.concat([new_elems, overlapping_elems])

    # Remove temp files
    os.remove(original_bedops)
    os.remove(new_bedops)

    return result
