import pandas as pd
import os
import time

from modules.bedops import get_bedops_bash_file, bedops_contrast

def plus_and_minus_dataframe_splitter(data_input, column):
    plus_data = data_input[data_input[column] == 'plus'].copy()
    minus_data = data_input[data_input[column] == 'minus'].copy()

    return plus_data, minus_data

def set_overlapping_status(data_input, contrast_data_bedops):
    for index, row in data_input.iterrows():
        # Set row as a dataframe of one line
        row_df = pd.DataFrame(row).T
        sequence_for_bedops = get_bedops_bash_file(row_df)

        # We know that row_df will match with `contrast_data_bedops`. But we want to know in which strand.
        # Get the elements from `contrast_data_bedops` that overlap with `row_df`
        contrast_overlaps_with_row_df = bedops_contrast(contrast_data_bedops, sequence_for_bedops, 'coincidence')

        # If they are in the same strand, should extend.
        # Compare the pandas.Series
        if contrast_overlaps_with_row_df.shape[0] == 1:
            if row_df['sstrand'].iloc[0] == contrast_overlaps_with_row_df['sstrand'].iloc[0]:
                contrast_overlaps_with_row_df_bedops = get_bedops_bash_file(contrast_overlaps_with_row_df)
                extended_sequence = bedops_contrast(sequence_for_bedops,
                                                    contrast_overlaps_with_row_df_bedops,
                                                    'merge')
                data_input.loc[index, ['sstart', 'send']] = extended_sequence.iloc[0][['sstart', 'send']]
                os.remove(contrast_overlaps_with_row_df_bedops)
            else:  # They are not in the same strand
                data_input.loc[index, :] = pd.NA
        else: # if `contrast_overlaps_with_row_df.shape[0] > 1`
            # ---------------------------------------------------------------
            dict_counter = {} # TODO: remove in future, it's only a checker
            for _, elem in contrast_overlaps_with_row_df.iterrows():
                if elem['sstrand'] in dict_counter:
                    dict_counter[elem['sstrand']] += 1
                else:
                    dict_counter[elem['sstrand']] = 1
            dict_len = len(dict_counter)
            # ---------------------------------------------------------------
            if dict_len == 1:
                if row_df['sstrand'].iloc[0] == list(dict_counter.keys())[0]:
                    contrast_overlaps_with_row_df_bedops = get_bedops_bash_file(contrast_overlaps_with_row_df)
                    extended_sequence = bedops_contrast(sequence_for_bedops,
                                                        contrast_overlaps_with_row_df_bedops,
                                                        'merge')
                    data_input.loc[index, ['sstart', 'send']] = extended_sequence.iloc[0][['sstart', 'send']]
                    os.remove(contrast_overlaps_with_row_df_bedops)
                else:  # They are not in the same strand
                    data_input.loc[index, :] = pd.NA
            elif dict_len > 1:  # In this case, the match will be with different strands
                # TODO: review with JMR
                # We will add the data to the strand it matches.
                # The one with inverse match, will be ignored
                for _, elem in contrast_overlaps_with_row_df.iterrows():
                    if elem['sstrand'] == row_df['sstrand'].iloc[0]:
                        # Transform elem from pandas.core.series.Series to pandas.core.frame.DataFrame
                        elem = pd.DataFrame(elem).T
                        contrast_overlaps_with_row_df_bedops = get_bedops_bash_file(elem) # Only the element that matches strand
                        extended_sequence = bedops_contrast(sequence_for_bedops,
                                                            contrast_overlaps_with_row_df_bedops,
                                                            'merge')
                        data_input.loc[index, ['sstart', 'send']] = extended_sequence.iloc[0][['sstart', 'send']]
                        os.remove(contrast_overlaps_with_row_df_bedops)
                    else: # if the element doesn't match the strand
                        pass
            else: # if dict_len == 0
                pass

        os.remove(sequence_for_bedops)

    data_input.dropna(inplace=True)
    data_input['sstart'] = data_input['sstart'].astype(int)
    data_input['send'] = data_input['send'].astype(int)

    # Remove temp files

    return data_input


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

     6. **Cleanup Temporary Files**:
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
    new_data_no_overlaps_contrast = bedops_contrast(new_bedops, original_bedops, 'opposite')

    # `new_data_no_overlaps_contrast` are considered new elements. Let's save them correctly.
    # From 'new_data', extract the elements that have the same 'sseqid, 'sstart', 'send' and 'sstrand' as in
    # `new_data_no_overlaps_contrast`
    new_elements = pd.merge(new_data, new_data_no_overlaps_contrast,
                            on=['sseqid', 'sstart', 'send', 'sstrand'],
                            how='inner')
    print(f"\t\t\t- New elements: {new_elements.shape[0]}")

    # And the rest elements, which for sure, overlap
    leftover_elements = pd.merge(new_data, new_data_no_overlaps_contrast,
                                 on=['sseqid', 'sstart', 'send', 'sstrand'],
                                 how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)
    print(f"\t\t\t- Overlapping elements: {leftover_elements.shape[0]}")

    # Now we need to check if the overlapping is in the same strand as the original data
    if not leftover_elements.empty:
        tic = time.perf_counter()
        print(f"\t\t\t- Checking strand orientation in overlapping elements:")
        leftover_elements = set_overlapping_status(leftover_elements, original_bedops)
        toc = time.perf_counter()
        print(f"\t\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    # Combine new elements and processed leftover elements
    result = pd.concat([new_elements, leftover_elements])

    # Remove temp files
    os.remove(original_bedops)
    os.remove(new_bedops)


    return result
