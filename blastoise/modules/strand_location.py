import pandas as pd
import os

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
        if row_df['sstrand'].name == contrast_overlaps_with_row_df['sstrand'].name:
            contrast_overlaps_with_row_df_bedops = get_bedops_bash_file(contrast_overlaps_with_row_df)
            extended_sequence = bedops_contrast(sequence_for_bedops,
                                                contrast_overlaps_with_row_df_bedops,
                                                'merge')
            data_input.loc[index, ['sstart', 'send']] = extended_sequence.iloc[0][['sstart', 'send']]
            os.remove(contrast_overlaps_with_row_df_bedops)
        else:  # They are not in the same strand
            data_input.loc[index, :] = pd.NA


    data_input.dropna(inplace=True)

    # Remove temp files
    # noinspection PyUnboundLocalVariable
    os.remove(sequence_for_bedops)
    
    return data_input


def set_strand_direction(data_input):
    # First, let's take the original data with coordinates
    original_contrast_data = data_input[['og_sseqid', 'og_sstart', 'og_send', 'og_sstrand',]].copy()
    original_contrast_data.columns = ['sseqid', 'sstart', 'send', 'sstrand'] # Change column names
    original_contrast_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    original_contrast_data.drop_duplicates(inplace=True)
    


    # Let's split it in "plus" and "minus"
    # original_contrast_plus, original_contrast_minus = plus_and_minus_dataframe_splitter(original_contrast_data, 'sstrand')

    # Do the same with the new data
    new_data = data_input[['sseqid', 'sstart', 'send', 'sstrand']].copy()
    new_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
    new_data.drop_duplicates(inplace=True)
    # new_plus, new_minus = plus_and_minus_dataframe_splitter(new_data, 'sstrand')

    # Let's check if the new sequences overlap with the original
    original_bedops = get_bedops_bash_file(original_contrast_data)
    new_bedops = get_bedops_bash_file(new_data)

    # And where there is not a coincidence
    new_data_no_overlaps_contrast = bedops_contrast(new_bedops, original_bedops, 'opposite')
    

    # `new_data_no_overlaps_contrast` are considered new elements. Let's save them correctly.
    # From 'new_data', extract the elements that have the same 'sseqid, 'sstart', 'send' and 'sstrand' as in
    # `new_data_no_overlaps_contrast`
    new_elements = pd.merge(new_data, new_data_no_overlaps_contrast,
                            on=['sseqid', 'sstart', 'send', 'sstrand'],
                            how='inner')

    # And the rest elements, which for sure, overlap
    leftover_elements = pd.merge(new_data, new_data_no_overlaps_contrast,
                                 on=['sseqid', 'sstart', 'send', 'sstrand'],
                                 how='left', indicator=True).query('_merge == "left_only"').drop('_merge', axis=1)

    # Now we need to check if the overlapping is in the same strand as the original data
    if not leftover_elements.empty:
        leftover_elements = set_overlapping_status(leftover_elements, original_bedops)

    # Combine new elements and processed leftover elements
    result = pd.concat([new_elements, leftover_elements])

    # Remove temp files
    os.remove(original_bedops)
    os.remove(new_bedops)


    return result
