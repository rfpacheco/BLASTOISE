# Modules needed
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os

from modules.files_manager import columns_to_numeric, end_always_greater_than_start

def get_bedops_bash_file(data):
    """
    Generate a string representing a BEDOPS-compatible bash file.

    This function constructs a string formatted as a BEDOPS-compatible input,
    which contains genomic regions information extracted from the input data.
    The data is processed and formatted into a string representation that
    can be used in subsequent BEDOPS operations.

    Args:
        data (pd.DataFrame): A pandas DataFrame where each row contains genomic
            region information. The required columns are 'sseqid', 'sstart',
            and 'send'.

    Returns:
        str: A formatted string representing the BEDOPS-compatible input.
    """

# Create a temporary file to store the BEDOPS-compatible data
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        # Write each row of the DataFrame into the temporary file in BEDOPS format
        for _, row in data.iterrows():
            temp_file.write(f"{row['sseqid']}\t{row['sstart']}\t{row['send']}\t{row['sstrand']}\n")

        # Store the name of the temporary file for use in the BEDOPS command
        bedops_file = temp_file.name

    # TODO: remember to os.remove every object created with this function
    return bedops_file

def bedops_contrast(base_df_path, contrast_df_path, bedops_mode):
    """
    Executes a BEDOPS command to compare genomic regions between two BED files and returns a DataFrame.

        Parameters:
        base_df_path (str): Path to the base BED file.
        contrast_df_path (str): Path to the contrast BED file.
        bedops_mode (str): Mode of operation for the BEDOPS command.
                           Supported modes are 'coincidence' (checks which elements in the base are in contrast),
                           'opposite' (checks which elements in the base are not in contrast),
                           and 'merge' (merges overlapping intervals).

        Returns:
        pd.DataFrame: DataFrame containing the results of the BEDOPS operation with columns 'sseqid', 'sstart', and 'send'.
    """
    bedops_mode_map = {'coincidence': '--element-of 1',
                       'opposite': '--not-element-of 1',
                       'merge': '--merge'}
    cmd_mode = bedops_mode_map.get(bedops_mode)

    # Shows the elements in `base_df` that overlap with `contrast_df`
    cmd_coincidence = f"bedops {cmd_mode} {base_df_path} {contrast_df_path}"
    check_coincidence = subprocess.run(cmd_coincidence, shell=True, capture_output=True, text=True,
                                       universal_newlines=True,
                                       executable='/usr/bin/bash').stdout.strip()
    if bedops_mode == 'merge':
        columns_needed = ['sseqid', 'sstart', 'send']
    else:
        columns_needed = ['sseqid', 'sstart', 'send', 'sstrand']

    check_coincidence = pd.DataFrame([x.split("\t") for x in check_coincidence.split("\n") if x],
                                     columns=columns_needed)
    check_coincidence[['sstart', 'send']] = check_coincidence[['sstart', 'send']].astype(int)

    return check_coincidence


def bedops_coincidence(main_data, data_for_contrast, strand, genome_fasta):
    """
    Processes and identifies overlapping and non-overlapping genomic data between
    two datasets using BEDOPS tools, sorts and merges the results, and retrieves
    their sequences based on the genomic strand and a given genome FASTA file.

    Parameters:
    main_data (pd.DataFrame):
        The primary dataset containing genomic data to be analyzed.
    data_for_contrast (pd.DataFrame): 
        The secondary dataset used for comparison against the main dataset.
    strand (str):
        The genomic strand orientation for which sequences need to be extracted
        ('+' or '-').
    genome_fasta (str):
        Path to the genome FASTA file for sequence extraction.

    Returns:
    tuple
        A tuple containing:
        - coincidence_data (pd.DataFrame):
            Data overlapping in both `main_data` and `data_for_contrast`
            after merging.
        - new_data (pd.DataFrame):
            Data unique to `main_data` (not found in `data_for_contrast`). 
        - only_in_contrast_data (pd.DataFrame):
            Data unique to `data_for_contrast` (not found in `main_data`).
    """
    main_data = main_data.sort_values(by=['sseqid', 'sstart'])  # Sort the data frame by the start coordinate
    data_for_contrast = data_for_contrast.sort_values(by=['sseqid', 'sstart'])  # Sort the data frame by the start coordinate
    main_data_len = main_data.shape[0] # Get length of the main data
    data_contrast_len = data_for_contrast.shape[0] # Get length of the contrast data

    main_data_bedops = get_bedops_bash_file(main_data) # get a tmp bash file for the main data
    contrast_data_bedops = get_bedops_bash_file(data_for_contrast) # get a tmp bash file for contrast data

    # Check elements in the "main data" that overlap in the "contrast" data
    main_exists_in_contrast_data = bedops_contrast(main_data_bedops, contrast_data_bedops, 'coincidence')
    print("")
    print("\t\t\t- Coincidence data:")
    if main_data_len > 0:
        print(f"\t\t\t\t- New data in Previous data: {main_exists_in_contrast_data.shape[0]}/{main_data_len} - {main_exists_in_contrast_data.shape[0]/main_data_len*100:.2f}%")
    else: # In this case `main_data_len == 0`
        print(f"\t\t\t\t- New data in Previous data: {main_exists_in_contrast_data.shape[0]}/{main_data_len}")

    contrast_exists_in_main_data = bedops_contrast(contrast_data_bedops, main_data_bedops, 'coincidence')
    if data_contrast_len > 0:
        print(f"\t\t\t\t- Previous data in New data: {contrast_exists_in_main_data.shape[0]}/{data_contrast_len} - {contrast_exists_in_main_data.shape[0]/data_contrast_len*100:.2f}%")
    else: # Then `data_contrast_len == 0`
        print(f"\t\t\t\t- Previous data in New data: {contrast_exists_in_main_data.shape[0]}/{data_contrast_len}")
    # -----------------------------------------------------------------------------
    # There would be elements that are in both datas. The next step is to merge them.
    main_exists_in_contrast_data_bedops = get_bedops_bash_file(main_exists_in_contrast_data)
    contrast_exists_in_main_data_bedops = get_bedops_bash_file(contrast_exists_in_main_data)

    merged_data = bedops_contrast(main_exists_in_contrast_data_bedops, contrast_exists_in_main_data_bedops, 'merge')
    print(f"\t\t\t\t- Merged data: {merged_data.shape[0]}")

    coincidence_data = merged_data.copy()
    coincidence_data['sstrand'] = strand
    # ----------------- ------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Check elements from `main_data` that DO NOT overlap with `data_for_contrast`-
    # Because these elements are not in the contrast data, they will be novice elements.
    print("")
    print("\t\t\t- NOT coincidence data:")
    # noinspection DuplicatedCode
    main_not_exists_in_contrast_data = bedops_contrast(main_data_bedops, contrast_data_bedops, 'opposite')
    if main_data_len > 0:
        print(f"\t\t\t\t- New data NOT in Previous data: {main_not_exists_in_contrast_data.shape[0]}/{main_data_len} - {main_not_exists_in_contrast_data.shape[0]/main_data_len*100:.2f}%")
    else:
        print(f"\t\t\t\t- New data NOT in Previous data: {main_not_exists_in_contrast_data.shape[0]}/{main_data_len}")

    if not main_not_exists_in_contrast_data.empty:  # If the data frame has data
        new_data = main_not_exists_in_contrast_data.copy()
    else:  # If the data frame is empty
        new_data = pd.DataFrame()
    # -----------------------------------------------------------------------------
    # Now check the elements in Old that are not in Last
    # noinspection DuplicatedCode
    contrast_not_exists_in_main_data = bedops_contrast(contrast_data_bedops, main_data_bedops, 'opposite')
    if data_contrast_len > 0:
        print(f"\t\t\t\t- Previous data NOT in New data: {contrast_not_exists_in_main_data.shape[0]}/{data_contrast_len} - {contrast_not_exists_in_main_data.shape[0]/data_contrast_len*100:.2f}%")
    else:
        print(f"\t\t\t\t- Previous data NOT in New data: {contrast_not_exists_in_main_data.shape[0]}/{data_contrast_len}")

    if not contrast_not_exists_in_main_data.empty:  # If the data frame has lines
        only_in_contrast_data = contrast_not_exists_in_main_data.copy()
    else:  # If the data frame is empty
        only_in_contrast_data = pd.DataFrame()

    # -----------------------------------------------------------------------------
    # Remove temp files from "get_bedops_bash_file()" function
    os.remove(main_data_bedops)
    os.remove(contrast_data_bedops)
    os.remove(main_exists_in_contrast_data_bedops)
    os.remove(contrast_exists_in_main_data_bedops)

    return coincidence_data, new_data, only_in_contrast_data