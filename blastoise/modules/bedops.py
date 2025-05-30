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


def bedops_main(data_input, genome_fasta):
    """
    Processes genomic data to create strand-specific sorted BEDOPS files, merges ranges using
    the BEDOPS tool, retrieves genomic sequences, and formats the final output as a DataFrame
    with relevant columns and data.

    Arguments:
        data_input (pd.DataFrame): Input data containing genome coordinates and strand
            information. Expected to have specific columns like 'sstrand', 'sseqid',
            'sstart', 'send', etc.
        genome_fasta (str): Path to the genome FASTA file used for sequence retrieval.

    Returns:
        pd.DataFrame: A processed DataFrame containing genomic ranges with associated
            sequence data and additional calculated attributes such as sequence length.
    """
    # -----------------------------------------------------------------------------
    # 1) Filter and sort data
    # -----------------------------------------------------------------------------
    print("\t\t\t\t- Getting diff. strand data.")
    columns_ids = data_input.columns  # gets the column names
    df_plus = data_input[data_input['sstrand'] == 'plus']  # filters the "+" strand
    df_minus = data_input[data_input['sstrand'] == 'minus']  # filters the "-" strand

    # Sort the data by the start coordinate
    print("\t\t\t\t- Sorting data.")
    df_plus = df_plus.sort_values(by=['sseqid', 'sstart'])  # sorts the "+" strand by the start coordinate
    df_minus = df_minus.sort_values(by=['sseqid', 'sstart'])  # sorts the "-" strand by the start coordinate

    # -----------------------------------------------------------------------------
    # 2) BEDOPS files creation in tmp BASH:
    # -----------------------------------------------------------------------------
    #  row[1] == Chromosome ID, row[10] == Start coordinate, row[11] == End coordinate
    print("\t\t\t\t- Creating BEDOPS files.")
    plus_bedops_bash = get_bedops_bash_file(df_plus)
    minus_bedops_bash = get_bedops_bash_file(df_minus)

    # -----------------------------------------------------------------------------
    # 3) BEDOPS function call with subprocess.
    # -----------------------------------------------------------------------------
    # BEDOPS call to "plus.bed" and "minus.bed" files
    # Using subprocess to call BEDOPS.
    # Using subprocess.check_output() to get the output from the command.
    # shell=True to have the input of .check_output() as a string.
    # universal_newlines=True to have the EoL character as "\n" (Unix-like systems).
    # The output will be a variable of strings.
    print("\t\t\t\t- Calling BEDOPS for plus strand.")
    cmd = f"bedops --merge {plus_bedops_bash}"
    result_plus = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                                 universal_newlines=True, executable="/usr/bin/bash")
    df_plus_bedops = result_plus.stdout

    print("\t\t\t\t- Calling BEDOPS for minus strand.")
    cmd = f"bedops --merge {minus_bedops_bash}"
    result_minus = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                                  universal_newlines=True, executable="/usr/bin/bash")
    df_minus_bedops = result_minus.stdout

    # Now let's transform then into Data Frames
    print("\t\t\t\t- Converting BEDOPS plus output to Data Frames.")
    df_plus_bedops = pd.DataFrame([x.split("\t") for x in df_plus_bedops.split("\n") if x],
                                  columns=['sseqid', 'sstart', 'send'])  # transforms the "+" strand BEDOPS output into a Data Frame
    df_plus_bedops = columns_to_numeric(df_plus_bedops, ['sstart', 'send'])

    print("\t\t\t\t- Converting BEDOPS minus output to Data Frames.")
    df_minus_bedops = pd.DataFrame([x.split("\t") for x in df_minus_bedops.split("\n") if x],
                                   columns=['sseqid', 'sstart', 'send'])  # transforms the "-" strand BEDOPS output into a Data Frame
    df_minus_bedops = columns_to_numeric(df_minus_bedops, ['sstart', 'send'])

    # -----------------------------------------------------------------------------
    # Add column "sstrand" in both data frames.
    # df_plus_bedops with "plus" and df_minus_bedops with "minus"
    df_plus_bedops['sstrand'] = 'plus'
    df_minus_bedops['sstrand'] = 'minus'

    # Join both data frames
    all_data = pd.concat([df_plus_bedops, df_minus_bedops], ignore_index=True)

    # Add "len" column
    all_data['len'] = all_data['send'] - all_data['sstart'] + 1

    # remove temp files
    os.remove(plus_bedops_bash)
    os.remove(minus_bedops_bash)

    return all_data