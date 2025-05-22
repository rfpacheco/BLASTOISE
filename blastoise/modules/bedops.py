# Modules needed
import pandas as pd
import numpy as np
import subprocess
import os

from modules.files_manager import columns_to_numeric

def get_data_sequence(data, strand, genome_fasta):
    """
    Fetches nucleotide sequences from a specified genome database using BLAST commands.

    This function retrieves specified sequence ranges from a genome database in a fasta format.
    The input strand direction is specified, and BLAST commands are executed to obtain the
    sequences based on the provided data.

    Args:
        data (DataFrame): A pandas DataFrame containing the sequence details. It must
        contain the columns 'sseqid', 'sstart', and 'send'.

        strand (str): A string indicating the strand direction for sequence retrieval, typically
        'plus' or 'minus'.

        genome_fasta (str): The file path to the genome database in fasta format to query against.

    Returns:
        DataFrame: A pandas DataFrame containing the retrieved sequences. The DataFrame includes
        columns for 'sseqid', 'sstart', 'send', 'sstrand', and 'sseq'.

    Raises:
        subprocess.CalledProcessError: If the `blastdbcmd` tool encounters an error during execution.
    """
    sequences = []
    for _, row in data.iterrows():
        sseqid = row['sseqid']
        start = row['sstart']
        end = row['send']
        cmd = f"blastdbcmd -db {genome_fasta} -entry {sseqid} -range {start}-{end} -strand {strand} -outfmt %s"

        sequence = subprocess.check_output(cmd, universal_newlines=True).replace('\n', '')

        sequences.append({
            'sseqid': sseqid,
            'sstart': start,
            'send': end,
            'sstrand': strand,
            'sseq': sequence
        })

    sequences_df = pd.DataFrame(sequences)

    return sequences_df

def bedops_contrast(base_df_path, contrast_df_path, bedops_mode):
    """
    Executes a BEDOPS command to compare genomic regions between two BED files and returns a DataFrame.

        Parameters:
        base_df_path (str): Path to the base BED file.
        contrast_df_path (str): Path to the contrast BED file.
        bedops_mode (str): Mode of operation for the BEDOPS command.
                           Supported modes are 'coincidence' (checks which elements in base are in contrast),
                           'opposite' (checks which elements in base are not in contrast),
                           and 'merge' (merges overlapping intervals).

        Returns:
        pd.DataFrame: DataFrame containing the results of the BEDOPS operation with columns 'sseqid', 'sstart', and 'send'.
    """
    bedops_mode_map = {'coincidence': '--element-of 1',
                       'opposite': '--not-element-of 1',
                       'merge': '--merge'}
    cmd_mode = bedops_mode_map.get(bedops_mode)

    # Check which elements in 'base_df' are inside 'contrast_df'
    cmd_coincidence = f"bedops {cmd_mode} {base_df_path} {contrast_df_path}"
    check_coincidence = subprocess.check_output(cmd_coincidence, shell=True, universal_newlines=True)
    check_coincidence = pd.DataFrame([x.split("\t") for x in check_coincidence.split("\n") if x],
                             columns=["sseqid", "sstart", "send"])
    check_coincidence = columns_to_numeric(check_coincidence, ["sstart", "send"])
    return check_coincidence


def bedops_main(data_input, genome_fasta):
    """
    Processes genomic data to generate and manipulate BEDOPS files, retrieve sequences,
    and structure the resulting data into a standardized DataFrame format.

    Parameters:
        data_input (pd.DataFrame): The input DataFrame containing genomic data to process. Must include specific
            columns such as 'sstrand', 'sseqid', 'sstart', 'send', etc.
        genome_fasta (str): Path to the genome FASTA file to extract sequences using BEDOPS outputs.

    Returns:
        pd.DataFrame: A new DataFrame with processed genomic data comprising 15 columns. Returns
            an empty DataFrame if the input data or intermediate computation yields no result.

    Raises:
        None
    """
    # -----------------------------------------------------------------------------
    # 1) Filter and sort data
    # -----------------------------------------------------------------------------
    columns_ids = data_input.columns  # gets the column names
    df_plus = data_input[data_input['sstrand'] == 'plus']  # filters the "+" strand
    df_minus = data_input[data_input['sstrand'] == 'minus']  # filters the "-" strand

    # Sort the data by the start coordinate
    df_plus = df_plus.sort_values(by=['sseqid', 'sstart'])  # sorts the "+" strand by the start coordinate
    df_minus = df_minus.sort_values(by=['sseqid', 'sstart'])  # sorts the "-" strand by the start coordinate

    # -----------------------------------------------------------------------------
    # 2) BEDOPS files creation in tmp BASH:
    # -----------------------------------------------------------------------------
    #  row[1] == Chromosome ID, row[10] == Start coordinate, row[11] == End coordinate
    plus_bedops_bash = "<(echo -e '" + '\n'.join(
        [f"{row['sseqid']}\t{row['sstart']}\t{row['send']}" for _, row in df_plus.iterrows()]) + "')"
    minus_bedops_bash = "<(echo -e '" + '\n'.join(
        [f"{row['sseqid']}\t{row['send']}\t{row['sstart']}" for _, row in df_minus.iterrows()]) + "')"
    
    # -----------------------------------------------------------------------------
    # 3) BEDOPS function call with subprocess.
    # -----------------------------------------------------------------------------
    # BEDOPS call to "plus.bed" and "minus.bed" files
    # Using subprocess to call BEDOPS.
    # Using subprocess.check_output() to get the output from the command.
    # shell=True to have the input of .check_output() as a string.
    # universal_newlines=True to have the EoL character as "\n" (Unix-like systems).
    # The output will be a variable of strings.
    cmd = f"bedops --merge {plus_bedops_bash}"
    result_plus = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                                 universal_newlines=True, executable="/usr/bin/bash")
    df_plus_bedops = result_plus.stdout

    cmd = f"bedops --merge {minus_bedops_bash}"
    result_minus = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                                  universal_newlines=True, executable="/usr/bin/bash")
    df_minus_bedops = result_minus.stdout

     # Now let's transform then into Data Frames
    df_plus_bedops = pd.DataFrame([x.split("\t") for x in df_plus_bedops.split("\n") if x],
                                  columns=['sseqid', 'sstart', 'send'])  # transforms the "+" strand BEDOPS output into a Data Frame
    df_minus_bedops = pd.DataFrame([x.split("\t") for x in df_minus_bedops.split("\n") if x],
                                   columns=['sseqid', 'sstart', 'send'])  # transforms the "-" strand BEDOPS output into a Data Frame
    # -----------------------------------------------------------------------------
    # 4) Call `blastdbcmd` to get the sequences with the function get_data_sequence()
    # -----------------------------------------------------------------------------
    if df_plus_bedops.empty:  # In case the original data is empty, the code needs to keep going
        df_plus_bedops_seq = pd.DataFrame()  # creates an empty Data Frame
    else:  # If the original data is not empty, tit uses get_data_sequence
        df_plus_bedops_seq = get_data_sequence(df_plus_bedops, 'plus', genome_fasta)

    # The same for the minus strand:
    if df_minus_bedops.empty:
        df_minus_bedops_seq = pd.DataFrame()
    else:   
        df_minus_bedops_seq = get_data_sequence(df_minus_bedops, 'minus', genome_fasta)

 
    # Let's reorder the `df_minus_bedops_seq` data frame:
    df_minus_bedops_seq[['sstart', 'send']] = df_minus_bedops_seq[['send', 'sstart']].copy()  # swap only values

    # -----------------------------------------------------------------------------
    # 5) Processing data
    # -----------------------------------------------------------------------------
    # Join both data frames
    all_data = pd.concat([df_plus_bedops_seq, df_minus_bedops_seq], ignore_index=True)  # joins both Data Frames

    # Adding sequence length to the DataFrame:
    if not all_data.empty:
        new_column = np.array([len(x) for x in all_data.loc[:, 'sseq']], dtype=np.int64)  # Creates a list with the length of each sequence
        all_data.insert(1, 'length', new_column)  # Inserts the new column with the sequence length. Column index is shifted.
        # -----------------------------------------------------------------------------
        # 6) Correctly modeling the output Data Frame to 15 columns and output as CSV file.
        # -----------------------------------------------------------------------------
        new_data = pd.DataFrame(index=range(all_data.shape[0]), columns=columns_ids)  # Creates a new Data Frame with 15 columns. The rows depend on the .shape[0]

        columns_needed = ['sseqid', 'sstart', 'send', 'sstrand', 'sseq', 'length']
        new_data.loc[:, columns_needed] = all_data.loc[:, columns_needed].copy()
        new_data = columns_to_numeric(new_data, ['length', 'sstart', 'send'])
        
        return new_data  # returns the new Data Frame
    else:
        return pd.DataFrame()  # returns an empty Data Frame


def bedops_coincidence(last_df, old_df, folder_path, strand, genome_fasta):
    """
    Will tell the elements from old_df that are in last_df.
    last_df is the last data frame that we have.
    old_df is the first input data frame.
    """
    last_df = last_df.sort_values(by=["sseqid", "sstart"])  # Sort the data frame by the start coordinate
    old_df = old_df.sort_values(by=["sseqid", "sstart"])  # Sort the data frame by the start coordinate
    last_length = last_df.shape[0]
    old_length = old_df.shape[0]

    last_df_path = os.path.join(f"{folder_path}_1.1_last.bed")
    old_df_path = os.path.join(f"{folder_path}_1.2_old.bed")

    last_df[["sseqid", "sstart","send"]].to_csv(last_df_path, sep="\t", header=False, index=False)
    old_df[["sseqid", "sstart","send"]].to_csv(old_df_path, sep="\t", header=False, index=False)

    last_in_old = bedops_contrast(last_df_path, old_df_path, 'coincidence')
    print("")
    print("\t\t\t- Coincidence data:")
    print(f"\t\t\t\t- New data in Previous data: {last_in_old.shape[0]}/{last_length} - {last_in_old.shape[0]/last_length*100:.2f}%")

    old_in_last = bedops_contrast(old_df_path, last_df_path, 'coincidence')
    print(f"\t\t\t\t- Previous data in New data: {old_in_last.shape[0]}/{old_length} - {old_in_last.shape[0]/old_length*100:.2f}%")
    # -----------------------------------------------------------------------------
    # Let's merge


    old_in_last_path = os.path.join(f"{folder_path}_2.1_Old_in_Last.bed")
    last_in_old_path = os.path.join(f"{folder_path}_2.2_Last_in_Old.bed")

    old_in_last[["sseqid", "sstart", "send"]].to_csv(old_in_last_path, sep="\t", header=False, index=False)
    last_in_old[["sseqid", "sstart", "send"]].to_csv(last_in_old_path, sep="\t", header=False, index=False)

    merged_last_old = bedops_contrast(old_in_last_path, last_in_old_path, 'merge')
    print(f"\t\t\t\t- Merged data: {merged_last_old.shape[0]}")

    # Now recapture the elements with the genome
    coincidence_data = get_data_sequence(merged_last_old, strand, genome_fasta)
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Now let's check the elements that are not in df2 (the first input). They would be the new elements.
    print("")
    print("\t\t\t- NOT coincidence data:")
    last_notin_old = bedops_contrast(last_df_path, old_df_path, 'opposite')
    print(f"\t\t\t\t- New data NOT in Previous data: {last_notin_old.shape[0]}/{last_length} - {last_notin_old.shape[0]/last_length*100:.2f}%")

    if not last_notin_old.empty:  # If the data frame is not empty
        new_data = get_data_sequence(last_notin_old, strand, genome_fasta)
    else:  # If the data frame is empty
        new_data = pd.DataFrame()
    # -----------------------------------------------------------------------------
    # Now check the elements in Old that are not in Last
    old_notin_last = bedops_contrast(old_df_path, last_df_path, 'opposite')
    print(f"\t\t\t\t- Previous data NOT in New data: {old_notin_last.shape[0]}/{old_length} - {old_notin_last.shape[0]/old_length*100:.2f}%")

    if not old_notin_last.empty:  # If the data frame is not empty
        old_data_exclusive = get_data_sequence(old_notin_last, strand, genome_fasta)
    else:  # If the data frame is empty
        old_data_exclusive = pd.DataFrame()

    # -----------------------------------------------------------------------------
    return coincidence_data, new_data, old_data_exclusive


