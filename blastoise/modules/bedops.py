# Modules needed
import pandas as pd
import numpy as np
import subprocess
import os

from modules.files_manager import columns_to_numeric

def get_data_sequence(data, strand, genome_fasta):
    """
    This function gets the sequence of the data from the fasta file. It will keep the Chromosome ID, start coordinate, end coordinate and strand.

    :param data: A pandas data frame with the data read of the BED files.
    :type data: pandas.core.frame.DataFrame

    :param strand: The strand of the sequence. It can be "plus" or "minus".
    :type strand: str

    :param genome_fasta: Path to the whole genome sequence in FASTA format.
    :type genome_fasta: string
    """
    sequences = []
    for _, row in data.iterrows():
        sseqid = row["sseqid"]
        start = row["sstart"]
        end = row["send"]
        cmd = [
            "blastdbcmd",
            "-db", genome_fasta,
            "-entry", sseqid,
            "-range", f"{start}-{end}",
            "-strand", strand,
            "-outfmt", "%s"
        ]

        sequence = subprocess.check_output(cmd, universal_newlines=True).replace('\n', '')

        sequences.append({
            "sseqid": sseqid,
            "sstart": start,
            "send": end,
            "sstrand": strand,
            "sseq": sequence
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


def bedops_main(data_input, genome_fasta, writing_path_input):
    """
    Processes genomic data by filtering, sorting, and merging BEDOPS files, and obtaining sequences.
    The resulting data is restructured and returned as a formatted DataFrame.

    Args:
        data_input (pandas.DataFrame): The input data frame containing genomic information.
        genome_fasta (str): Path to the genome FASTA file.
        writing_path_input (str): Directory where BEDOPS files will be written.

    Returns:
        pandas.DataFrame: A DataFrame with processed and formatted genomic data.
    """
    # -----------------------------------------------------------------------------
    # 1) Filter and sort data
    # -----------------------------------------------------------------------------
    columns_ids = data_input.columns  # gets the columns names
    df_plus = data_input[data_input["sstrand"] == "plus"]  # filters the "+" strand
    df_minus = data_input[data_input["sstrand"] == "minus"]  # filters the "-" strand

    # Sort the data by the start coordinate
    df_plus = df_plus.sort_values(by=["sseqid", "sstart"])  # sorts the "+" strand by the start coordinate
    df_minus = df_minus.sort_values(by=["sseqid", "sstart"])  # sorts the "-" strand by the start coordinate

    # -----------------------------------------------------------------------------
    # 2) BEDOPS files creation:
    # -----------------------------------------------------------------------------
    #  row[1] == Chromosome ID, row[10] == Start coordinate, row[11] == End coordinate
    plus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_plus.bed")
    minus_path = os.path.join(writing_path_input, os.path.basename(writing_path_input) + "_minus.bed")
    
    df_plus[["sseqid", "sstart", "send"]].to_csv(plus_path, sep="\t", header=False, index=False)  # creates a BED file for the "+" strand
    # In the minus strand it's important to change the order of the coordinates, because "bedops" reads them like "plus" strand.
    # If not, it will not merge them.
    df_minus[["sseqid", "send", "sstart"]].to_csv(minus_path, sep="\t", header=False, index=False)  # creates a BED file for the "-" strand

    # -----------------------------------------------------------------------------
    # 3) BEDOPS function call with subprocess.
    # -----------------------------------------------------------------------------
    # BEDOPS call to "plus.bed" and "minus.bed" files
    # Using subprocess to call BEDOPS.
    # Using subprocess.check_output() to get the output from the command.
    # shell=True to have the input of .check_output() as a string.
    # universal_newlines=True to have the EoL character as "\n" (Unix-like systems).
    # The output will be a variable of strings.
    df_plus_bedops = subprocess.check_output(f"bedops --merge {plus_path}", shell=True, universal_newlines=True)  # merges the "+" strand BED file
    df_minus_bedops = subprocess.check_output(f"bedops --merge {minus_path}", shell=True, universal_newlines=True)  # merges the "-" strand BED file

     # Now let's transform then into Data Frames
    df_plus_bedops = pd.DataFrame([x.split("\t") for x in df_plus_bedops.split("\n") if x],
                                    columns=["sseqid", "sstart", "send"])  # transforms the "+" strand BEDOPS output into a Data Frame
    df_minus_bedops = pd.DataFrame([x.split("\t") for x in df_minus_bedops.split("\n") if x],
                                    columns=["sseqid", "sstart", "send"])  # transforms the "-" strand BEDOPS output into a Data Frame
    # -----------------------------------------------------------------------------
    # 4) Call `blastdbcmd` to get the sequences with the function get_data_sequence()
    # -----------------------------------------------------------------------------
    if df_plus_bedops.empty:  # In case the original data is empty, the code needs to keep going
        df_plus_bedops_wseq = pd.DataFrame()  # creates an empty Data Frame
    else:  # If the original data is not empty, tit uses get_data_sequence
        df_plus_bedops_wseq = get_data_sequence(df_plus_bedops, "plus", genome_fasta)

    # The same for the minus strand:
    if df_minus_bedops.empty:
        df_minus_bedops_wseq = pd.DataFrame()
    else:   
        df_minus_bedops_wseq = get_data_sequence(df_minus_bedops, "minus", genome_fasta)

 
    # Let's reorderthe `df_minus_bedps_wseq` data frame:
    df_minus_bedops_wseq[["sstart", "send"]] = df_minus_bedops_wseq[["send", "sstart"]].copy()  # swap only values
    # -----------------------------------------------------------------------------
    # 5) Processing data
    # -----------------------------------------------------------------------------
    # Join both data frames
    all_data = pd.concat([df_plus_bedops_wseq, df_minus_bedops_wseq], ignore_index=True)  # joins both Data Frames

    # Adding sequence length to the DataFrame:
    if not all_data.empty:
        new_column = np.array([len(x) for x in all_data.loc[:, "sseq"]], dtype=np.int64)  # creates a list with the length of each sequence
        all_data.insert(1, "length", new_column)  # inserts the new column with the sequence length. Column index are shifted.
        # -----------------------------------------------------------------------------
        # 6) Correctly modeling the output Data Frame to 15 columns and output as CSV file.
        # -----------------------------------------------------------------------------
        new_data = pd.DataFrame(index=range(all_data.shape[0]), columns=columns_ids)  # creates a new Data Frame with 15 columns. The rows depend on the .shape[0]

        new_data.loc[:,["sseqid", "length", "sstart", "send", "sstrand", "sseq"]] = all_data.loc[:,["sseqid", "length", "sstart", "send", "sstrand", "sseq"]].copy()
        new_data = columns_to_numeric(new_data, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
        
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


