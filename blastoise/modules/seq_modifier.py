import numpy as np
import pandas as pd
import subprocess
import re
import math
import os

from modules.bedops import bedops_main

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def specific_sequence_1000nt(data_input, main_folder_path, genome_fasta, extend_number):
    """
    Expands sequences in a given DataFrame to a specific length (1000 nucleotides).

    Parameters:
    data_input: DataFrame
        The DataFrame containing sequence information.
    main_folder_path: str
        Path to the main folder.
    genome_fasta: str
        Path to the genome FASTA file.
    extend_number: int
        The length to which sequences need to be extended.

    Returns:
    DataFrame
        The modified DataFrame with sequences expanded to the specified length.
    """
    # -----------------------------------------------------------------------------
    for index2, (index, element) in enumerate(data_input.iterrows()):
        if "plus" in element["sstrand"]:
            lower_coor = int(element["sstart"])  # We get the start of the sequence
            upper_coor = int(element["send"])  # We get the end of the sequence
        else:  # If it's the "-" strand
            lower_coor = int(element["send"])
            upper_coor = int(element["sstart"])
        
        subject_length = upper_coor - lower_coor + 1
        if subject_length < extend_number:  # If the sequence is less than 1000nt, we'll expand it.
            leftover_length = extend_number - subject_length  # We get the difference between 1000 and the length of the sequence
            leftover_length_halved = leftover_length / 2  # We divide it by 2 to get the half of the difference
            # leftover_length_halved = math.ceil(leftover_length_halved)  # We round up the number
            lower_coor = lower_coor - leftover_length_halved  # We subtract the half of the difference to the start
            upper_coor = upper_coor + leftover_length_halved  # We add the half of the difference to the end

            # Check if the coordinates are not out of the genome
            if lower_coor <= 0:
                lower_coor = 1
                upper_coor += leftover_length_halved  # Since `new_start` is 1, we add the half of the difference to the end
            
            # We get the sequence from the whole genome
            cmd = f"blastdbcmd -db {genome_fasta} -entry {element['sseqid']} -range {int(lower_coor)}-{int(upper_coor)} -strand {element['sstrand']} -outfmt %s"
            seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

            # Now with the data let's modify the data frame
            data_input.loc[index, "pident"] = np.nan
            data_input.loc[index, "length"] = np.nan
            data_input.loc[index, "qstart"] = np.nan
            data_input.loc[index, "qend"] = np.nan

            if "plus" in element["sstrand"]:
                data_input.loc[index, "sstart"] = int(lower_coor)
                data_input.loc[index, "send"] = int(upper_coor)
            else:
                data_input.loc[index, "sstart"] = int(upper_coor)
                data_input.loc[index, "send"] = int(lower_coor)

            data_input.loc[index, "evalue"] = np.nan
            data_input.loc[index, "bitscore"] = np.nan
            data_input.loc[index, "qlen"] = np.nan
            data_input.loc[index, "slen"] = len(seq)
            data_input.loc[index, "sseq"] = seq
        else:  # If the sequence is already 1000nt we just pass
            pass
        
    return data_input


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# 3) Corrector se secuencias, obtendra las secuencias originales.


def specific_sequence_corrected(data_input, nucleotides1000_df, first_data_input, main_folder_path, chromosome_ID, genome_fasta, run_phase):
    """
    The main use of this function is to get the real "coordinates" of the sequence.


    1. First we read the 1000nt CSV file. And we extract from row[0] (i.e., Query row) every *different* query (not repeated.
    2. Then we open the 1000nt_BLASTER file from the 1000nt against each other. To understand it:

       - In the 1000nt file there are "x" number of rows. Each row will be given an index like "Seq_z_LinJ.01_plus", were "z" is the index.
       - So the first row (row[0]) in the 1000nt file will be given "Seq_1.." as a name. And every row with the "Seq_1.." subject in the 1000nt_BLASTER file will the result from the first row in the 1000nt file (row[0]).
       - There will be a maximum of "x" "Seq_x.." in the 1000nt_BLASTER file.

    3. In the 1000nt_BLASTER we search every result from a specific "Seq_z..". And we get the maximum alignment and it's coordinates.
    4. We'll get the "z" number and search it's corresponding sequence in the 1000nt file. And with that and the coordinates, we use ``blastcmd`` to get the correct coordinates from the 1000nt sequence. We need to be careful if the sequence is "plus" or "minus".

    :param data_input: Path to the CSV file we'll use to filter data. It's the result from a BLASTn made between the expanded 1000nt sequences.
    :type data_input: pandas DataFrame

    :param nucleotides1000_directory: *return* result from the function :func:`~specific_sequence_1000nt`
    :type nucleotides1000_directory: string

    :param main_folder_path: Path where we'll save the CSV data file to.
    :type main_folder_path: string

    :param chromosome_ID: Identification of the cromosome, e.g., "LinJ.07"
    :type chromosome_ID: string

    :return: CSV File with the corrected coordinates of our sequences.
    :rtype: CSV file
    """

    # First from the BLASTn (one againts each other), we get the IDs of the sequences.
    folder_corrected_sequences_modified_path = os.path.join(main_folder_path, "modified")
    folder_corrected_sequences_not_modified_path = os.path.join(main_folder_path, "not_modified")
    os.makedirs(folder_corrected_sequences_modified_path, exist_ok=True)
    os.makedirs(folder_corrected_sequences_not_modified_path, exist_ok=True)


    data_input_selected = data_input[data_input["qseqid"] != data_input["sseqid"]].copy()  # We filter the data to get only the sequences that are different from each other.

    selected_data_sseqid = data_input_selected["sseqid"].unique()  # We get the unique qseqid
    data_input_rejected = data_input[data_input["qseqid"] == data_input["sseqid"]].copy()  # We filter the data to get only the sequences that are the same.   
    data_input_rejected = data_input_rejected[~data_input_rejected["sseqid"].isin(selected_data_sseqid)]  # Now get the rows in "data_input_rejected" that are not in "data_input_selected"

    selected_grouped = data_input_selected.groupby("sseqid")  # We group the data by the subject ID
    coor_dic = {} # We create an empty data frame to store the data
    for name, group in selected_grouped:
        group["difference"] = abs(group["send"].astype(int) - group["sstart"].astype(int))  # We get the difference between the start and end coordinates. Important the "abs" since the coordinates are differente depending if it's the "+" or "-" strand
        max_index = group["difference"].idxmax()  # We get the index of the maximum difference

        # Now, we need to differentiate between the "+" and "-" strand coordinates. But not in the data yet.
        if int(group.loc[max_index, "sstart"]) < int(group.loc[max_index, "send"]):  # like in "plus" strand
            lower_coor = group.loc[max_index, "sstart"]
            upper_coor = group.loc[max_index, "send"]
        else:  # If it's the "-" strand
            lower_coor = group.loc[max_index, "send"]
            upper_coor = group.loc[max_index, "sstart"]

        index = int(re.search(r'\d+', name).group())
        # Now we add all the data to the dictionary
        coor_dic[name] = [index, lower_coor, upper_coor]
    
    # Loop through "coor_dic" to get the sequences in the 1000nt data frame
    tmp_list = []
    for _, value in coor_dic.items():
        index = value[0]
        lower_coor = int(value[1])
        upper_coor = int(value[2])
        len_real_coor = abs(int(first_data_input.iloc[index,:]["sstart"]) - int(first_data_input.iloc[index,:]["send"])) + 1
        strand = first_data_input.iloc[index,:]["sstrand"]
        qseqid = first_data_input.iloc[index,:]["qseqid"]

        # We get the sequence from the whole genome, But first, we need to adjust the coordinates.
        if strand == "plus":
            new_lower_coor = int(first_data_input.iloc[index,:]["sstart"]) + lower_coor - 1  # -1 because if lower_coor == 1, then we need to get the first nucleotide
            upper_coor_diff = len_real_coor - upper_coor
            new_upper_coor = int(first_data_input.iloc[index,:]["send"]) - upper_coor_diff # We subtract the difference between the length of the sequence and the upper coordinate
        else:  # If it's the "-" strand
            new_lower_coor = int(first_data_input.iloc[index,:]["send"]) + lower_coor - 1  # -1 because if lower_coor == 1, then we need to get the first nucleotide
            upper_coor_diff = len_real_coor - upper_coor
            new_upper_coor = int(first_data_input.iloc[index,:]["sstart"]) - upper_coor_diff
        
        cmd = f"blastdbcmd -db {genome_fasta} -entry {chromosome_ID} -range {new_lower_coor}-{new_upper_coor} -strand {strand} -outfmt %s"
        seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()
        # Now we add the data to the data frame
        # We reset the sstart and send position for minus strand
        if strand == "minus":
            new_lower_coor, new_upper_coor = new_upper_coor, new_lower_coor
        # And now we add the data to the data frame
        tmp_list.append([qseqid,  # qseqid
                         chromosome_ID,  # sseqid
                         "",  # pident
                         "",  # length
                         "",  # qstart 
                         "",  # qend
                         new_lower_coor,  # sstart
                         new_upper_coor,  # send
                         "",  # evalue
                         "",  # bitscore
                         "",  # qlen,
                         len(seq),  # slen
                         strand,  # sstrand
                         seq])  # sseq
    
    data_new_coordinates = pd.DataFrame(tmp_list, columns=["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand", "sseq"])
    data_new_coordinates_path = os.path.join(folder_corrected_sequences_modified_path, f"{run_phase}_modified.csv")
    data_new_coordinates.to_csv(data_new_coordinates_path, index=False, header=True, sep=",")

    # Now for the rejected data. This one only hits with themselves
    tmp_list2 = []
    for _, row in data_input_rejected.iterrows():  # We iterate through the rejected data
        index = int(re.search(r'\d+', row["sseqid"]).group())  # We get the index of the sequence
        recatched_row = first_data_input.iloc[index]  # We get the original data from the first data input
        tmp_list2.append(recatched_row)
    
    data_rejected = pd.DataFrame(tmp_list2)
    data_rejected_path = os.path.join(folder_corrected_sequences_not_modified_path, f"{run_phase}_not_modified.csv")
    data_rejected.to_csv(data_rejected_path, index=False, header=True, sep=",")

    
    all_data = pd.concat([data_new_coordinates, data_rejected], ignore_index=True)  # We concatenate the data

    return all_data
    
    