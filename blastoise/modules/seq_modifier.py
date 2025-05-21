import numpy as np
import subprocess

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def specific_sequence_1000nt(data_input, genome_fasta, extend_number, limit_len):
    """
    Adjusts subject sequences in a data frame to meet a specified length requirement.

    Iterates through a data frame containing sequence alignment details, evaluates the length
    of each subject sequence, and adjusts the coordinates if the sequence is shorter than
    the specified limit. If the sequence is shorter, it extends the sequence by adjusting
    start and end positions equally on both sides while avoiding coordinates outside genome
    bounds. Also updates relevant columns in the data frame with the new sequence details.

    Parameters:
    data_input (pd.DataFrame): A data frame containing sequence alignment details. Expected to have the following 
        columns: 'sstrand', 'sstart', 'send', 'sseqid', and alignment metrics such as
        'pident', 'length', 'qstart', etc.
    genome_fasta (str): The path to the genome FASTA file used for extracting sequence data.
    extend_number (int): The minimum length required for the sequence. If a sequence is shorter than this
        value, it's extended to this length.
    limit_len (int): A specific length limit that if met, the process of extension is skipped.

    Returns:
    pd.DataFrame: A modified data frame where sequences having lengths below the specified minimum
        are extended. Updated metrics and extracted sequence are reflected in the modified
        data frame.

    Raises:
    subprocess.CalledProcessError: If an error occurs when running the system command to retrieve the genome sequence.
    KeyError: If any required column is missing in the provided data frame.

    Notes:
    The function assumes that the input 'data_input' follows a specific structure with
    the necessary columns to process coordinate adjustments and sequence extraction.
    Coordinates are adjusted considering both strands ('plus' and '-') with careful attention
    to avoid exceeding genome boundaries.
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
        if subject_length < extend_number:  # If the sequence is less than 1000 nt, we'll expand it.
            leftover_length = extend_number - subject_length  # We get the difference between 1000 and the length of the sequence
            leftover_length_halved = leftover_length / 2  # We divide it by 2 to get the half of the difference
            # leftover_length_halved = math.ceil(leftover_length_halved) # We round up the number
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
        else:  # If the sequence is already 1000 nt, we just pass
            pass
        
    return data_input
