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
            lower_coor = int(element["sstart"])
            upper_coor = int(element["send"])
        else:  # If it's the "-" strand
            lower_coor = int(element["send"])
            upper_coor = int(element["sstart"])
        
        subject_len = upper_coor - lower_coor + 1 # Calculate the nucleotide size of the element
        if subject_len < limit_len:  # If the sequence is less than 1000 nt, there's room to expand it
            lower_coor = lower_coor - extend_number # Extends the lower sequences by `extend_number`
            upper_coor = upper_coor + extend_number # Extends the upper sequence by `extend_number`

            # Take into account that the numbers cannot be negative neither go against the limit length of the sequence
            if lower_coor <= 0:  # Because in BLASTn, 0 does not exist, it starts in 1
                lower_coor = 1

            # `upper_coor` does not need to change, because BLASTn won't extract longer values that the maximum sequence length.
            # but let's calculate, to take into account the sizes.
            cmd = f"blastdbcmd -db {genome_fasta} -entry {element['sseqid']} -outfmt '%l'" # Command to extract the max len
            chrom_max_len = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip() # Extract the max len
            chrom_max_len = int(chrom_max_len) # Convert to int
            if upper_coor > chrom_max_len:
                upper_coor = chrom_max_len
            
            new_subject_len = upper_coor - lower_coor + 1 # Calculate the nucleotide size of the element after extension.
            if new_subject_len < limit_len:  # If the sequence is still less than `limit_len` after the adjustment, continue the extension
               # Get the sequence, but extended
                cmd = f"blastdbcmd -db {genome_fasta} -entry {element['sseqid']} -range {int(lower_coor)}-{int(upper_coor)} -strand {element['sstrand']} -outfmt %s"
                seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

                # Now with the data let's modify the data frame
                # IMPORTANT: it will modify the `data_input` data frame as a "pointer".
                # Fill with NaN values where there's no data
                data_input.loc[index, "pident"] = np.nan
                data_input.loc[index, "length"] = np.nan
                data_input.loc[index, "qstart"] = np.nan
                data_input.loc[index, "qend"] = np.nan

                # Now we can update the coordinates and the sequence. Depending on the strand
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
            else:  # If the new seq is going to reach the `limit_len`. Don't extend it. It won't modify the data.
                pass

    return data_input
