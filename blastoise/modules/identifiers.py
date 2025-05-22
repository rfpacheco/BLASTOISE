import os
import time
import pandas as pd

from modules.files_manager import fasta_creator, columns_to_numeric
from modules.seq_modifier import sequence_extension
from modules.filters import global_filters_main

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def genome_specific_chromosome_main(data_input, main_folder_path, genome_fasta, identity_1, run_phase, word_size, min_length, extend_number, limit_len, coincidence_data=None):
    """
    Processes genomic data by extending sequences, creating FASTA files, and filtering BLASTn results. The function handles
    genomic sequence analysis, including extending sequences to a defined length, executing BLASTn comparisons, and filtering
    results based on provided parameters.

    Parameters:
        data_input (pd.DataFrame): Input data detailing sequences and relevant genomic annotations.
        main_folder_path (str): The main directory path used to save outputs and intermediate files.
        genome_fasta (str): Path to the reference genome FASTA file used for BLASTn queries.
        identity_1 (float): Minimum percentage identity threshold for BLASTn comparisons.
        run_phase (int): Indicator for the current run phase used in naming folders.
        word_size (int): Word size parameter for BLASTn, influencing sequence alignment sensitivity.
        min_length (int): Minimum alignment length required to retain BLASTn results.
        extend_number (int): Number of nucleotides to extend each sequence.
        limit_len (int): Maximum length of the sequence.
        coincidence_data (pd.DataFrame, optional): Previously computed data to be combined with the current input. Defaults to None.

    Returns:
        pd.DataFrame: Filtered data from BLASTn results after sequence extension and processing.
    """
    from modules.blaster import blastn_dic, blastn_blaster  # Delayed import --> to break the circular import. Need to be at the start of function.

    run_phase_extension_path = os.path.join(main_folder_path, f"run_{str(run_phase)}")
    os.makedirs(run_phase_extension_path, exist_ok=True)  # Folder
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    # Extend a sequence to 1000 nt. Saved into a pandas' Data Frame. It modifies the original data_input
    # If coincidence_data is not None. We add the previous data to the new one, so we don't lose the previous data.
    if coincidence_data is not None:
        data_input = pd.concat([data_input, coincidence_data], ignore_index=True).copy()
        data_input.sort_values(by=["sstrand", "sseqid", "sstart"], inplace=True)
    else:
        pass

    sequences_extended = sequence_extension(data_input=data_input,
                                            genome_fasta=genome_fasta,
                                            extend_number=extend_number,
                                            limit_len=limit_len)
    sequences_extended_fasta_path = os.path.join(run_phase_extension_path, f"run_{extend_number}nt.fasta")  # Path to the output FASTA file
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.1. Sequence extension to {extend_number} nt:\n",
          f"\t\t\t- Data row length: {sequences_extended.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------s
    # If coincidence_data is not None. We add the previous data to the new one, so we don't lose the previous data.
    # if coincidence_data is not None:
    #     sequences_extended = pd.concat([sequences_extended, coincidence_data], ignore_index=True).copy()
    #     sequences_extended.sort_values(by=["sstrand", "sseqid", "sstart"], inplace=True) # Sort the data frame by the start coordinate
    # else:
    #     pass
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    fasta_creator(sequences_extended, sequences_extended_fasta_path)
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.2. Fasta {extend_number} nt file creation:\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")  
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    second_blaster = blastn_blaster(query_path=sequences_extended_fasta_path,
                                    dict_path=genome_fasta,
                                    perc_identity=identity_1,
                                    word_size=word_size)
    second_blaster = columns_to_numeric(second_blaster, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    toc = time.perf_counter()
    print("")
    print(f"\t\t2.3. BLASTn against genome:\n",
          f"\t\t\t- Data row length: {second_blaster.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Removing extended coordinates in `second_blaster` from the data of `sequence_1000`
    second_blaster_not_extended = second_blaster.copy()  # Copy data from the blaster
    # Split sequence name to get original coordinates
    split_cols = second_blaster_not_extended['qseqid'].str.split('-', expand=True)
    second_blaster_not_extended['og_sseqid'] = split_cols[0].str.split('_').str[-1]
    second_blaster_not_extended['og_sstart'] = pd.to_numeric(split_cols[1])
    second_blaster_not_extended['og_send'] = pd.to_numeric(split_cols[2])
    second_blaster_not_extended['og_sstrand'] = split_cols[3]

    # Remove matches where original and current coordinates are identical
    matches_mask = (second_blaster_not_extended['og_sseqid'] == second_blaster_not_extended['sseqid']) & \
                   (second_blaster_not_extended['og_sstart'] == second_blaster_not_extended['sstart']) & \
                   (second_blaster_not_extended['og_send'] == second_blaster_not_extended['send']) & \
                   (second_blaster_not_extended['og_sstrand'] == second_blaster_not_extended['sstrand'])

    removed_count = matches_mask.sum()
    second_blaster_not_extended = second_blaster_not_extended[~matches_mask]
    print(f"\t\t\t- Removed {removed_count} self-matches from extended sequences")

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    filtered_data = global_filters_main(data_input=second_blaster_not_extended,
                                        genome_fasta=genome_fasta,
                                        writing_path=run_phase_extension_path,
                                        min_length=min_length)

    toc = time.perf_counter()
    print("")
    print("\t\t2.4. Filtering BLASTn against genome:\n",
          f"\t\t\t- Data row length: {len(filtered_data)}\n",  # Not .shape[0] in case it's empty
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    return filtered_data  # Returns the data frame