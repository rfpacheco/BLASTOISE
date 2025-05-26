import os
import time
import pandas as pd

from modules.files_manager import fasta_creator, columns_to_numeric, end_always_greater_than_start, get_data_sequence
from modules.seq_modifier import sequence_extension
from modules.filters import global_filters_main
from modules.strand_location import set_strand_direction

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def genome_specific_chromosome_main(data_input, main_folder_path, genome_fasta, identity_1, run_phase, word_size, min_length, extend_number, limit_len, coincidence_data=None):
    """
    Executes the main function to process genome-specific chromosomes, extending sequences, performing BLASTn
    analysis, filtering results, and returning a filtered DataFrame.

    Arguments:
        data_input (pd.DataFrame): The input data containing sequence information for processing.
        main_folder_path (str): Path to the folder where output files will be stored.
        genome_fasta (str): Path to the genome FASTA file used for sequence alignment.
        identity_1 (float): The minimum percentage of identity required for BLASTn matches.
        run_phase (int): An integer representing the current processing phase.
        word_size (int): The BLASTn word size configuration for sequence alignment.
        min_length (int): The minimum alignment length for filtering BLASTn results.
        extend_number (int): The number of nucleotides by which the sequences will be extended.
        limit_len (int): The maximum allowable length for the extended sequences.
        coincidence_data (pd.DataFrame, optional): Data that incorporates the coincidence elements between the run "n" and run "n-1". Defaults to None.

    Returns:
        pd.DataFrame: A pandas DataFrame containing filtered results after sequence extension, BLASTn analysis,
        and post-processing.
    """
    from modules.blaster import blastn_blaster  # Delayed import --> to break the circular import. Need to be at the start of function.

    # Output run fase number
    run_phase_extension_path = os.path.join(main_folder_path, f"run_{str(run_phase)}")
    os.makedirs(run_phase_extension_path, exist_ok=True)  # Create the run phase folder
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    if coincidence_data is not None: # If coincidence data exist with actual information.
        print("")
        print('\t\t2.0. Input data information')
        print(f"\t\t\t- Input data row length: {data_input.shape[0]}")
        print(f"\t\t\t- Coincidence data row length: {coincidence_data.shape[0]}")
        # Join the data input coming from run "n" with the coincidence data between run "n" and "n-1"
        data_input = pd.concat([data_input, coincidence_data], ignore_index=True).copy()
        data_input.sort_values(by=['sstrand', 'sseqid', 'sstart'], inplace=True)
        print(f"\t\t\t- Data row length after joining: {data_input.shape[0]}")
    else: # If coincidence_data == None, just keep going with `input_data` as it is.
        pass

    data_to_extend = data_input.copy() # Copy `data_input` to not modify the original data

    print("")
    print(f"\t\t2.1. Sequence extension in both directions to {extend_number} nt:\n",
          f"\t\t\t- Data row length: {data_to_extend.shape[0]}")

    # Extend `data_to_extend` in both directions to `extend_number` nt
    sequences_extended = sequence_extension(data_input=data_to_extend,
                                            genome_fasta=genome_fasta,
                                            extend_number=extend_number,
                                            limit_len=limit_len)
    sequences_extended_fasta_path = os.path.join(run_phase_extension_path, f"run_{extend_number}nt.fasta")  # Path to the output FASTA file
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # Here he have `data_input` with the original coordinates and `sequences_extended` with the extended coordinates.
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print(f"\t\t2.2. Fasta {extend_number} nt file creation:")
    # Get a fasta file from the `extended_sequences`.
    # Save the coordinates from the `extended_sequences` and the original coordinates before the extension in `data_input` to the fasta file.
    sequences_extended_with_sseq = get_data_sequence(sequences_extended, genome_fasta)
    fasta_creator(sequences_extended_with_sseq, sequences_extended_fasta_path, id_names=data_input)
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    # Launch the `fasta_file` formed with the `extended_sequences` to the genomes
    second_blaster = blastn_blaster(query_path=sequences_extended_fasta_path,
                                    dict_path=genome_fasta,
                                    perc_identity=identity_1,
                                    word_size=word_size)
    second_blaster = columns_to_numeric(second_blaster) # Transform to an int type
    second_blaster = end_always_greater_than_start(second_blaster) # Make sure "send" > "sstart"

    toc = time.perf_counter()
    print("")
    print(f"\t\t2.3. BLASTn against genome:\n",
          f"\t\t\t- Data row length: {second_blaster.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    # Removing extended coordinates in `second_blaster` using the `extended_sequences` coordinates.
    second_blaster_not_extended = second_blaster.copy()  # Copy data to not modify the original data

    # Split sequence name to get original coordinates
    ## First split will be the index number
    ## Second split will be the extended coordinates
    ## Third split will be the original coordinates
    split_cols = second_blaster_not_extended['qseqid'].str.split('_', expand=True) # Split sequence identifier into components by hyphen delimiter

    # Note information about EXTENDED coordinates
    second_blaster_not_extended['ext_sseqid'] = split_cols[1].str.split('-').str[0] # Extract the extended sequence ID - take the last part after underscore split
    second_blaster_not_extended['ext_sstart'] = pd.to_numeric(split_cols[1].str.split('-').str[1]) # Convert start position string to numeric, store in the extended start column
    second_blaster_not_extended['ext_send'] = pd.to_numeric(split_cols[1].str.split('-').str[2]) # Convert end position string to numeric, store in the extended end column
    second_blaster_not_extended['ext_sstrand'] = split_cols[1].str.split('-').str[3] # Get strand orientation from the last component, store in the extended strand column

    # Note information about ORIGINAL coordinates
    second_blaster_not_extended['og_sseqid'] = split_cols[2].str.split('-').str[0]
    second_blaster_not_extended['og_sstart'] = pd.to_numeric(split_cols[2].str.split('-').str[1])
    second_blaster_not_extended['og_send'] = pd.to_numeric(split_cols[2].str.split('-').str[2])
    second_blaster_not_extended['og_sstrand'] = split_cols[2].str.split('-').str[3]

    # Remove matches where original and current coordinates are identical
    matches_mask = (second_blaster_not_extended['ext_sseqid'] == second_blaster_not_extended['sseqid']) & \
                   (second_blaster_not_extended['ext_sstart'] == second_blaster_not_extended['sstart']) & \
                   (second_blaster_not_extended['ext_send'] == second_blaster_not_extended['send']) & \
                   (second_blaster_not_extended['ext_sstrand'] == second_blaster_not_extended['sstrand'])

    # noinspection PyUnresolvedReferences
    removed_count = matches_mask.sum() # Count the number of elements removed
    second_blaster_not_extended = second_blaster_not_extended[~matches_mask] # Remove the elements
    print(f"\t\t\t- Removed {removed_count} self-matches from extended sequences")
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Set the strand orientation
    print("")
    print(f"\t\t 2.4. Setting strand orientation")
    print(f"\t\t\t- Data row length: {second_blaster_not_extended.shape[0]}")
    tic = time.perf_counter()
    second_blaster_not_extended_oriented = set_strand_direction(data_input=second_blaster_not_extended)
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print("\t\t2.4. Filtering data:\n")
    filtered_data = global_filters_main(data_input=second_blaster_not_extended_oriented,
                                        genome_fasta=genome_fasta,
                                        writing_path=run_phase_extension_path,
                                        min_length=min_length)

    toc = time.perf_counter()
    print(f"\t\t\t- Data row length: {len(filtered_data)}\n",  # Not .shape[0] in case it's empty
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    return filtered_data  # Returns the data frame