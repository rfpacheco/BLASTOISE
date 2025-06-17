import os
import time
from typing import Optional

import pandas as pd

from modules.files_manager import fasta_creator, get_data_sequence
from modules.seq_modifier import sequence_extension
from modules.filters import global_filters_main
from modules.strand_location import set_strand_direction

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def genome_specific_chromosome_main(
        data_input: pd.DataFrame,
        main_folder_path: str,
        genome_fasta: str,
        identity_1: float,
        run_phase: int,
        word_size: int,
        min_length: int,
        extend_number: int,
        limit_len: int,
        coincidence_data: Optional[pd.DataFrame] | None = None
) -> pd.DataFrame:
    """
    Processes genome-specific chromosome data by performing several operations such as sequence
    extensions, BLAST alignments, and data filtering. The function handles optional input for
    coincidence data and manages intermediate and final results in a structured manner.

    Parameters
    ----------
    data_input: pd.DataFrame
        The primary input data containing sequence information.
    main_folder_path: str
        The root folder path where intermediate and final results will be stored.
    genome_fasta: str
        Path to the genome FASTA file used for sequence processing.
    identity_1: float
        Minimum sequence identity required for BLAST alignment.
    run_phase: int
        Current phase of the processing pipeline for directory and log management.
    word_size: int
        Word size parameter for BLAST alignment.
    min_length: int
        Minimum sequence length required to pass filters.
    extend_number: int
        Number of nucleotides by which sequences are extended in both directions.
    limit_len: int
        Maximum allowable sequence length for extended sequences.
    coincidence_data: pd.DataFrame | None, optional
        Data frame containing information on coinciding sequences between previous and current runs. Default is None.

    Returns
    -------
    pd.DataFrame
        A filtered data frame representing the processed result of genome-specific
        chromosome sequences while adhering to all defined filters and extensions.
    """
    from modules.blaster import blastn_blaster  # Delayed import --> to break the circular import. Need to be at the start of function.

    # Output run fase number
    run_phase_extension_path = os.path.join(main_folder_path, f"run_{str(run_phase)}")
    os.makedirs(run_phase_extension_path, exist_ok=True)  # Create the run phase folder
    # -----------------------------------------------------------------------------
    if coincidence_data is not None: # If coincidence data exist with actual information.
        print("")
        print('\t\t2.1. Input data information')
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
    print(f"\t\t2.2. Sequence extension in both directions to {extend_number} nt:\n",
          f"\t\t\t- Data row length: {data_to_extend.shape[0]}")

    tic = time.perf_counter()
    # Extend `data_to_extend` in both directions to `extend_number` nt
    sequences_extended = sequence_extension(
        data_to_extend,
        genome_fasta,
        extend_number,
        limit_len
    )
    sequences_extended['len'] = sequences_extended['send'] - sequences_extended['sstart'] + 1
    sequences_extended_fasta_path = os.path.join(run_phase_extension_path, f"run_{extend_number}nt.fasta")  # Path to the output FASTA file
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # Here he have `data_input` with the original coordinates and `sequences_extended` with the extended coordinates.
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print(f"\t\t2.3. Fasta {extend_number} nt file creation:")
    # Get a fasta file from the `extended_sequences`.
    # Save the coordinates from the `extended_sequences` and the original coordinates before the extension in `data_input` to the fasta file.
    sequences_extended_with_sseq = get_data_sequence(sequences_extended, genome_fasta)
    fasta_creator(sequences_extended_with_sseq, sequences_extended_fasta_path, id_names=data_input)
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    # Launch the `fasta_file` formed with the `extended_sequences` to the genomes
    print("")
    print(f"\t\t2.4. BLASTn against genome:")
    second_blaster = blastn_blaster(
        sequences_extended_fasta_path,
        genome_fasta,
        identity_1,
        word_size
    )
    # Get elems only with a seq length >= `min_length`
    second_blaster = second_blaster[second_blaster['len'] >= min_length]
    toc = time.perf_counter()
    print(f"\t\t\t- Data row length: {second_blaster.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print(f"\t\t2.5. Setting all coordinates information")
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
    print(f"\t\t 2.6. Setting strand orientation")
    print(f"\t\t\t- Data row length: {second_blaster_not_extended.shape[0]}")
    tic = time.perf_counter()
    second_blaster_not_extended_oriented = set_strand_direction(
        second_blaster_not_extended,
        run_phase,
        main_folder_path
    )
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print("\t\t2.7. Filtering data:")
    filtered_data = global_filters_main(
        second_blaster_not_extended_oriented,
        min_length
    )

    toc = time.perf_counter()
    print(f"\t\t\t- Data row length: {len(filtered_data)}\n",  # Not .shape[0] in case it's empty
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    return filtered_data  # Returns the data frame