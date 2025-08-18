"""
BLASTOISE Module: Sequence Identification and Processing
=======================================================

This module provides functionality for identifying and processing genomic sequences
in the BLASTOISE pipeline. It handles sequence extension, BLAST alignment, coordinate
tracking, and filtering to identify relevant genomic regions.

The module contains one main function:
1. `sequence_identifier`: Processes genomic data through a series of steps including
   sequence extension, BLAST alignment, coordinate tracking, and filtering.

This function serves as a critical component in the BLASTOISE pipeline, transforming
raw genomic coordinates into properly extended, aligned, and filtered sequences
that can be used for downstream analysis of repetitive elements.

Author: R. Pacheco
"""

import os
import time
from typing import Optional

import pandas as pd

from .files_manager import fasta_creator, get_data_sequence
from .seq_extension import sequence_extension
from .filters import global_filters_main
from .strand_location import set_strand_direction

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def sequence_identifier(
        data_input: pd.DataFrame,
        main_folder_path: str,
        genome_fasta: str,
        identity_1: float,
        run_phase: int,
        word_size: int,
        min_length: int,
        extend_number: int,
        limit_len: int,
        coincidence_data: Optional[pd.DataFrame] | None = None,
        n_jobs: int = -1
) -> pd.DataFrame:
    """
    Processes genome-specific chromosome data by performing several operations such as sequence
    extensions, BLAST alignments, and data filtering. The function handles optional input for
    coincidence data and manages intermediate and final results in a structured manner.

    This function implements a multistep process to identify and process genomic sequences:
    1. Merges input data with coincidence data from previous runs (if provided)
    2. Extends sequences in both directions to capture flanking regions
    3. Creates FASTA files from the extended sequences
    4. Performs BLAST alignment against reference genome
    5. Extracts and processes coordinate information from alignment results
    6. Sets strand orientation for each sequence
    7. Applies filters to remove unwanted sequences

    Each step is timed and progress information is printed to the console.

    Parameters
    ----------
    data_input: pd.DataFrame
        The primary input data containing sequence information with columns 'sseqid', 
        'sstart', 'send', and 'sstrand'.
    main_folder_path: str
        The root folder path where intermediate and final results will be stored.
    genome_fasta: str
        Path to the genome FASTA file used for sequence processing.
    identity_1: float
        Minimum sequence identity required for BLAST alignment (0-100).
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
        Data frame containing information on coinciding sequences between previous and current runs. 
        Must have the same structure as data_input. Default is None.
    n_jobs: int, optional
        Number of jobs for parallel processing. -1 means using all processors. Default is -1.

    Returns
    -------
    pd.DataFrame
        A filtered data frame representing the processed result of genome-specific
        chromosome sequences while adhering to all defined filters and extensions.

    Raises
    ------
    FileNotFoundError
        If the genome FASTA file cannot be found or accessed.
    ValueError
        If required columns are missing from the input data frames.
    """

    # Import blastn_blaster here to avoid circular import issues
    from blastoise.modules.blaster import run_blastn_alignment  # Delayed import to break the circular import

    # -----------------------------------------------------------------------------
    # STEP 1: Setup output directory for this run phase
    # -----------------------------------------------------------------------------
    run_phase_extension_path = os.path.join(main_folder_path, f"run_{str(run_phase)}")
    os.makedirs(run_phase_extension_path, exist_ok=True)  # Create the run phase folder

    # -----------------------------------------------------------------------------
    # STEP 2: Process and merge input data with coincidence data (if provided)
    # -----------------------------------------------------------------------------
    if coincidence_data is not None:  # If coincidence data exists with actual information
        print("")
        print('\t\t2.1. Input data information')
        print(f"\t\t\t- Input data row length: {data_input.shape[0]}")
        print(f"\t\t\t- Coincidence data row length: {coincidence_data.shape[0]}")

        # Join the data input from current run "n" with the coincidence data between runs "n" and "n-1"
        # This ensures we process both new sequences and those that overlap with previous runs
        data_input = pd.concat([data_input, coincidence_data], ignore_index=True).copy()
        data_input.sort_values(by=['sstrand', 'sseqid', 'sstart'], inplace=True)
        print(f"\t\t\t- Data row length after joining: {data_input.shape[0]}")
    # If coincidence_data is None, we just continue with the original data_input

    # Create a copy of the input data to avoid modifying the original
    data_to_extend = data_input.copy()

    # -----------------------------------------------------------------------------
    # STEP 3: Extend sequences in both directions to capture flanking regions
    # -----------------------------------------------------------------------------
    print("")
    print(f"\t\t2.2. Sequence extension in both directions to {extend_number} nt:\n",
          f"\t\t\t- Data row length: {data_to_extend.shape[0]}")

    tic = time.perf_counter()
    # Extend each sequence by the specified number of nucleotides in both directions
    # This helps capture flanking regions that might be part of the same genomic element
    sequences_extended = sequence_extension(
        data_to_extend,
        genome_fasta,
        extend_number,
        limit_len,
        n_jobs=n_jobs
    )

    # Calculate the length of each extended sequence
    sequences_extended['len'] = sequences_extended['send'] - sequences_extended['sstart'] + 1

    # Define the path where the FASTA file with extended sequences will be saved
    sequences_extended_fasta_path = os.path.join(run_phase_extension_path, f"run_{extend_number}nt.fasta")

    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # At this point we have `data_input` with the original coordinates and `sequences_extended` with the extended coordinates

    # -----------------------------------------------------------------------------
    # STEP 4: Create FASTA file from extended sequences for BLAST analysis
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print(f"\t\t2.3. Fasta {extend_number} nt file creation:")

    # Retrieve the actual sequence data for each extended coordinate set
    sequences_extended_with_sseq = get_data_sequence(sequences_extended, genome_fasta)

    # Create a FASTA file that includes both extended coordinates and original coordinates
    # This allows us to track the relationship between extended and original sequences
    fasta_creator(sequences_extended_with_sseq, sequences_extended_fasta_path, id_names=data_input)

    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # STEP 5: Perform BLAST alignment against the reference genome
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print(f"\t\t2.4. BLASTn against genome:")

    # Execute BLAST search using the extended sequences against the reference genome
    # This identifies all regions in the genome that match our extended sequences
    second_blaster = run_blastn_alignment(
        sequences_extended_fasta_path,
        genome_fasta,
        identity_1,
        word_size
    )

    # Filter out sequences that are shorter than the minimum length threshold
    second_blaster = second_blaster[second_blaster['len'] >= min_length]

    toc = time.perf_counter()
    print(f"\t\t\t- Data row length: {second_blaster.shape[0]}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    # STEP 6: Extract and process coordinate information from BLAST results
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print(f"\t\t2.5. Setting all coordinates information")

    # Create a copy of the BLAST results to avoid modifying the original data
    second_blaster_not_extended = second_blaster.copy()

    # Parse the sequence identifiers to extract coordinate information
    # The sequence IDs in the FASTA file were formatted to contain both extended and original coordinates
    # Format: index_extendedCoordinates_originalCoordinates
    # Where coordinates are formatted as: seqid-start-end-strand
    split_cols = second_blaster_not_extended['qseqid'].str.split('_', expand=True)

    # Extract and convert EXTENDED coordinates from the second part of the split
    second_blaster_not_extended['ext_sseqid'] = split_cols[1].str.split('-').str[0]  # Sequence ID
    second_blaster_not_extended['ext_sstart'] = pd.to_numeric(split_cols[1].str.split('-').str[1])  # Start position
    second_blaster_not_extended['ext_send'] = pd.to_numeric(split_cols[1].str.split('-').str[2])  # End position
    second_blaster_not_extended['ext_sstrand'] = split_cols[1].str.split('-').str[3]  # Strand orientation

    # Extract and convert ORIGINAL coordinates from the third part of the split
    second_blaster_not_extended['og_sseqid'] = split_cols[2].str.split('-').str[0]  # Sequence ID
    second_blaster_not_extended['og_sstart'] = pd.to_numeric(split_cols[2].str.split('-').str[1])  # Start position
    second_blaster_not_extended['og_send'] = pd.to_numeric(split_cols[2].str.split('-').str[2])  # End position
    second_blaster_not_extended['og_sstrand'] = split_cols[2].str.split('-').str[3]  # Strand orientation

    # Remove self-matches (sequences that match their own extended coordinates)
    # This prevents counting the same sequence twice in different forms
    matches_mask = (second_blaster_not_extended['ext_sseqid'] == second_blaster_not_extended['sseqid']) & \
                   (second_blaster_not_extended['ext_sstart'] == second_blaster_not_extended['sstart']) & \
                   (second_blaster_not_extended['ext_send'] == second_blaster_not_extended['send']) & \
                   (second_blaster_not_extended['ext_sstrand'] == second_blaster_not_extended['sstrand'])

    # Count and remove the self-matches
    removed_count = matches_mask.sum()
    second_blaster_not_extended = second_blaster_not_extended[~matches_mask]

    print(f"\t\t\t- Removed {removed_count} self-matches from extended sequences")
    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # STEP 7: Set strand orientation for each sequence
    # -----------------------------------------------------------------------------
    print("")
    print(f"\t\t 2.6. Setting strand orientation")
    print(f"\t\t\t- Data row length: {second_blaster_not_extended.shape[0]}")
    tic = time.perf_counter()

    # Determine the correct strand orientation for each sequence
    # This ensures consistent strand assignment and handles overlapping sequences properly
    second_blaster_not_extended_oriented = set_strand_direction(
        second_blaster_not_extended,
        run_phase,
        main_folder_path,
        n_jobs=n_jobs
    )

    toc = time.perf_counter()
    print(f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    # STEP 8: Apply filters to remove unwanted sequences
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print("\t\t2.7. Filtering data:")

    # Apply global filters to remove sequences that don't meet criteria
    # This includes length filters and other quality control measures
    filtered_data = global_filters_main(
        second_blaster_not_extended_oriented,
        min_length
    )

    toc = time.perf_counter()
    # Use len() instead of .shape[0] to handle empty DataFrames gracefully
    print(f"\t\t\t- Data row length: {len(filtered_data)}\n",
          f"\t\t\t- Execution time: {toc - tic:0.2f} seconds")

    # Return the final filtered data frame with properly oriented and filtered sequences
    return filtered_data
