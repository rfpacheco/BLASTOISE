"""
BLASTOISE Module: BLAST Operations and Sequence Processing
=========================================================

This module provides core functionality for the BLASTOISE pipeline, handling BLAST database
creation, sequence alignment, and iterative sequence discovery. It serves as the engine
for identifying repetitive genomic elements through a series of BLAST operations and
data processing steps.

The module contains three main functions:
1. `blastn_dic`: Creates a BLAST database from genome FASTA file
2. `blastn_blaster`: Performs BLASTn alignment and returns results as a DataFrame
3. `repetitive_blaster`: Executes the iterative process of sequence extension, 
   re-alignment, and comparison to discover all instances of repetitive elements

These functions work together to implement the core sequence discovery algorithm
of BLASTOISE, progressively identifying and refining the set of repetitive sequences
in the target genome.

Author: R. Pacheco
"""

import os
import pandas as pd
import subprocess
import time
import shutil
import logging
from datetime import datetime

from .aesthetics import print_message_box
from .seq_identifier import sequence_identifier
from .compare import compare_main
from .files_manager import end_always_greater_than_start, get_data_sequence
from .strand_location import del_last_overlapping_elem
from .genomic_ranges import get_merge_stranded
from blastoise.extra.utils.csv_to_gff import csv_to_gff


def blastn_dic(path_input: str, path_output: str) -> None:
    """
    Create a BLAST-compatible nucleotide database from a FASTA file.

    This function executes the NCBI 'makeblastdb' command-line utility to create a 
    nucleotide database that can be used for subsequent BLAST searches. The database 
    is configured to preserve sequence IDs for proper reference in search results.

    Parameters
    ----------
    path_input : str
        Path to the input FASTA file to be used for building the BLAST database.
    path_output : str
        Path where the database files will be stored.

    Raises
    ------
    Exception
        If the BLAST database creation fails, the error is logged but not raised.

    Notes
    -----
    The function suppresses standard output and error streams from the makeblastdb
    command to avoid cluttering the console. Errors are captured and logged using
    the logging module.
    """
    
    try:
        # "parse_seqids" is used to keep the sequence ID in the output.
        cmd = f"makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        logging.error(f"Error: Blast Dictionary couldn't be created: {e}", exc_info=True)


def blastn_blaster(
        query_path: str,
        dict_path: str,
        perc_identity: float,
        word_size: int = 15
) -> pd.DataFrame:
    """
    Execute BLASTn alignment and return results as a structured DataFrame.

    This function performs a nucleotide BLAST (BLASTn) search using the provided query
    sequences against a pre-built BLAST database. It configures the search with the
    specified identity threshold and word size, then parses the results into a pandas
    DataFrame with standardized column names and data types.

    The function also performs post-processing on the results:
    1. Converts coordinate and e-value columns to appropriate data types
    2. Ensures that 'send' is always greater than 'sstart' for consistent orientation
    3. Calculates sequence length and adds it as a column
    4. Reorders columns for better readability

    Parameters
    ----------
    query_path : str
        Path to the FASTA file containing query nucleotide sequences.
    dict_path : str
        Path to the pre-built BLAST database (created with blastn_dic).
    perc_identity : float
        a Minimum percentage identity threshold for reporting matches (0-100).
    word_size : int, optional
        Size of the word used for seeding alignments in the BLAST algorithm.
        Default is 15.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the parsed and processed BLAST results with columns:
        'qseqid' (query sequence ID), 'sseqid' (subject sequence ID), 
        'sstart' (subject start position), 'send' (subject end position),
        'sstrand' (subject strand), 'evalue' (expectation value),
        'sseq' (aligned subject sequence), and 'len' (length of alignment).

    Raises
    ------
    subprocess.CalledProcessError
        If the BLAST command fails to execute properly.

    Notes
    -----
    The BLAST output format is configured to return only specific fields needed
    for downstream analysis, which may differ from standard BLAST output formats.
    """

    cmd = (
        f"blastn -word_size {word_size} "
        f"-query {query_path} "
        f"-db {dict_path} "
        f"-perc_identity {perc_identity} "
        f"-outfmt '10 qseqid sseqid sstart send sstrand evalue sseq'"
    )
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    data.columns = ['qseqid', 'sseqid', 'sstart', 'send', 'sstrand', 'evalue', 'sseq']
    # Get 'sstart', 'send' to 'int' type
    data[['sstart', 'send']] = data[['sstart', 'send']].astype(int)
    # get 'evalue' as a 'float' type
    data['evalue'] = data['evalue'].astype(float)
    # Make sure 'send' > 'sstart'
    data = end_always_greater_than_start(data)
    # Create 'len' colum
    data['len'] = data['send'] - data['sstart'] + 1
    # Place it between 'sent' and 'sstrand' column
    data = data[['qseqid', 'sseqid', 'sstart', 'send', 'sstrand', 'evalue', 'sseq', 'len']]
    return data


def repetitive_blaster(
        data_input: pd.DataFrame,
        genome_fasta: str,
        folder_path: str,
        numbering: int,
        start_time: str,
        identity_1: float,
        tic_start: float,
        word_size: int,
        min_length: int,
        extend_number: int,
        limit_len: int,
        coincidence_data: pd.DataFrame | None = None,
        mask: pd.DataFrame | None = None,
        n_jobs: int = -1
) -> None:
    """
    Execute the iterative repetitive sequence discovery process.

    This function is the core of the BLASTOISE pipeline, implementing an iterative
    approach to discovering all instances of repetitive elements in a genome. For each
    iteration (run), it performs the following steps:

    1. Sorts and prepares input data
    2. Identifies and extends sequences using sequence_identifier function
    3. Applies optional masking to filter out predefined genomic regions
    4. Compares results with previous runs to identify new sequences
    5. Either terminates if no new sequences are found or continues to the next iteration
       with the newly discovered sequences

    The function maintains detailed logs of its progress, including timing information
    and statistics about the number of sequences found in each iteration. Results from
    each run are saved to disk for traceability and potential reuse.

    Parameters
    ----------
    data_input : pd.DataFrame
        A DataFrame containing sequence data from previous iterations or initial BLAST results.
        Must include columns: 'sseqid', 'sstart', 'send', and 'sstrand'.
    genome_fasta : str 
        Path to the BLAST database of the reference genome.
    folder_path : str
        Directory path for saving intermediate files and results.
    numbering : int
        The identifier for the current iteration/run number.
    start_time : str
        Formatted timestamp marking when the overall analysis began.
    identity_1 : float
        Percentage identity threshold for BLAST searches.
    tic_start : float
        Initial timer value to measure total program execution time.
    word_size : int
        Word size parameter for BLAST searches.
    min_length : int
        Minimum sequence length to retain in the analysis.
    extend_number : int
        Number of nucleotides to extend sequences in both directions.
    limit_len : int
        Length threshold that triggers sequence extension.
    coincidence_data : pd.DataFrame | None, optional
        Sequences found in common between previous iterations. Default is None.
    mask : pd.DataFrame | None, optional
        DataFrame defining regions to exclude from the analysis. Default is None.
    n_jobs : int, optional
        Number of jobs for parallel processing. -1 means using all processors. Default is -1.

    Returns
    -------
    None
        Results are saved to disk rather than returned.

    Notes
    -----
    This function is recursive, calling itself with updated parameters until the
    termination condition (no new sequences found) is met. This design allows the
    algorithm to continue discovering sequences until convergence.

    The function creates several output files:
    - CSV files containing sequence data for each run
    - GFF files for visualization in genome browsers
    - A final consolidated CSV file with all discovered sequences
    """

    # Call the aesthetics function RUN identifier.
    print_message_box('RUN ' + str(numbering))
    tic_main = time.perf_counter()  # Start the timer

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    data_ordered = data_input.sort_values(by=["sseqid", "sstart"])
    toc = time.perf_counter()
    print('1. Initial data:\n',
          f"\t- Data row length: {data_input.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")
    print("")
    terminal_width = shutil.get_terminal_size().columns  # Get the terminal width

    print("")
    print('2. Individual searching and cleaning:')
    tic = time.perf_counter()
    now_time = datetime.now()
    formatted_now_time = now_time.strftime('%Y %B %d at %H:%M')
    print("")
    print(f"{' ' * 7}{'-' * 74}")
    start_time_text = f"Program started: {start_time}"
    end_time_text = f"Program time now: {formatted_now_time}"
    run_text = f"RUN {numbering}"
    print(f"{run_text:>{terminal_width}}")
    print(f"{start_time_text:>{terminal_width}}")
    print(f"{end_time_text:>{terminal_width}}")

    whole_group = sequence_identifier(
        data_input=data_ordered,
        main_folder_path=folder_path,
        genome_fasta=genome_fasta,
        identity_1=identity_1,
        run_phase=numbering,
        coincidence_data=coincidence_data,
        word_size=word_size,
        min_length=min_length,
        limit_len=limit_len,
        extend_number=extend_number,
        n_jobs=n_jobs
    )

    # -----------------------------------------------------------------------------
    # STEP 3: Apply optional masking to filter out predefined genomic regions
    # -----------------------------------------------------------------------------
    if mask is not None:
        from blastoise.modules.filters import remove_masking_zone
        # Remove sequences that overlap with masked regions
        whole_group = remove_masking_zone(whole_group, mask)

    # Report statistics after sequence identification and optional masking
    toc = time.perf_counter()
    print(f"\t- Data row length: {whole_group.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    # STEP 4: Compare results with previous runs
    # -----------------------------------------------------------------------------
    print("")
    print(f"3. Comparison VS Previous Run:")

    if coincidence_data is not None:
        # Not the first run - we have data from previous iterations
        print("")
        print(f"\t- Previous Run data:\n",
              f"\t\t- Coincidence data row length: {coincidence_data.shape[0]}\n",
              f"\t\t- New data row length: {data_input.shape[0]}")

        # Combine previous coincidence data with current input data for comprehensive comparison
        # This ensures we compare against all previously discovered sequences, not just the new ones
        data_input = pd.concat([coincidence_data, data_input], ignore_index=True)
        data_input.sort_values(by=["sseqid", "sstart"], inplace=True)
        print(f"\t\t- Total data row length: {data_input.shape[0]}")
    else:
        # First run - no previous data to compare against
        print(f"\t- Previous Run data:\n",
              f"\t\t- First Run row length: {data_input.shape[0]}")

    # -----------------------------------------------------------------------------
    # Remove overlapping sequences to ensure clean data
    # -----------------------------------------------------------------------------
    # When sequences overlap, keep only the largest one to avoid redundancy
    # whole_group = current run results (n)
    # data_input = combined previous results (n-1)
    whole_group = del_last_overlapping_elem(whole_group)

    # -----------------------------------------------------------------------------
    # Compare current run results with previous runs to identify new sequences
    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    print("")
    print(f"\t- Results in this RUN:")

    # The compare_main function splits the data into three categories:
    # 1. coincidence_data: Sequences found in both current and previous runs
    # 2. new_data: Sequences found only in the current run
    # 3. old_data_exclusive: Sequences found only in previous runs
    coincidence_data, new_data, old_data_exclusive = compare_main(whole_group, data_input)
    toc = time.perf_counter()

    # Report comparison results
    print("")
    print(f"\t\t- Coincidence data from run 'n' and 'n-1': {coincidence_data.shape[0]}\n",
          f"\t\t- New data detected only from run 'n': {new_data.shape[0]}\n",
          f"\t\t- Previous data detected only from 'n-1': {old_data_exclusive.shape[0]}")

    # -----------------------------------------------------------------------------
    # Process sequences from previous runs that weren't found in current run
    # -----------------------------------------------------------------------------
    # Initialize variable to hold sequences below a minimum length threshold
    old_data_exclusive_less_than_min_length = None

    # Calculate sequence lengths for filtering
    old_data_exclusive['len'] = abs(old_data_exclusive['send'] - old_data_exclusive['sstart']) + 1

    # Split sequences from previous runs into two groups based on length
    if not old_data_exclusive.empty and (old_data_exclusive['len'] < min_length).sum() > 0:
        # Extract sequences shorter than the minimum length
        old_data_exclusive_less_than_min_length = old_data_exclusive[old_data_exclusive['len'] < min_length]
        # Keep only sequences meeting or exceeding the minimum length
        old_data_exclusive = old_data_exclusive[old_data_exclusive['len'] >= min_length]

        # Remove length column to avoid conflicts in later concatenation operations
        old_data_exclusive_less_than_min_length.drop(columns=['len'], inplace=True)
        old_data_exclusive.drop(columns=['len'], inplace=True)

        # Report statistics on the filtering
        print(f"\t\t\t- Sequences < {min_length} bp: {old_data_exclusive_less_than_min_length.shape[0]}")
        print(f"\t\t\t- Sequences >= {min_length} bp: {old_data_exclusive.shape[0]}")

    # -----------------------------------------------------------------------------
    # Prepare data for next iteration
    # -----------------------------------------------------------------------------
    # Combine newly discovered sequences with short sequences from previous runs
    # Short sequences are included because they might extend in the next iteration
    if old_data_exclusive_less_than_min_length is not None:
        new_data_and_old = pd.concat([new_data, old_data_exclusive_less_than_min_length], ignore_index=True)
        new_data_and_old.sort_values(by=["sseqid", "sstart"], inplace=True)
        print('\t' * 3 + f"- New data + less than {min_length}: {new_data_and_old.shape[0]}")
    else:
        new_data_and_old = new_data

    # Update coincidence data to include sequences from previous runs that weren't found
    # in the current run but meet the minimum length requirement
    if not coincidence_data.empty and not old_data_exclusive.empty:
        # Add sequences from previous runs to coincidence data
        coincidence_data = pd.concat([coincidence_data, old_data_exclusive], ignore_index=True)
        print(f"\t\t- Coincidence data + Previous data: {coincidence_data.shape[0]}")

        # Merge overlapping sequences and remove duplicates
        coincidence_data = get_merge_stranded(coincidence_data)
        coincidence_data = del_last_overlapping_elem(coincidence_data)
        coincidence_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
        print(f"\t\t\t - After merging: {coincidence_data.shape[0]}")

    # Report timing information
    print(f"\t\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    # STEP 5: Determine whether to terminate or continue to next iteration
    # -----------------------------------------------------------------------------
    if new_data.shape[0] == 0:
        # -----------------------------------------------------------------------------
        # TERMINATION CONDITION: No new sequences found
        # -----------------------------------------------------------------------------
        # Prepare final output by selecting only essential columns
        coincidence_data = coincidence_data[['sseqid', 'sstart', 'send', 'sstrand']].copy()

        # Ensure consistent coordinate orientation (sstart < send)
        coincidence_data = end_always_greater_than_start(coincidence_data)

        # Retrieve the actual sequence data for each genomic region
        coincidence_data = get_data_sequence(coincidence_data, genome_fasta)

        # Save final results to CSV and GFF formats
        final_csv_path = os.path.join(folder_path, "blastoise_df.csv")
        coincidence_data.to_csv(final_csv_path, index=False, header=True, sep=",")
        csv_to_gff(final_csv_path)

        # Report termination and final statistics
        print("")
        print(f"4. Stopping:")
        print(f"\t- No new data found.")
        print(f"\t- BLASTOISE final data row length: {coincidence_data.shape[0]}")

        # Exit the recursive process
        return
    else:
        # -----------------------------------------------------------------------------
        # CONTINUE: New sequences found, prepare for next iteration
        # -----------------------------------------------------------------------------
        # Report timing for current run
        toc_main = time.perf_counter()
        print("")
        print(f"RUN {numbering} finished:\n",
              f"\t- Execution time: {toc_main - tic_main:0.2f} seconds")

        # Save intermediate results for this iteration
        # Combine all sequences for comprehensive tracking
        save_run_file = pd.concat([new_data_and_old, coincidence_data], ignore_index=True)
        save_run_file.sort_values(by=['sseqid', 'sstart'], inplace=True)

        # Remove any overlapping elements to ensure clean data
        save_run_file = del_last_overlapping_elem(save_run_file)

        # Calculate sequence lengths
        save_run_file['len'] = save_run_file['send'] - save_run_file['sstart'] + 1

        # Create a directory for storing run-specific results
        runs_folder = os.path.join(folder_path, "RUNS")
        os.makedirs(runs_folder, exist_ok=True)

        # Save this run's results to CSV and GFF formats
        run_saver_path = os.path.join(runs_folder, f"run_{numbering}.csv")
        save_run_file.to_csv(run_saver_path, sep=",", header=True, index=False)
        csv_to_gff(run_saver_path)

        # Report statistics and file location
        print(f"\t- Coincidence data + 'n' exclusive data + 'n-1' exclusive data: {save_run_file.shape[0]}")
        print(f"\t- Data file saved at {run_saver_path}")

        # -----------------------------------------------------------------------------
        # RECURSIVE CALL: Continue to the next iteration with updated data
        # -----------------------------------------------------------------------------
        # Increment run counter
        numbering += 1

        # Call this function recursively with:
        # - new_data_and_old: Newly discovered sequences to extend in the next iteration
        # - coincidence_data: Accumulated sequences from all previous iterations
        repetitive_blaster(
            data_input=new_data_and_old,
            genome_fasta=genome_fasta,
            folder_path=folder_path,
            numbering=numbering,
            start_time=start_time,
            identity_1=identity_1,
            tic_start=tic_start,
            word_size=word_size,
            min_length=min_length,
            extend_number=extend_number,
            limit_len=limit_len,
            coincidence_data=coincidence_data,
            mask=mask,
            n_jobs=n_jobs
        )
