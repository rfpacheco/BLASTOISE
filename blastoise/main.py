"""
BLASTOISE: A Tool for Repetitive Sequence Discovery in Genomic Data
====================================================================

This script serves as the main entry point for the BLASTOISE pipeline, designed to
identify and analyze repetitive sequences, such as SIDERs (Short Interspersed
Degenerated Retroposons), in genomic data.

The pipeline executes the following major steps:
1.  Sets up a dedicated working directory for the analysis.
2.  Performs an initial BLASTn search to identify seed sequences from the
    provided input data against a reference genome.
3.  Merges the initial findings to consolidate overlapping genomic regions.
4.  Enters an iterative process (`repetitive_blaster`) where sequences are
    extended, re-aligned, and compared against previous findings until no new
    sequences are discovered.
5.  Optionally applies a mask to filter out predefined genomic regions.
6.  Cleans up intermediate files and saves the final results to
    'blastoise_df.csv' in the main output directory.

Author: R. Pacheco
Version: 0.4.2
License: MIT
"""

import argparse
import os
import shutil
import time
import subprocess
import tempfile
from datetime import datetime
from typing import Tuple
# noinspection PyPackageRequirements
import pandas as pd

from blastoise.modules.blaster import blastn_dic, blastn_blaster
from blastoise.modules.aesthetics import print_message_box, blastoise_art
from blastoise.modules.genomic_ranges import get_merge_stranded, get_overlapping_info, get_interval_overlap
from blastoise.modules.seq_extension import sequence_extension
from blastoise.modules.strand_location import match_data_and_remove


def parse_arguments() -> argparse.Namespace:
    """
    Parses command-line arguments for the BLASTOISE script.

    This function sets up the argument parser with all required and optional parameters
    for running the BLASTOISE pipeline, including input/output paths, sequence identity
    thresholds, and various filtering parameters.

    Returns
    -------
    argparse.Namespace
        An object containing the parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        prog='BLASTOISE',
        description='BLASTOISE: A Tool for Repetitive Sequence Discovery in Genomic Data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-d', '--data', type=str, required=True, 
                        help='Path to the input data file.')
    parser.add_argument('-g', '--genome', type=str, required=True, 
                        help='Path to the reference genome file.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path for the output working directory. The specified directory will be created if '
                             'it does not exist.')
    parser.add_argument('-i', '--identity', type=int, default=60, 
                        help='Identity percentage for the first BLASTn step.')
    parser.add_argument('-ws', '--word_size', type=int, default=15, 
                        help='Word size for BLASTn.')
    parser.add_argument('-min', '--min_length', type=int, default=100, 
                        help='Minimum sequence length for filtering.')
    parser.add_argument('-ext', '--extend', type=int, default=100, 
                        help='Number of nucleotides for sequence extension.')
    parser.add_argument('-lim', '--limit', type=int, default=1000, 
                        help='Length limit to trigger sequence extension.')
    parser.add_argument('-m', '--mask', type=str, required=False, 
                        help='Path to an optional mask file (CSV format).')
    parser.add_argument('-j', '--jobs', type=int, default=-1, 
                        help='Number of jobs for parallel processing. -1 means using all processors.')

    return parser.parse_args()


def setup_workspace(
        output_dir: str, 
        data_file: str, 
        genome_file: str
) -> Tuple[str, str, str]:
    """
    Creates the main output directory and a subdirectory for original data,
    then copies the input files into it.

    This function prepares the workspace by creating the necessary directory structure
    and copying the input files to a dedicated location within the output directory.
    This ensures that the original files are preserved and all operations are performed
    on copies.

    Parameters
    ----------
    output_dir : str
        The root directory for all output.
    data_file : str
        Path to the user's input data file.
    genome_file : str
        Path to the user's reference genome file.

    Returns
    -------
    Tuple[str, str, str]
        A tuple containing the path to the created output directory, the new path 
        for the data file, and the new path for the genome file.

    Raises
    ------
    FileNotFoundError
        If the source data file or genome file cannot be found.
    """
    # Create the main output directory
    os.makedirs(output_dir, exist_ok=True)
    print(f"{'.'*20} Output folder created at: {output_dir}")

    # Create a subdirectory for original data and copy files
    original_data_dir = os.path.join(output_dir, 'original_data')
    os.makedirs(original_data_dir, exist_ok=True)

    # Define destination paths
    dest_data_path = os.path.join(original_data_dir, os.path.basename(data_file))
    dest_genome_path = os.path.join(original_data_dir, os.path.basename(genome_file))

    # Copy files, exiting if a source file doesn't exist
    try:
        shutil.copy(data_file, dest_data_path)
        shutil.copy(genome_file, dest_genome_path)
    except FileNotFoundError as e:
        print(f"Error: Source file not found - {e}")
        exit(1)

    return output_dir, dest_data_path, dest_genome_path


def run_initial_blast(
    data_path: str,
    genome_path: str,
    output_dir: str,
    identity: int,
    word_size: int
) -> Tuple[pd.DataFrame, str]:
    """
    Performs the initial BLASTn search and merges the results.

    This function creates a BLAST database from the reference genome, performs the
    initial BLASTn search using the provided parameters, and merges overlapping
    regions in the results using PyRanges.

    Parameters
    ----------
    data_path : str
        Path to the input data file.
    genome_path : str
        Path to the reference genome file.
    output_dir : str
        Path to the main output directory.
    identity : int
        Identity percentage threshold for BLASTn matches.
    word_size : int
        Word size parameter for BLASTn algorithm.

    Returns
    -------
    Tuple[pd.DataFrame, str]
        A tuple containing the merged BLASTn results as a DataFrame and 
        the path to the created BLAST database.
    """
    # Create BLASTn database
    blast_db_dir = os.path.join(output_dir, 'blast_database')
    os.makedirs(blast_db_dir, exist_ok=True)
    blast_db_path = os.path.join(blast_db_dir, os.path.basename(genome_path))
    blastn_dic(path_input=genome_path, path_output=blast_db_path)

    # Run first BLASTn
    print_message_box(message='First BLASTn step initiated')
    tic = time.perf_counter()
    initial_blast_results = blastn_blaster(
        query_path=data_path,
        dict_path=blast_db_path,
        perc_identity=identity,
        word_size=word_size
    )
    # Keep only essential columns
    initial_blast_results = initial_blast_results[['qseqid', 'sseqid', 'sstart', 'send', 'sstrand']]
    toc = time.perf_counter()
    print(f"1. Initial data:\n"
          f"\t- Data row length: {len(initial_blast_results)}\n"
          f"\t- Execution time: {toc - tic:0.2f} seconds")

    # Merge results using PyRanges
    print('\t- Filtering and merging data:')
    tic = time.perf_counter()
    merged_results = get_merge_stranded(data_input=initial_blast_results)
    toc = time.perf_counter()
    print(f"\t\t- Merged data row length: {len(merged_results)}\n"
          f"\t\t- Execution time: {toc - tic:0.2f} seconds")

    return merged_results, blast_db_path


def repetitive_sider_searcher(
        data_input: pd.DataFrame,
        genome_path: str,
        extend_number: int,
        word_size: int,
        identity: int,
        min_length: int,
        limit_len: int,
) -> pd.DataFrame:

    # Extend data using the recursive method
    # -------------------------------------------------------------------
    # STEP 1: Extend the input data using the recursive method
    # -------------------------------------------------------------------
    data_extended = sequence_extension(
        data_input=data_input,
        genome_fasta=genome_path,
        extend_number=extend_number,
        limit_len=limit_len,
        identity=identity,
        word_size=word_size,
        min_length=min_length,
        n_jobs=-1
    )

    # -------------------------------------------------------------------
    # STEP 2: Check for the existence of overlapping data
    # -------------------------------------------------------------------
    # Check for overlapping data
    # noinspection PyArgumentList
    overlapping_info = get_overlapping_info(data_extended)
    ## Print overlapping information
    same_strand = len(overlapping_info.get("same_strand", []))
    opposite_strand = len(overlapping_info.get("opposite_strand", []))
    print(f"Overlapping elements found: same strand = {same_strand}, different strand = {opposite_strand}")

    # -------------------------------------------------------------------
    # STEP 3: Launch second Blast to find new elements
    # -------------------------------------------------------------------
    iteration = 1
    are_there_new_elems = True
    while are_there_new_elems:
        print(f"Iteration {iteration}:")
        iteration += 1
        new_elems = pd.DataFrame()  # Will hold the newly discovered elements. Will be used to reject overlaps
        for i, row in data_extended.iterrows():
            # -------------------------------------------------------------------
            # STEP 3.1: Perform the BLASTn to find new elements
            # -------------------------------------------------------------------
            # Create a temporary FAST file with the sequence
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
                temp_fasta.write(f">{row.i}\n{row.sequence}\n")
                temp_fasta_path = temp_fasta.name

            # Perform the BLAST search
            blast_results = blastn_blaster(
                query_path=temp_fasta_path,
                dict_path=genome_path,
                perc_identity=identity,
                word_size=word_size
            )

            # -------------------------------------------------------------------
            # STEP 3.2: Filter the data
            # -------------------------------------------------------------------
            # Remove elements with a len < min_length
            blast_results = blast_results[blast_results['len'] >= min_length]

            # Remove sequences that overlap with `data_extended` or `new_elems`
            ## With `data_extended`
            rm_from_data_extended = get_interval_overlap(blast_results, data_extended)
            match_data_and_remove(blast_results, rm_from_data_extended)

            ## With `new_elems`  # TODO: I am not so sure about this one
            rm_from_new_elems = get_interval_overlap(blast_results, new_elems)
            match_data_and_remove(blast_results, rm_from_new_elems)

            # -------------------------------------------------------------------
            # STEP 3.3: Add data to the collection of `new_elems`
            # -------------------------------------------------------------------
            if not blast_results.empty:
                new_elems = pd.concat([new_elems, blast_results], ignore_index=True)

        # -------------------------------------------------------------------
        # STEP 4a: No new elements found
        # -------------------------------------------------------------------
        # If nothing was found, exit while loop
        if new_elems.empty:
            are_there_new_elems = False
        else:
            # -------------------------------------------------------------------
            # STEP 4b: New elements found
            # -------------------------------------------------------------------
            new_elems_extended = sequence_extension(
                data_input=new_elems,
                genome_fasta=genome_path,
                extend_number=extend_number,
                limit_len=limit_len,
                identity=identity,
                word_size=word_size,
                min_length=min_length,
                n_jobs=-1
            )
            # Check overlapping data inside `new_elems_extended`
            overlapping_info = get_overlapping_info(new_elems_extended)
            same_strand = len(overlapping_info.get("same_strand", []))
            opposite_strand = len(overlapping_info.get("opposite_strand", []))
            print(f"Overlapping elements found: same strand = {same_strand}, different strand = {opposite_strand}")

            # Check overlapping data of `new_elems_extended` with `data_extended`
            overlapping_info = get_overlapping_info(new_elems_extended, data_extended)
            same_strand = len(overlapping_info.get("same_strand", []))
            opposite_strand = len(overlapping_info.get("opposite_strand", []))
            print(f"Overlapping elements found: same strand = {same_strand}, different strand = {opposite_strand}")

            # Merge newly extended elems with the accumulated dataset
            total_data = pd.concat([data_extended, new_elems_extended], ignore_index=True)
            total_data.sort_values(by=['sseqid', 'sstart'], inplace=True)  # Sort by 'sseqid' and 'sstart'
            total_data.reset_index(drop=True, inplace=True)  # Reset index

            # Replace `data_extended` that entered the while loop with our new `total_data
            data_extended = total_data.copy()

    # -------------------------------------------------------------------
    # STEP 5: Outside while-loop
    # -------------------------------------------------------------------
    print("Repetitive SIDE-eR search completed.")
    return data_extended


def finalize_results(output_dir: str, df: pd.DataFrame) -> None:
    """
    Saves a DataFrame as a CSV file in the specified output directory.

    Parameters
    ----------
    output_dir : str
        The directory where the CSV file will be saved.
    df : pd.DataFrame
        The DataFrame to be saved.
    """
    output_file = os.path.join(output_dir, "blastoise_df.csv")
    df.to_csv(output_file, index=False)


def main() -> None:
    """
    Main function to orchestrate the BLASTOISE pipeline.

    This function coordinates the entire BLASTOISE workflow by:
    1. Parsing command-line arguments
    2. Setting up the workspace
    3. Running the initial BLAST search
    4. Applying optional masking
    5. Executing the iterative repetitive sequence discovery process
    6. Finalizing and organizing results
    7. Reporting execution statistics

    The function handles exceptions and ensures proper cleanup of resources.
    """
    args = parse_arguments()
    start_time = datetime.now()
    tic_main = time.perf_counter()
    formatted_start_time = start_time.strftime('%Y %B %d at %H:%M')
    print(f"{'.'*20} Program started: {formatted_start_time}")

    try:
        # 1. Setup workspace and copy input files
        output_dir, data_path, genome_path = setup_workspace(
            output_dir=os.path.expanduser(args.output),
            data_file=os.path.expanduser(args.data),
            genome_file=os.path.expanduser(args.genome)
        )

        # 2. Run initial BLAST and process results
        initial_data, blast_db_path = run_initial_blast(
            data_path=data_path,
            genome_path=genome_path,
            output_dir=output_dir,
            identity=args.identity,
            word_size=args.word_size
        )

        # 3. Run extension algorithm
        extended_data = sequence_extension(
            initial_data,
            blast_db_path,
            extend_number=args.extend,
            limit_len=args.limit,
            identity=args.identity,
            word_size=args.word_size,
            min_length=args.min_length,
            extension_direction='both',
            n_jobs=args.jobs
        )

        # 4. Run the iterative part
        final_data = repetitive_sider_searcher(
            data_input=extended_data,
            genome_path=genome_path,
            extend_number=args.extend,
            word_size=args.word_size,
            identity=args.identity,
            min_length=args.min_length,
            limit_len=args.limit
        )

        # 5. Finalize and clean up
        finalize_results(output_dir, final_data)

    except Exception as e:
        print_message_box(f"An unexpected error occurred: {e}")
        exit(1)

    # 6. Print final summary
    toc_main = time.perf_counter()
    end_time = datetime.now()
    formatted_end_time = end_time.strftime("%Y %B %d at %H:%M")
    final_data_path = os.path.join(output_dir, "blastoise_df.csv")

    print_message_box(message="END OF THE PROGRAM")
    print(f"\t- Total execution time: {toc_main - tic_main:0.2f} seconds\n"
          f"\t- Program started: {formatted_start_time}\n"
          f"\t- Program ended: {formatted_end_time}\n"
          f"\t- Blastoise final file saved at: {final_data_path}")

    # Set file permissions for the output directory
    try:
        subprocess.run(["chmod", "-R", "a+w", output_dir], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Warning: Could not set final permissions on {output_dir}. Error: {e}")


if __name__ == "__main__":
    """
    Entry point for the BLASTOISE application.

    When this script is executed directly (not imported as a module),
    the main function is called to run the BLASTOISE pipeline.
    """
    main()
    print_message_box("BLASTING OUT!")
    blastoise_art()
