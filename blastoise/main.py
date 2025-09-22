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
import warnings

# Silence the pkg_resources deprecation warning emitted by sorted_nearest
warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    module=r"sorted_nearest(\.__init__)?",
    message=r"pkg_resources is deprecated as an API\."
)

import argparse
import os
import shutil
import time
import subprocess
import tempfile
from datetime import datetime
from typing import Tuple, Optional
import pandas as pd
import logging

logger = logging.getLogger(__name__)  # Create logger instance for the current file

from blastoise.modules.aesthetics import print_message_box, blastoise_art
from blastoise.modules.blaster import create_blast_database, run_blastn_alignment
from blastoise.modules.filters import match_data_and_remove
from blastoise.modules.formats import format_output_dataframe, write_gff_from_formatted
from blastoise.modules.genomic_ranges import fetch_overlapping_intervals, merge_overlapping_intervals
from blastoise.modules.seq_extension import sequence_extension


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
    parser.add_argument('-ws', '--word_size', type=int, default=11,
                        help='Word size for BLASTn.')
    parser.add_argument('-min', '--min_length', type=int, default=100, 
                        help='Minimum sequence length for filtering.')
    parser.add_argument('-ext', '--extend', type=int, default=100, 
                        help='Number of nucleotides for sequence extension.')
    parser.add_argument('-lim', '--limit', type=int, default=1000, 
                        help='Length limit to trigger sequence extension.')
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
    print(f"- Output folder created at: {output_dir}")

    # Create a subdirectory for original data and copy files
    original_data_dir: str = os.path.join(output_dir, 'original_data')
    os.makedirs(original_data_dir, exist_ok=True)

    # Define destination paths
    dest_data_path: str = os.path.join(original_data_dir, os.path.basename(data_file))
    dest_genome_path: str = os.path.join(original_data_dir, os.path.basename(genome_file))

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
    blast_db_dir: str = os.path.join(output_dir, 'blast_database')
    os.makedirs(blast_db_dir, exist_ok=True)
    blast_db_path: str = os.path.join(blast_db_dir, os.path.basename(genome_path))
    create_blast_database(path_input=genome_path, path_output=blast_db_path)

    # Run first BLASTn
    print_message_box(message='First BLASTn step initiated')
    tic: float = time.perf_counter()
    initial_blast_results: pd.DataFrame = run_blastn_alignment(
        query_path=data_path,
        dict_path=blast_db_path,
        perc_identity=identity,
        word_size=word_size
    )
    # Keep only essential columns
    initial_blast_results = initial_blast_results[['qseqid', 'sseqid', 'sstart', 'send', 'sstrand']]
    toc: float = time.perf_counter()
    print(f"1. Initial data:\n"
          f"\t- Data row length: {len(initial_blast_results)}\n"
          f"\t- Execution time: {toc - tic:0.2f} seconds")

    # Merge results using PyRanges
    print('\t- Filtering and merging data:')
    tic = time.perf_counter()
    merged_results: pd.DataFrame = merge_overlapping_intervals(initial_blast_results, strand=True)
    merged_results.sort_values(by=['sseqid', 'sstart'], inplace=True)
    merged_results.reset_index(drop=True, inplace=True)
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
        n_jobs: int
) -> pd.DataFrame:
    """
    Finds and extends repetitive sider regions from a genomic dataset iteratively.

    This function leverages BLAST alignment and recursive extension to identify and iteratively extend repetitive sider
    regions within a given genomic dataset. The process involves extending sequences, performing BLASTn-based searches
    to find new regions of interest, and filtering overlapping sequences to accumulate a final set of unique, extended
    repetitive sider regions.

    Parameters
    ----------
    data_input : pd.DataFrame
        Initial DataFrame containing sequences to start the extension and search process.
    genome_path : str
        Path to the reference genome in FASTA format used for sequence extension and BLASTn alignment.
    extend_number : int
        The number of base pairs to extend during each recursive step.
    word_size : int
        The BLASTn parameter indicating seed word size for alignment.
    identity : int
        Minimum percentage identity required for BLAST alignment matches.
    min_length : int
        Minimum length of sequence matches to retain after filtering.
    limit_len : int
        Maximum allowable length for extended sequences during processing.
    n_jobs : int
        Number of parallel jobs/threads to use for computational steps such as BLAST alignment.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the final, accumulated set of extended repetitive sider regions, with overlaps resolved
        and filtered sequences meeting the defined criteria.
    """

    # Extend data using the recursive method
    # -------------------------------------------------------------------
    # STEP 1: Extend the input data using the recursive method
    # -------------------------------------------------------------------
    print(f"Extending initial {len(data_input)} sequences")
    logger.info(f"Extending initial {len(data_input)} sequences")
    data_extended: pd.DataFrame = sequence_extension(
        data_input=data_input,
        genome_fasta=genome_path,
        extend_number=extend_number,
        limit_len=limit_len,
        identity=identity,
        word_size=word_size,
        min_length=min_length,
        n_jobs=n_jobs
    )

    # -------------------------------------------------------------------
    # STEP 2: Launch second Blast to find new elements
    # -------------------------------------------------------------------
    accumulated_data: pd.DataFrame = data_extended.copy()
    iteration: int = 0
    are_there_new_elems: bool = True
    while are_there_new_elems:
        iteration += 1
        print_message_box(f"Iteration {iteration}:")
        new_elems: pd.DataFrame = pd.DataFrame()  # Will hold the newly discovered elements. Will be used to reject overlaps
        i: int; row: pd.Series
        for i, row in data_extended.iterrows():
            print(f"\n-Analyzing row {i+1}/{len(data_extended)}:")
            # -------------------------------------------------------------------
            # STEP 2.1: Perform the BLASTn to find new elements
            # -------------------------------------------------------------------
            # Create a temporary FASTA file with the sequence and ensure cleanup
            temp_fasta_path: Optional[str] = None
            try:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
                    seq_stripped: str = row.sseq.replace('-', '')  # Remove hyphens
                    temp_fasta.write(f">{i}_{row.sseqid}_{row.sstart}_{row.send}\n{seq_stripped}\n")
                    temp_fasta_path = temp_fasta.name

                # Perform the BLAST search
                blast_results: pd.DataFrame = run_blastn_alignment(
                    query_path=temp_fasta_path,
                    dict_path=genome_path,
                    perc_identity=identity,
                    word_size=word_size
                )
                if not blast_results.empty:
                    blast_results = blast_results[['sseqid', 'sstart', 'send', 'sstrand', 'len', 'sseq']]
                blast_results.sort_values(by=['sseqid', 'sstart'], inplace=True)
            finally:
                if temp_fasta_path and os.path.exists(temp_fasta_path):
                    os.unlink(temp_fasta_path)

            # -------------------------------------------------------------------
            # STEP 2.2: Filter the data
            # -------------------------------------------------------------------
            # Remove elements with a len < min_length
            if not blast_results.empty:
                blast_results = blast_results[blast_results['len'] >= min_length]

            # Remove sequences that overlap with `accumulated_data` or `new_elems`
            ## With `data_extended`
            rm_from_accumulated_data: pd.DataFrame = fetch_overlapping_intervals(blast_results, accumulated_data)
            if not blast_results.empty and not rm_from_accumulated_data.empty:
                print(f"\t- Overlaps with accumulated_data removed: {len(rm_from_accumulated_data)}/{len(blast_results)}")
                blast_results = match_data_and_remove(blast_results, rm_from_accumulated_data)

            ## With `new_elems`
            if not new_elems.empty:
                rm_from_new_elems: pd.DataFrame = fetch_overlapping_intervals(blast_results, new_elems)
                if not rm_from_new_elems.empty:
                    print(f"\t- Overlaps with new_elems removed: {len(rm_from_new_elems)}/{len(blast_results)}")
                    blast_results = match_data_and_remove(blast_results, rm_from_new_elems)

            # -------------------------------------------------------------------
            # STEP 2.3: Add data to the collection of `new_elems`
            # -------------------------------------------------------------------
            if not blast_results.empty:
                new_elems = pd.concat([new_elems, blast_results], ignore_index=True)
                new_elems.sort_values(by=['sseqid', 'sstart'], inplace=True)

        if new_elems.empty:
            # -------------------------------------------------------------------
            # STEP 3a: No new elements found
            # -------------------------------------------------------------------
            # If nothing was found, exit while loop
            are_there_new_elems = False
        else:
            # -------------------------------------------------------------------
            # STEP 3b: New elements found
            # -------------------------------------------------------------------
            print('='*40)
            print(f"\n- Extending {len(new_elems)} new sequences")
            new_elems_extended: pd.DataFrame = sequence_extension(
                data_input=new_elems,
                genome_fasta=genome_path,
                extend_number=extend_number,
                limit_len=limit_len,
                identity=identity,
                word_size=word_size,
                min_length=min_length,
                n_jobs=n_jobs,
                prune_enabled=True,
                prune_against_df=accumulated_data
            )

            # Check overlapping data of `new_elems_extended` with `accumulated_data`
            overlapping_elems_extended: pd.DataFrame = fetch_overlapping_intervals(new_elems_extended, accumulated_data)
            if not overlapping_elems_extended.empty:
                print(f"\t- Remove overlaps with accumulated_data: {len(overlapping_elems_extended)}")
                new_elems_extended = match_data_and_remove(new_elems_extended, overlapping_elems_extended)
                new_elems_extended.sort_values(by=['sseqid', 'sstart'], inplace=True)
                new_elems_extended.reset_index(drop=True, inplace=True)

            accumulated_data = pd.concat([accumulated_data, new_elems_extended], ignore_index=True)
            accumulated_data.sort_values(by=['sseqid', 'sstart'], inplace=True)
            accumulated_data.reset_index(drop=True, inplace=True)

            # Replace `data_extended` that entered the while loop with the new elements discovered
            data_extended = new_elems_extended.copy()
            data_extended.sort_values(by=['sseqid', 'sstart'], inplace=True)
            data_extended.reset_index(drop=True, inplace=True)

    # -------------------------------------------------------------------
    # STEP 4: Outside while-loop
    # -------------------------------------------------------------------
    return accumulated_data


def finalize_results(
        output_dir: str,
        df: pd.DataFrame,
        input_file: str,
        genome_file: str
) -> Tuple[str, str]:
    """
    Finalize and save the analysis results to specific file formats.

    This function finalizes the analysis results by formatting the supplied
    DataFrame and writing it to both CSV and GFF file formats. The output
    filenames are determined using the given input and genome filenames, with
    the target directory specified as `output_dir`. The formatted results are
    written to these files, and their file paths are returned for further
    processing or validation.

    Parameters
    ----------
    output_dir : str
        The directory where the results will be saved.
    df : pandas.DataFrame
        The DataFrame containing the analysis results to be formatted and saved.
    input_file : str
        The file path of the input dataset used in the analysis.
    genome_file : str
        The file path of the genome data used in the analysis.

    Returns
    -------
    Tuple[str, str]
        A tuple containing the file paths of the saved CSV and GFF files. The
        first element is the path to the CSV file, and the second element is the
        path to the GFF file.

    """
    # Build output filenames using basenames
    input_base: str = os.path.splitext(os.path.basename(input_file))[0]
    genome_base: str = os.path.splitext(os.path.basename(genome_file))[0]
    csv_path: str = os.path.join(output_dir, f"BLASTOISE--{input_base}--{genome_base}.csv")
    gff_path: str = os.path.join(output_dir, f"BLASTOISE--{input_base}--{genome_base}.gff")

    # Format DataFrame for CSV
    formatted: pd.DataFrame = format_output_dataframe(df)

    # Write CSV (with header)
    formatted.to_csv(csv_path, index=False)

    # Write GFF from formatted
    write_gff_from_formatted(formatted, gff_path)

    return csv_path, gff_path


def main() -> None:
    """
    Main entry point of the program that orchestrates multiple steps for processing genomic and data files.
    The workflow includes parsing arguments, setting up the workspace, running initial and iterative BLAST searches,
    processing results, and performing final cleanup and output. It also handles timing, error reporting, and logging.

    Raises
    ------
    SystemExit
        Raised when there is an unhandled exception during the program execution, causing an exit with error message.

    Notes
    -----
    The program measures the total execution time, logs the program start and end times, and summarizes the
    results paths for CSV and GFF files. Proper file permissions are attempted for the output directory at
    the end of the program.
    """
    # Configure root logger to avoid console output; only attach file handler later
    root_logger: logging.Logger = logging.getLogger()
    # Remove any pre-existing handlers (e.g., from environments that pre-configure logging)
    root_logger.handlers.clear()
    root_logger.setLevel(logging.INFO)
    # Add a NullHandler to suppress lastResort console output until file handler is attached
    if not any(isinstance(h, logging.NullHandler) for h in root_logger.handlers):
        root_logger.addHandler(logging.NullHandler())

    args: argparse.Namespace = parse_arguments()
    start_time: datetime = datetime.now()
    tic_main: float = time.perf_counter()
    formatted_start_time: str = start_time.strftime('%Y %B %d at %H:%M')
    print(f"- Program started: {formatted_start_time}")

    try:
        # 1. Setup workspace and copy input files
        output_dir: str; data_path: str; genome_path: str
        output_dir, data_path, genome_path = setup_workspace(
            output_dir=os.path.expanduser(args.output),
            data_file=os.path.expanduser(args.data),
            genome_file=os.path.expanduser(args.genome)
        )
        # Add a file handler to log to the output directory
        try:
            log_file_path: str = os.path.join(output_dir, 'blastoise.log')
            file_handler: logging.FileHandler = logging.FileHandler(log_file_path)
            file_handler.setLevel(logging.INFO)
            file_handler.setFormatter(
                logging.Formatter('%(asctime)s [%(levelname)s] %(name)s:%(funcName)s:%(lineno)d: %(message)s')
            )
            root_logger: logging.Logger = logging.getLogger()
            # Check for existing FileHandler with same base filename
            has_file_handler: bool = any(
                isinstance(h, logging.FileHandler) and
                getattr(h, 'baseFilename', None) == file_handler.baseFilename
                for h in root_logger.handlers
            )

            # Add handler if not already present
            if not has_file_handler:
                root_logger.addHandler(file_handler)
            logger.info("Logging to file: %s", log_file_path)
            # Now that the file handler is attached, record the start message to the log file
            logger.info("Program started: %s", formatted_start_time)
        except Exception as log_e:
            logger.warning("Failed to attach file handler for logging: %s", log_e)

        logger.info(
            "Workspace ready:\n\tData: %s\n\tGenome: %s\n\tOutput Path: %s", data_path, genome_path, output_dir
        )

        # 2. Run initial BLAST and process results
        initial_data: pd.DataFrame; blast_db_path: str
        initial_data, blast_db_path = run_initial_blast(
            data_path=data_path,
            genome_path=genome_path,
            output_dir=output_dir,
            identity=args.identity,
            word_size=args.word_size
        )
        logger.info(f"Initial BLAST completed. \n\tRows: {len(initial_data)} \n\tDB: {blast_db_path}")

        # Early exit if nothing to process
        csv_path: str; gff_path: str
        if initial_data.empty:
            print_message_box("No data to process after initial BLAST. Writing empty results.")
            logger.info("No data after initial BLAST. Writing empty results.")
            csv_path, gff_path = finalize_results(output_dir, initial_data, args.data, args.genome)
        else:
            # 3. Run the iterative part
            print_message_box("Running iterative search")
            logger.info(f"Running iterative search for {len(initial_data)} sequences.")
            final_data: pd.DataFrame = repetitive_sider_searcher(
                data_input=initial_data,
                genome_path=blast_db_path,
                extend_number=args.extend,
                word_size=args.word_size,
                identity=args.identity,
                min_length=args.min_length,
                limit_len=args.limit,
                n_jobs=args.jobs,
            )

            # 4. Finalize and clean up
            csv_path, gff_path = finalize_results(output_dir, final_data, args.data, args.genome)
            logger.info(f"Final results written. \n\tCSV: {csv_path} \n\tGFF: {gff_path}")

    except Exception as e:
        logger.exception("Unhandled exception during execution")
        print_message_box(f"An unexpected error occurred: {e}")
        exit(1)

    # 5. Print final summary
    toc_main: float = time.perf_counter()
    end_time: datetime = datetime.now()
    formatted_end_time: str = end_time.strftime("%Y %B %d at %H:%M")

    print_message_box(message="END OF THE PROGRAM")
    print(f"\t- Total execution time: {toc_main - tic_main:0.2f} seconds\n"
          f"\t- Program started: {formatted_start_time}\n"
          f"\t- Program ended: {formatted_end_time}\n"
          f"\t- CSV saved at: {csv_path if 'csv_path' in locals() else 'N/A'}\n"
          f"\t- GFF saved at: {gff_path if 'gff_path' in locals() else 'N/A'}")

    logger.info(f"Program finished. \n\tDuration: {toc_main - tic_main:.2f} seconds | \n\tEnd: {formatted_end_time}")

    # Set file permissions for the output directory
    try:
        subprocess.run(["chmod", "-R", "a+w", output_dir], check=True)
    except subprocess.CalledProcessError as e:
        logger.warning(f"Could not set final permissions on {output_dir}: {e}")
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
