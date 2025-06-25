"""
Coordinate Corrector Functions Module

This module contains utility functions for coordinate correction operations.
These functions support the main coordinate correction pipeline by providing functionality for:
- Command-line argument parsing
- Directory and database setup
- Sequence coordinate correction using BLASTN
- Filtering BLASTN results
- Creating output DataFrames
- Adding sequences to DataFrames
- Saving results to disk

This module is designed to be imported by the main coordinate_corrector.py script.
"""

import argparse
import os
import sys
import logging
import pandas as pd
import pyranges as pr
import subprocess
from typing import Dict, List, Any

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from .extra_functions import fetch_dna_sequence, general_blastn_blaster
from blastoise.modules.aesthetics import print_message_box


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the coordinate corrector script.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments with the following attributes:
        - file: Path to the input CSV file
        - dict_path: Path to the genome FASTA file
        - word_size: Word size parameter for BLASTN
        - min_length: Minimum length of the sequence to be considered
        - identity: Percentage identity threshold for BLASTN
        - evalue: E-value threshold for BLASTN
    """
    parser = argparse.ArgumentParser(
        description="Correct coordinates for sequences in a CSV file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-f", "--file", 
        type=str, 
        required=True,
        help="Path to the input CSV file containing sequence data."
    )

    parser.add_argument(
        "-d", "--dict_path", 
        type=str, 
        required=True,
        help="Path to the genome FASTA file for BLASTN database creation."
    )

    parser.add_argument(
        "-ws", "--word_size", 
        type=int, 
        default=15,
        help="Word size parameter for BLASTN."
    )

    parser.add_argument(
        "-min", "--min_length", 
        type=int, 
        default=100,
        help="Minimum length of the sequence to be considered after correction."
    )

    parser.add_argument(
        "-max", "--max_length",
        type=int,
        default=1000,
        help="Maximum length of the sequence to be considered after correction."
    )

    parser.add_argument(
        "-i", "--identity",
        type=int,
        default=60,
        help="Percentage identity threshold for BLASTN search."
    )

    parser.add_argument(
        "-e", "--evalue",
        type=float,
        default=1.0E-09,
        help="E-value threshold for BLASTN search."
    )

    return parser.parse_args()


def merge_intervals_qrange(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping or adjacent genomic intervals within a DataFrame using qstart and qend columns.

    This function takes a DataFrame containing genomic interval information with 'qstart'
    and 'qend' columns, converts it to PyRanges format for merging, and returns the
    processed intervals with the original qstart/qend column structure preserved.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing genomic intervals to be merged. Must contain columns
        'sseqid', 'qstart', and 'qend'.

    Returns
    -------
    pd.DataFrame
        A DataFrame with merged genomic intervals, maintaining the original
        structure with columns 'sseqid', 'qstart', and 'qend'. If the input DataFrame
        is empty or no intervals can be merged, an empty DataFrame with these columns
        is returned.

    Raises
    ------
    KeyError
        If any of the required columns are missing from the input DataFrame.
    ValueError
        If the input DataFrame cannot be converted to a PyRanges object.
    """
    if df.empty:
        return pd.DataFrame(columns=['sseqid', 'qstart', 'qend'])

    # Create a mapping for qstart/qend to PyRanges format
    df_for_pyranges = df.rename(columns={
        "qstart": "Start",
        "qend": "End"
    })

    # Add a dummy chromosome column since PyRanges requires it
    df_for_pyranges["Chromosome"] = "dummy"

    # Convert to PyRanges, merge, and convert back
    pr_df = pr.PyRanges(df_for_pyranges)
    merged = pr_df.merge()

    # Handle an empty result
    if hasattr(merged, 'df'):
        # Convert back to original column names
        result_df = merged.df.rename(columns={
            "Start": "qstart",
            "End": "qend"
        })
        return result_df
    else:
        return pd.DataFrame(columns=['sseqid', 'qstart', 'qend'])


def correct_coordinates(
    df: pd.DataFrame, 
    blastn_db_path: str, 
    word_size: int,
    min_length: int,
    max_length: int,
    identity: int,
    evalue: float
) -> Dict[str, List[List[Any]]]:
    """
    Adjust sequence coordinates using BLASTN and PyRanges processing.

    This function takes a DataFrame containing sequence information, uses BLASTN to process sequences,
    filters sequences based on coordinate overlap and strand validity, merges the results with PyRanges,
    and computes updated start and end coordinates.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing sequence data with columns for sseqid, sstrand, sstart, send, and sseq.
    blastn_db_path : str
        Path to the BLASTN database used for sequence alignment.
    word_size : int
        Word size parameter for BLASTN search (-W parameter).
    min_length : int
        Minimum length of the sequence to be considered after correction.
    max_length: int
        Maximum length of the sequence to be considered after correction.
    identity : int, optional
        Percentage identity threshold for BLASTN search. Default is 60.
    evalue : float, optional
        E-value threshold for BLASTN search. The default is 1.0E-09.

    Returns
    -------
    Dict[str, List[List[Any]]]
        Dictionary where keys are sequence identifiers and values are lists of lists,
        each inner list containing [chromosome, strand, start, end] for a corrected sequence.
    """
    logger = logging.getLogger('coordinate_corrector')
    logger.info("Starting coordinate correction process")

    # Create a dictionary to store results
    results_dict = {}

    # Process each row in the DataFrame
    for pos, (index, row) in enumerate(df.iterrows(), start=1):
        # Prepare the query data
        name_id = f"{row.sseqid}_{row.sstrand}_{row.sstart}-{row.send}"
        seq = row.sseq
        query = f"<(echo -e '>{name_id}\\n{seq}')"
        start_coor = row.sstart
        end_coor = row.send
        strand_seq = row.sstrand
        name_chr = row.sseqid
        seq_len = end_coor - start_coor + 1

        logger.info(f"Analyzing row {pos}/{df.shape[0]} with name_id {name_id}")

        # Run BLASTN
        blastn_df = general_blastn_blaster(
            query_path=query,
            dict_path=blastn_db_path,
            word_size=word_size,
            perc_identity=identity,
            evalue=evalue
        )

        # Filter the BLASTN results
        # Remove rows with coordinates that overlap with the original coordinates
        blastn_df: pd.DataFrame = filter_blastn_results(
            blastn_df,
            start_coor,
            end_coor,
            name_chr
        )

        # Add len column for the subject
        blastn_df['len'] = abs(blastn_df['qend'] - blastn_df['qstart']) + 1

        # Only when the length is >= min_length
        blastn_df = blastn_df.loc[blastn_df['len'] >= min_length]
        blastn_df = blastn_df.loc[blastn_df['len'] <= max_length]

        # Sort and merge the filtered results
        blastn_df.sort_values('qstart', inplace=True)

        # Merge intervals using PyRanges
        merged_df = merge_intervals_qrange(blastn_df)

        # Initialize the results for this sequence
        results_dict[name_id] = []

        # Calculate new coordinates
        for _, merged_row in merged_df.iterrows():
            # Calculate new start and end coordinates
            new_start = start_coor + merged_row["qstart"] - 1
            new_end = start_coor + merged_row["qend"] - 1

            # Check if the sequence meets the minimum length requirement
            if abs(new_end - new_start) + 1 >= min_length:
                results_dict[name_id].append([name_chr, strand_seq, new_start, new_end])

        # Log the results
        print(f"Found {len(results_dict[name_id])} sequences for {name_id} ==> len {end_coor - start_coor + 1}")
        for seq in results_dict[name_id]:
            # Check if the coordinates changed:
            if seq[2] != start_coor or seq[3] != end_coor:
                print(f"\t{seq[0]}:{seq[1]}:{seq[2]}-{seq[3]} ==> len {seq[3] - seq[2] + 1}")
            else:
                print("\tNo changes made")

    return results_dict


def filter_blastn_results(
    blastn_df: pd.DataFrame, 
    start_coor: int, 
    end_coor: int, 
    name_chr: str
) -> pd.DataFrame:
    """
    Filter BLASTN results to remove rows with coordinates that overlap with the original coordinates.

    Parameters
    ----------
    blastn_df : pd.DataFrame
        A DataFrame containing BLASTN results with columns for sstart, send, sseqid, and sstrand.
    start_coor : int
        Start coordinate of the original sequence.
    end_coor : int
        End coordinate of the original sequence.
    name_chr : str
        Chromosome name of the original sequence.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with overlapping regions removed.
    """
    # Filter for plus strand
    # Remove rows where:
    # (sstart is within start_coor and end_coor, OR send is within start_coor and end_coor),
    # AND sseqid matches name_chr, AND sstrand is "plus"
    plus_filter = (
        (((blastn_df["sstart"] >= start_coor) & (blastn_df["sstart"] <= end_coor)) |
         ((blastn_df["send"] <= end_coor) & (blastn_df["send"] >= start_coor))) &
        (blastn_df["sseqid"] == name_chr) &
        (blastn_df["sstrand"] == "plus")
    )
    blastn_df = blastn_df[~plus_filter].copy()

    # Filter for minus strand
    # Remove rows where:
    # (sstart is within start_coor and end_coor, OR send is within start_coor and end_coor),
    # AND sseqid matches name_chr, AND sstrand is "minus"
    minus_filter = (
        (((blastn_df["sstart"] <= end_coor) & (blastn_df["sstart"] >= start_coor)) |
         ((blastn_df["send"] >= start_coor) & (blastn_df["send"] <= end_coor))) &
        (blastn_df["sseqid"] == name_chr) &
        (blastn_df["sstrand"] == "minus")
    )
    blastn_df = blastn_df[~minus_filter].copy()

    return blastn_df


def create_output_dataframe(results_dict: Dict[str, List[List[Any]]]) -> pd.DataFrame:
    """
    Create a DataFrame from the dictionary of corrected coordinates.

    Parameters
    ----------
    results_dict : Dict[str, List[List[Any]]]
        Dictionary where keys are sequence identifiers and values are lists of lists,
        each inner list containing [chromosome, strand, start, end] for a corrected sequence.

    Returns
    -------
    pd.DataFrame
       A DataFrame containing corrected coordinates with columns for sseqid, sstrand, sstart, and send.
    """
    logger = logging.getLogger('coordinate_corrector')
    logger.info("Creating output DataFrame")

    # Extract data from the dictionary and create a list of rows
    rows = []
    for key, value in results_dict.items():
        for element in value:
            rows.append({
                'sseqid': element[0],
                'sstrand': element[1],
                'sstart': element[2],
                'send': element[3]
            })

    # Create DataFrame
    result_df = pd.DataFrame(rows)

    # Reorder
    result_df.sort_values(['sseqid', 'sstart'], inplace=True)

    logger.info(f"Created DataFrame with {len(result_df)} rows")

    return result_df


def add_sequences(df: pd.DataFrame, blastn_db_path: str) -> pd.DataFrame:
    """
    Add sequence data to the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing coordinate data with columns for sseqid, sstrand, sstart, and send.
    blastn_db_path : str
        Path to the BLASTN database used to retrieve sequences.

    Returns
    -------
    pd.DataFrame
        A DataFrame with added 'sseq' column containing the nucleotide sequences.
    """
    logger = logging.getLogger('coordinate_corrector')
    if df.empty:
        logger.warning("DataFrame is empty, no sequences to add")
        return df

    logger.info("Adding sequences to DataFrame")

    # Add sequence data to each row
    df['sseq'] = df.apply(
        lambda x: fetch_dna_sequence(
            start_coor=x['sstart'],
            end_coor=x['send'],
            strand=x['sstrand'],
            chromosome=x['sseqid'],
            path_genome=blastn_db_path
        ), 
        axis=1
    )

    logger.info("Sequences added to DataFrame")

    return df


def save_results(df: pd.DataFrame, output_path: str, temp_dir: str) -> None:
    """
    Save the DataFrame to a CSV file and set permissions.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing the corrected coordinates and sequences to save.
    output_path : str
        Path to the output CSV file where results will be saved.
    temp_dir : str
        Path to the temporary directory that needs permission changes.

    Returns
    -------
    None
    """
    logger = logging.getLogger('coordinate_corrector')
    if df.empty:
        logger.warning("DataFrame is empty, no results to save")
        print_message_box("No sequences found after coordinate correction.")
        return

    logger.info(f"Saving results to {output_path}")

    # Save DataFrame to CSV
    df.to_csv(output_path, index=False)

    # Set file permissions
    try:
        subprocess.run(["chmod", "-R", "a+w", temp_dir], check=True)
        subprocess.run(["chmod", "a+w", output_path], check=True)
        logger.info("File permissions set")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error setting file permissions: {str(e)}")
        print(f"Error occurred while changing permissions: {e}")

    print_message_box(f"Coordinate correction completed. Results saved to {output_path}")
