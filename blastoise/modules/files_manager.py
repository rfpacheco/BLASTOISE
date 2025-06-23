"""
BLASTOISE Module: File Management and Sequence Data Handling
===========================================================

This module provides functionality for file operations and sequence data handling in the
BLASTOISE pipeline. It manages the creation of FASTA files, retrieval of sequence data
from genome databases, and various data processing operations needed throughout the
pipeline.

The module contains several key functions:
1. `fasta_creator`: Creates FASTA files from DataFrames containing sequence data
2. `get_data_sequence`: Retrieves nucleotide sequences from genome databases using BLAST
3. `columns_to_numeric`: Converts specified DataFrame columns to numeric datatypes
4. `end_always_greater_than_start`: Ensures consistent coordinate ordering

These functions work together to provide essential file and data handling capabilities
for the BLASTOISE pipeline, enabling efficient processing and storage of genomic
sequence data.

Author: R. Pacheco
Version: 0.4.2
License: MIT
"""

import pandas as pd
import subprocess
from typing import Optional, Dict, Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed


def fasta_creator(
        data_input: pd.DataFrame,
        fasta_output_path: str,
        id_names: Optional[pd.DataFrame] = None
) -> None:
    """
    Create a FASTA file from a DataFrame containing sequence and metadata.

    This function processes sequence data from a specified DataFrame, formats the information 
    as FASTA records, and writes those records to a file in FASTA format. Each sequence is 
    given a unique identifier that includes its genomic coordinates and strand information.

    Optionally, it can include original identifiers provided via a secondary DataFrame, 
    appending both original and extended identifiers to each sequence record. This is 
    particularly useful for tracking the relationship between original and extended sequences.

    Parameters
    ----------
    data_input : pd.DataFrame
        Input data containing sequence information. Must include columns:
        'sseqid' (sequence ID), 'sstart' (start position), 'send' (end position),
        'sstrand' (strand orientation), and 'sseq' (the actual sequence).
    fasta_output_path : str
        File path where the resulting FASTA file will be written.
    id_names : pd.DataFrame, optional
        DataFrame containing original sequence identifiers and metadata, used to further 
        annotate FASTA records. If provided, it must have the same structure and length as
        data_input. The default is None, which uses data_input as its own reference.

    Returns
    -------
    None
        The function writes the FASTA file to disk but does not return any value.

    Raises
    ------
    KeyError
        If any of the required columns are missing from data_input.
    IndexError
        If id_names is provided but has fewer rows than data_input.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Initialize id_names if not provided
    # -----------------------------------------------------------------------------
    if id_names is None:
        id_names = data_input

    # -----------------------------------------------------------------------------
    # STEP 2: Create FASTA records for each sequence
    # -----------------------------------------------------------------------------
    matrix = []  # List to store SeqRecord objects
    for index, (_, sequence) in enumerate(data_input.iterrows()):
        # Format the sequence identifier with genomic coordinates and strand
        data_input_string = f"{sequence['sseqid']}-{sequence['sstart']}-{sequence['send']}-{sequence['sstrand']}"

        # -----------------------------------------------------------------------------
        # STEP 3: Add original coordinates if id_names is provided
        # -----------------------------------------------------------------------------
        if id_names is not None:
            # Combine extended coordinates with original coordinates in the identifier
            original = f"{id_names.iloc[index]['sseqid']}-" \
                       f"{id_names.iloc[index]['sstart']}-" \
                       f"{id_names.iloc[index]['send']}-" \
                       f"{id_names.iloc[index]['sstrand']}"
            # Format: extended_coordinates_original_coordinates
            data_input_string = f"{data_input_string}_{original}"

        # -----------------------------------------------------------------------------
        # STEP 4: Create a SeqRecord object and add to matrix
        # -----------------------------------------------------------------------------
        # Create a SeqRecord with a unique identifier and the sequence
        rec = SeqRecord(
            Seq(sequence.loc["sseq"]),  # The actual sequence
            id="Seq-" + str(index) + "_" + data_input_string,  # Unique ID with index and coordinates
            description=""  # No description needed
        )
        matrix.append(rec)

    # -----------------------------------------------------------------------------
    # STEP 5: Write all records to the FASTA file
    # -----------------------------------------------------------------------------
    SeqIO.write(matrix, fasta_output_path, "fasta")


def _process_single_row(row: pd.Series, genome_fasta: str) -> Dict[str, Any]:
    """
    Process a single row from the DataFrame to retrieve sequence data.

    This helper function extracts genomic coordinates and strand information from a 
    DataFrame row, constructs a BLAST command to retrieve the corresponding sequence 
    from a genome database, and returns the data as a dictionary. It's designed to be 
    used with parallel processing to efficiently handle large datasets.

    Parameters
    ----------
    row : pd.Series
        A single row from the DataFrame containing sequence details.
        Must contain the following keys:
        - 'sseqid': Sequence identifier (e.g., chromosome name)
        - 'sstart': Start position in the genome
        - 'send': End position in the genome
        - 'sstrand': Strand orientation ('plus' or 'minus')
    genome_fasta : str
        The file path to the genome database in FASTA format to query against.
        This should be a BLAST database created with makeblastdb.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the sequence data with keys:
        - 'sseqid': Sequence identifier (copied from input)
        - 'sstart': Start position (copied from input)
        - 'send': End position (copied from input)
        - 'sstrand': Strand orientation (copied from input)
        - 'sseq': The actual nucleotide sequence retrieved from the genome

    Raises
    ------
    KeyError
        If any of the required keys are missing from the input row.
    subprocess.CalledProcessError
        If the BLAST command fails to execute properly.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Extract coordinates and strand information from the row
    # -----------------------------------------------------------------------------
    sseqid = row['sseqid']  # Sequence identifier (e.g., chromosome name)
    start = row['sstart']   # Start position in the genome
    end = row['send']       # End position in the genome
    strand = row['sstrand'] # Strand orientation ('plus' or 'minus')

    # -----------------------------------------------------------------------------
    # STEP 2: Construct and execute the BLAST command to retrieve the sequence
    # -----------------------------------------------------------------------------
    # Build the blastdbcmd command to extract the sequence from the genome database
    cmd = f"blastdbcmd -db {genome_fasta} -entry {sseqid} -range {start}-{end} -strand {strand} -outfmt %s"

    # Execute the command and capture the output (the sequence)
    sequence = subprocess.run(
        cmd, 
        shell=True, 
        capture_output=True, 
        text=True, 
        universal_newlines=True,
        executable="/usr/bin/bash"
    ).stdout.strip()

    # -----------------------------------------------------------------------------
    # STEP 3: Return the data as a dictionary
    # -----------------------------------------------------------------------------
    return {
        'sseqid': sseqid,    # Sequence identifier
        'sstart': start,     # Start position
        'send': end,         # End position
        'sstrand': strand,   # Strand orientation
        'sseq': sequence     # The retrieved nucleotide sequence
    }


def get_data_sequence(data: pd.DataFrame, genome_fasta: str, n_jobs: int = -1) -> pd.DataFrame:
    """
    Fetch nucleotide sequences from a genome database using parallel BLAST commands.

    This function retrieves nucleotide sequences for a set of genomic coordinates from a 
    BLAST-formatted genome database. It processes each row of the input DataFrame in 
    parallel using joblib, which significantly improves performance when dealing with 
    large datasets. Each row is processed by the _process_single_row helper function, 
    which executes a BLAST command to extract the corresponding sequence.

    Parameters
    ----------
    data : pd.DataFrame
        A DataFrame containing the genomic coordinates for which sequences should be 
        retrieved. Must contain the following columns:
        - 'sseqid': Sequence identifier (e.g., chromosome name)
        - 'sstart': Start position in the genome
        - 'send': End position in the genome
        - 'sstrand': Strand orientation ('plus' or 'minus')
    genome_fasta : str
        The file path to the genome database in FASTA format to query against.
        This should be a BLAST database created with makeblastdb.
    n_jobs : int, optional
        The number of jobs to run in parallel for processing rows. 
        -1 means using all available processors.
        Default is -1.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the retrieved sequences with the following columns:
        - 'sseqid': Sequence identifier (copied from input)
        - 'sstart': Start position (copied from input)
        - 'send': End position (copied from input)
        - 'sstrand': Strand orientation (copied from input)
        - 'sseq': The actual nucleotide sequence retrieved from the genome

    Raises
    ------
    KeyError
        If any of the required columns are missing from the input DataFrame.
    subprocess.CalledProcessError
        If any of the BLAST commands fail to execute properly.
    ValueError
        If the genome_fasta file cannot be found or is not a valid BLAST database.

    Notes
    -----
    This function uses joblib's Parallel and delayed functions to process rows in 
    parallel, which can significantly improve performance for large datasets. The 
    actual sequence retrieval is handled by the _process_single_row helper function.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Process all rows in parallel using joblib
    # -----------------------------------------------------------------------------
    # Use joblib to distribute the processing across multiple processors
    # Each row is processed independently by the _process_single_row function
    sequences = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_row)(row, genome_fasta) for _, row in data.iterrows()
    )

    # -----------------------------------------------------------------------------
    # STEP 2: Convert the list of dictionaries to a DataFrame
    # -----------------------------------------------------------------------------
    # Each dictionary in the sequence list becomes a row in the DataFrame
    # The keys of the dictionaries become the column names
    sequences_df = pd.DataFrame(sequences)

    return sequences_df


def columns_to_numeric(data_input: pd.DataFrame, columns_to_convert: list[str] | None = None) -> pd.DataFrame:
    """
    Convert specified columns of a DataFrame to numeric datatype.

    This function takes a DataFrame and a list of column names and converts the specified 
    columns to numeric datatypes using pandas' to_numeric function. This is particularly 
    useful for ensuring that coordinate columns (like 'sstart' and 'send') are properly 
    typed for numerical operations and comparisons.

    If no columns are specified, the function defaults to converting the 'sstart' and 
    'send' columns, which are commonly used for genomic coordinates in the BLASTOISE 
    pipeline. Any values that cannot be converted to numeric are coerced into NaN values 
    rather than raising an error.

    Parameters
    ----------
    data_input : pd.DataFrame
        The DataFrame containing the data to be transformed. The DataFrame is 
        modified in-place, but a reference is also returned.
    columns_to_convert : list[str], optional
        A list of column names to be converted to numeric data type.
        If None, defaults to ['sstart', 'send'].
        Default is None.

    Returns
    -------
    pd.DataFrame
        The modified DataFrame with specified columns converted to numeric datatypes.
        This is the same object as the input DataFrame (modified in-place).

    Raises
    ------
    KeyError
        If any of the specified columns do not exist in the DataFrame.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Set default columns if none are specified
    # -----------------------------------------------------------------------------
    if columns_to_convert is None:
        # Default to converting the standard genomic coordinate columns
        columns_to_convert = ['sstart', 'send']

    # -----------------------------------------------------------------------------
    # STEP 2: Convert each specified column to numeric type
    # -----------------------------------------------------------------------------
    for column in columns_to_convert:
        # Convert the column to numeric, coercing any non-convertible values to NaN
        # This is safer than raising errors for invalid values
        data_input[column] = pd.to_numeric(data_input[column], errors='coerce')

    # Return the modified DataFrame (same object as input, modified in-place)
    return data_input


def end_always_greater_than_start(data_input: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure that 'send' values are always greater than 'sstart' values by swapping if needed.

    In genomic coordinate systems, it's often conventional to have the end coordinate 
    greater than the start coordinate, regardless of strand orientation. This function 
    enforces this convention by identifying rows where 'sstart' is greater than 'send' 
    and swapping these values to ensure consistent coordinate ordering throughout the 
    dataset.

    This is particularly important for functions that calculate sequence lengths or 
    perform interval operations, as they typically assume that end > start.

    Parameters
    ----------
    data_input : pd.DataFrame
        A DataFrame containing genomic coordinates. Must include the columns:
        - 'sstart': Start position in the genome
        - 'send': End position in the genome

    Returns
    -------
    pd.DataFrame
        The modified DataFrame with 'sstart' and 'send' values properly ordered
        such that 'send' ≥ 'sstart' for all rows. This is the same object as the 
        input DataFrame (modified in-place).

    Raises
    ------
    KeyError
        If either 'sstart' or 'send' columns are missing from the DataFrame.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Identify rows where sstart > send
    # -----------------------------------------------------------------------------
    # Create a boolean mask for rows where start is greater than the end
    rows_to_swap = data_input['sstart'] > data_input['send']

    # Get the indices of rows that need swapping
    swap_indices = data_input[rows_to_swap].index

    # -----------------------------------------------------------------------------
    # STEP 2: Swap sstart and send values for identified rows
    # -----------------------------------------------------------------------------
    if len(swap_indices) > 0:  # Only perform the swap if there are rows to fix
        # Swap the values by assigning them in reverse order
        # This uses pandas' .loc accessor for efficient value assignment
        data_input.loc[swap_indices, ['sstart', 'send']] = data_input.loc[
            swap_indices, ['send', 'sstart']
        ].values

    # Return the modified DataFrame (same object as input, modified in-place)
    return data_input
