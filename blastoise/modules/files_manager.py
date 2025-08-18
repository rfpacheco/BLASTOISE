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
"""

import pandas as pd
import subprocess
from typing import Optional, Dict, Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed


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

