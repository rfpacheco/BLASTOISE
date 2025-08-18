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
from .files_manager import end_always_greater_than_start, get_data_sequence
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
        word_size: int = 15,
        query_coor: bool = False
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
    query_coor : bool, optional
        Whether to include query coordinates (qstart and qend) in the output.
        Default is True.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the parsed and processed BLAST results with columns:
        'qseqid' (query sequence ID), 'sseqid' (subject sequence ID),
        'qstart' (query start position, if query_coor=True), 'qend' (query end position, if query_coor=True),
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

    # Build the command with conditional query coordinates
    if query_coor:
        outfmt = "'10 qseqid sseqid qstart qend sstart send sstrand evalue sseq'"
        columns = ['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'evalue', 'sseq']
        numeric_columns = ['qstart', 'qend', 'sstart', 'send']
        final_column_order = ['qseqid', 'sseqid', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'evalue', 'sseq',
                              'len']
    else:
        outfmt = "'10 qseqid sseqid sstart send sstrand evalue sseq'"
        columns = ['qseqid', 'sseqid', 'sstart', 'send', 'sstrand', 'evalue', 'sseq']
        numeric_columns = ['sstart', 'send']
        final_column_order = ['qseqid', 'sseqid', 'sstart', 'send', 'sstrand', 'evalue', 'sseq', 'len']

    cmd = (
        f"blastn -word_size {word_size} "
        f"-query {query_path} "
        f"-db {dict_path} "
        f"-perc_identity {perc_identity} "
        f"-outfmt {outfmt}"
    )

    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    if not data:
        return pd.DataFrame(columns=columns)
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    data.columns = columns

    # Convert coordinate columns to int type 
    data[numeric_columns] = data[numeric_columns].astype(int)

    # get 'evalue' as a 'float' type
    data['evalue'] = data['evalue'].astype(float)

    # Make sure 'send' > 'sstart'
    data = end_always_greater_than_start(data)

    # Create 'len' column
    data['len'] = data['send'] - data['sstart'] + 1

    # Reorder columns for better readability
    data = data[final_column_order]

    return data
