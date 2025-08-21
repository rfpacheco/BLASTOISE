"""
BLASTOISE Module: BLAST Operations and Sequence Processing
=========================================================

This module provides core functionality for the BLASTOISE pipeline, handling BLAST database
creation and sequence alignment. It serves as the engine for identifying repetitive genomic
elements through a series of BLAST operations and data processing steps.

Main functions:
1. create_blast_database: Creates a BLAST database from a genome FASTA file
2. run_blastn_alignment: Performs BLASTn alignment and returns results as a DataFrame

Author: R. Pacheco
"""

import pandas as pd
import subprocess
import logging
from typing import List


def create_blast_database(path_input: str, path_output: str) -> None:
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
        Path where the database files will be stored (prefix for BLAST DB files).

    Raises
    ------
    RuntimeError
        If the BLAST database creation fails.

    Notes
    -----
    - Standard output and error from makeblastdb are suppressed to keep the console clean.
    - This function uses a non-shell invocation to avoid shell quoting issues and sets
      check=True to surface failures.
    """

    # "-parse_seqids" preserves the original sequence IDs in the database.
    cmd: List[str] = [
        "makeblastdb",
        "-in", path_input,
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", path_output,
    ]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        logging.error("makeblastdb failed to create BLAST database.")
        raise RuntimeError(f"Failed to create BLAST database for '{path_input}' at '{path_output}'.") from e


def run_blastn_alignment(
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
        Path to the pre-built BLAST database (created with create_blast_database).
    perc_identity : float
        a Minimum percentage identity threshold for reporting matches (0-100).
    word_size : int, optional
        Size of the word used for seeding alignments in the BLAST algorithm.
        Default is 15.
    query_coor : bool, optional
        Whether to include query coordinates (qstart and qend) in the output.
        Default is False.

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
    outfmt: str
    columns: List[str]; numeric_columns: List[str]; final_column_order: List[str]
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

    cmd: str = (
        f"blastn -word_size {word_size} "
        f"-query {query_path} "
        f"-db {dict_path} "
        f"-perc_identity {perc_identity} "
        f"-outfmt {outfmt}"
    )

    data: str = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    if not data:
        return pd.DataFrame(columns=columns)
    data_df: pd.DataFrame = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    data_df.columns = columns

    # Convert coordinate columns to int type 
    data_df[numeric_columns] = data_df[numeric_columns].astype(int)

    # get 'evalue' as a 'float' type
    data_df['evalue'] = data_df['evalue'].astype(float)

    # Make sure 'send' > 'sstart'
    mask: pd.Series = data_df['sstart'] > data_df['send']
    data_df.loc[mask, ['sstart', 'send']] = data_df.loc[mask, ['send', 'sstart']].values

    # Create 'len' column
    data_df['len'] = data_df['send'] - data_df['sstart'] + 1

    # Reorder columns for better readability
    data_df = data_df[final_column_order]

    return data_df
