"""
BLASTOISE Module: Sequence Extension and Modification
====================================================

This module provides functionality for extending and modifying genomic sequences
in the BLASTOISE pipeline. It handles the adjustment of sequence coordinates to
ensure sequences meet minimum length requirements while respecting genome boundaries.

The module contains two main functions:
1. `_process_single_row_extension`: A helper function that processes a single row
   for sequence extension.
2. `sequence_extension`: The main function that adjusts subject sequences in a data
   frame to meet specified length requirements.

These functions work together to ensure that sequences used in the BLASTOISE pipeline
are of sufficient length for reliable analysis, extending shorter sequences by
adjusting their coordinates and retrieving the extended sequences from the genome.

Author: R. Pacheco
"""

import subprocess
import pandas as pd
from typing import Dict, Any, Tuple
from joblib import Parallel, delayed


def _process_single_row_extension(
        row_data: Tuple[Any, pd.Series],
        genome_fasta: str,
        extend_number: int,
        limit_len: int
) -> Dict[str, Any]:
    """
    Process a single row for sequence extension.

    This helper function processes a single row from the DataFrame, evaluates the sequence length,
    and extends the sequence if it's shorter than the specified limit. The extension is done by
    adjusting the start and end coordinates equally on both sides, ensuring they don't go beyond 
    valid genome boundaries, and then retrieving the extended sequence using a BLAST command.

    The function performs several key steps:
    1. Extracts the index and data from the input tuple
    2. Calculates the current sequence length
    3. If the sequence is shorter than the limit, extends it by adjusting coordinates
    4. Ensures coordinates remain within valid genome boundaries
    5. Retrieves the extended sequence using blastdbcmd
    6. Returns the modified data or indicates no modification was needed

    Parameters
    ----------
    row_data : Tuple[Any, pd.Series]
        A tuple containing the index and the row data from the DataFrame. The row data
        must contain 'sstart', 'send', 'sseqid', and 'sstrand' columns.
    genome_fasta : str
        The path to the genome FASTA file used for extracting sequence data.
        This should be a BLAST database created with makeblastdb.
    extend_number : int
        The number of nucleotides to add on each side of the sequence.
        For example, if extend_number=100, the sequence will be extended by 
        100 nucleotides on both the 5' and 3' ends.
    limit_len : int
        A length threshold that triggers sequence extension. If the sequence
        is already longer than this value, no extension is performed.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the processed row data with keys:
        - 'index': The original row index
        - 'modified': Boolean flag indicating if the sequence was modified
        - 'length': The new sequence length (only if modified)
        - 'sstart': The new start coordinate (only if modified)
        - 'send': The new end coordinate (only if modified)
        - 'sseq': The extended sequence string (only if modified)

    Raises
    ------
    subprocess.CalledProcessError
        If the BLAST command to retrieve sequence information fails.
    KeyError
        If the required columns are missing from the input row data.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Extract data and calculate current sequence length
    # -----------------------------------------------------------------------------
    index, element = row_data
    lower_coor = element['sstart']
    upper_coor = element['send']

    # Calculate the current length of the sequence (inclusive of start and end positions)
    subject_len = upper_coor - lower_coor + 1

    # Initialize a result with default values (assuming no modification needed)
    result = {
        'index': index,
        'modified': False
    }

    # -----------------------------------------------------------------------------
    # STEP 2: Check if a sequence needs extension and adjust coordinates if necessary
    # -----------------------------------------------------------------------------
    if subject_len < limit_len:
        # Extend the sequence equally on both sides
        lower_coor = lower_coor - extend_number  # Extend at the 5' end
        upper_coor = upper_coor + extend_number  # Extend at the 3' end

        # -----------------------------------------------------------------------------
        # STEP 3: Ensure coordinates remain within valid genome boundaries
        # -----------------------------------------------------------------------------
        # Ensure lower coordinate is not negative (BLAST coordinates start at 1)
        if lower_coor <= 0:
            lower_coor = 1

        # Get the maximum length of the chromosome to avoid exceeding boundaries
        cmd = f"blastdbcmd -db {genome_fasta} -entry {element['sseqid']} -outfmt '%l'"
        chrom_max_len = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()
        chrom_max_len = int(chrom_max_len)

        # Ensure upper coordinate doesn't exceed chromosome length
        if upper_coor > chrom_max_len:
            upper_coor = chrom_max_len

        # Calculate the new sequence length after boundary adjustments
        new_subject_len = upper_coor - lower_coor + 1

        # -----------------------------------------------------------------------------
        # STEP 4: Check if the new length is < limit_len
        # -----------------------------------------------------------------------------
        if new_subject_len < limit_len:
            # -----------------------------------------------------------------------------
            # STEP 5: Retrieve the extended sequence using BLAST command
            # -----------------------------------------------------------------------------
            # Construct the BLAST command to extract the sequence
            cmd = (
                f"blastdbcmd -db {genome_fasta} "
                f"-entry {element['sseqid']} "
                f"-range {lower_coor}-{upper_coor} "
                f"-strand {element['sstrand']} "
                "-outfmt %s"
            )

            # Execute the command to get the sequence
            seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

            # -----------------------------------------------------------------------------
            # STEP 6: Update the result with the modified sequence information
            # -----------------------------------------------------------------------------
            result.update({
                'modified': True,
                'length': new_subject_len,
                'sstart': int(lower_coor),
                'send': int(upper_coor),
                'sseq': seq
            })

    # Return either the modified data or the default "not modified" result
    return result


def sequence_extension(
        data_input: pd.DataFrame,
        genome_fasta: str,
        extend_number: int,
        limit_len: int,
        n_jobs: int = -1
) -> pd.DataFrame:
    """
    Adjusts subject sequences in a data frame to meet a specified length requirement.

    This function processes a DataFrame of genomic sequences, identifying those that are
    shorter than a specified threshold and extending them to meet minimum length requirements.
    The extension is performed by adjusting the start and end coordinates equally on both
    sides while ensuring they remain within valid genome boundaries.

    The function uses parallel processing to efficiently handle large datasets, distributing
    the workload across multiple processors. Each sequence is processed independently, and
    the results are then integrated back into the original DataFrame.

    Parameters
    ----------
    data_input : pd.DataFrame
        A data frame containing sequence alignment details. Expected to have the following 
        columns: 'sstrand', 'sstart', 'send', 'sseqid', and alignment metrics such as
        'pident', 'length', 'qstart', etc.
    genome_fasta : str
        The path to the genome FASTA file used for extracting sequence data.
        This should be a BLAST database created with makeblastdb.
    extend_number : int 
        The number of nucleotides to add on each side of the sequence.
        For example, if extend_number=100, the sequence will be extended by 
        100 nucleotides on both the 5' and 3' ends.
    limit_len : int
        A length threshold that triggers sequence extension. If a sequence
        is already longer than this value, no extension is performed.
    n_jobs : int, optional
        The number of jobs to run in parallel. -1 means using all processors.
        Default is -1.

    Returns
    -------
    pd.DataFrame
        A modified data frame where sequences having lengths below the specified minimum
        are extended. The following columns are updated for extended sequences:
        - 'length': The new sequence length
        - 'sstart': The new start coordinate
        - 'send': The new end coordinate
        - 'sseq': The extended sequence string

    Raises
    ------
    subprocess.CalledProcessError
        If an error occurs when running the system command to retrieve the genome sequence.
    KeyError
        If any required column is missing in the provided data frame.

    Notes
    -----
    The function assumes that the input 'data_input' follows a specific structure with
    the necessary columns to process coordinate adjustments and sequence extraction.
    Coordinates are adjusted considering both strands ('plus' and 'minus') with careful 
    attention to avoid exceeding genome boundaries.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Process all rows in parallel using joblib
    # -----------------------------------------------------------------------------
    # Distribute the processing of each row across multiple processors
    # Each row is processed independently by the _process_single_row_extension function
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_row_extension)(row_data, genome_fasta, extend_number, limit_len)
        for row_data in data_input.iterrows()
    )

    # -----------------------------------------------------------------------------
    # STEP 2: Update the DataFrame with the results from parallel processing
    # -----------------------------------------------------------------------------
    # Iterate through the results and update only the rows that were modified
    for result in results:
        if result['modified']:
            index = result['index']
            # Update the relevant columns with the new values
            data_input.loc[index, 'length'] = result['length']
            data_input.loc[index, 'sstart'] = result['sstart']
            data_input.loc[index, 'send'] = result['send']
            data_input.loc[index, 'sseq'] = result['sseq']

    # Return the updated DataFrame with extended sequences
    return data_input
