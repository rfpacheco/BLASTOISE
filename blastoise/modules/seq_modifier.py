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
# noinspection PyPackageRequirements
import pandas as pd
import tempfile
import os
from typing import Dict, Any, Tuple
# noinspection PyPackageRequirements
from joblib import Parallel, delayed
from .genomic_ranges import merge_intervals, get_interval_overlap
from .strand_location import match_data_and_remove


def next_side_extension_checker(
        data_input: pd.DataFrame,
        query_data: pd.DataFrame,
        previous_extension: str = None
) -> tuple[pd.DataFrame, str|None]:
    """
    Check if query data satisfies the left and/or right sides of the original sequence coordinates.

    This function analyzes query sequences to determine if they completely cover the left side,
    right side, or both sides of an original sequence. Based on the coverage and the previous
    extension direction, it returns the modified coordinates and indicates which side(s) were satisfied.

    The previous_extension parameter constrains the output:
    - If previous_extension was 'right', even if both sides are satisfied, only 'right' can be returned
    - If previous_extension was 'left', even if both sides are satisfied, only 'left' can be returned

    Parameters
    ----------
    data_input : pd.DataFrame
        A DataFrame with one row containing the original coordinates of a sequence.
        Must have 'sstart' and 'send' columns.
    query_data : pd.DataFrame
        A DataFrame with other sequences containing query data.
        Must have 'qstart' and 'qend' columns.
    previous_extension : str, optional
        The direction of the previous extension ('left', 'right', or None).
        Constrains the possible return values when both sides are satisfied.

    Returns
    -------
    tuple[pd.DataFrame, str]
        A tuple containing:
        - Modified data_input DataFrame with updated coordinates if needed
        - String indicating satisfaction: "both", "left", "right", or None
    """

    # Extract start and end coordinates from the single row in data_input
    extended_start = data_input['sstart'].iloc[0]
    extended_end = data_input['send'].iloc[0]

    # If query_data is empty, return original data with None
    if query_data.empty:
        return data_input.copy(), None

    # Merge overlapping intervals in query data
    query_merge = merge_intervals(query_data, start_col='qstart', end_col='qend')

    # Get the range covered by query data
    query_left_value = query_merge['qstart'].min()
    query_right_value = query_merge['qend'].max()

    start_coordinate_with_query = extended_start + query_left_value - 1
    end_coordinate_with_query = extended_start + query_right_value - 1

    # Check coverage on both sides
    left_satisfied = start_coordinate_with_query == extended_start
    right_satisfied = end_coordinate_with_query == extended_end

    # Create a copy of the input data for modification
    result_data = data_input.copy()

    if left_satisfied and right_satisfied:
        # Both sides are completely satisfied
        # But check previous_extension constraint
        if previous_extension == 'right':
            # Previous extension was right, so left part was already discarded
            # Can only return 'right'
            return result_data, "right"
        elif previous_extension == 'left':
            # Previous extension was left, so right part was already discarded
            # Can only return 'left'
            return result_data, "left"
        else:
            # No previous extension constraint. Can return 'both'
            return result_data, "both"

    elif left_satisfied and not right_satisfied:
        # Only left side is satisfied, modify coordinates to where satisfaction ends
        # Check previous extension
        if previous_extension == 'right':
            # Modify right coordinate
            result_data.loc[result_data.index[0], 'send'] = end_coordinate_with_query
            return result_data, None
        else:
            result_data.loc[result_data.index[0], 'send'] = end_coordinate_with_query
            return result_data, "left"

    elif not left_satisfied and right_satisfied:
        # Only right side is satisfied, modify coordinates to where satisfaction begins
        # Check previous extensions
        if previous_extension == 'left':
            # Modify left coordinate
            result_data.loc[result_data.index[0], 'sstart'] = start_coordinate_with_query
            return result_data, None
        else:
            result_data.loc[result_data.index[0], 'sstart'] = start_coordinate_with_query
            return result_data, "right"

    else:
        # Neither side is completely satisfied
        # Modify coordinates to the range that was satisfied
        result_data.loc[result_data.index[0], 'sstart'] = start_coordinate_with_query
        result_data.loc[result_data.index[0], 'send'] = end_coordinate_with_query
        return result_data, None


def _process_single_row_extension(
        row_data: Tuple[Any, pd.Series],
        genome_fasta: str,
        extend_number: int,
        limit_len: int,
        extension_direction: str = "both",
        identity: int = 60,
        word_size: int = 15,
        min_length: int = 100,
) -> Dict[str, Any]:
    """
    Process a single row for sequence extension with recursive extension capability.

    This helper function processes a single row from the DataFrame, evaluates the sequence length,
    and extends the sequence if it's shorter than the specified limit. The extension can be done by
    adjusting the start and end coordinates on both sides, only on the left side, or only on the 
    right side, ensuring they don't go beyond valid genome boundaries. After extension, it performs 
    a BLAST search and uses the results to determine if further extension is possible through 
    recursive calls based on the next_side_extension_checker analysis.

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
    extension_direction : str, optional
        Direction for extension. Options are:
        - "both": Extend on both sides (default)
        - "left": Extend only on the left side (5' end)
        - "right": Extend only on the right side (3' end)
    identity : int, optional
        Identity percentage for BLAST search. Default is 60.
    word_size : int, optional
        Word size for BLAST search. Default is 15.
    min_length : int, optional
        Minimum sequence length for filtering BLAST results. Default is 100.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the processed row data with keys:
        - 'index': The original row index
        - 'modified': Boolean flag indicating if the sequence was modified
        - 'length': The final sequence length (only if modified)
        - 'sstart': The final start coordinate (only if modified)
        - 'send': The final end coordinate (only if modified)
        - 'sseq': The final extended sequence string (only if modified)

    Raises
    ------
    subprocess.CalledProcessError
        If the BLAST command to retrieve sequence information fails.
    KeyError
        If the required columns are missing from the input row data.
    ValueError
        If extension_direction is not one of "both", "left", or "right".
    """
    from .blaster import blastn_blaster

    # -----------------------------------------------------------------------------
    # STEP 1: Validate extension direction parameter
    # -----------------------------------------------------------------------------
    if extension_direction not in ["both", "left", "right"]:
        raise ValueError(f"Invalid extension_direction: {extension_direction}. Must be 'both', 'left', or 'right'.")

    # -----------------------------------------------------------------------------
    # STEP 2: Extract data and calculate current sequence length
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
    # STEP 3: Check if sequence already meets length requirement
    # -----------------------------------------------------------------------------
    if subject_len >= limit_len:
        return result

    # -----------------------------------------------------------------------------
    # STEP 4: Extend the sequence based on the specified direction
    # -----------------------------------------------------------------------------
    lower_coor_extended = int  # Initialize new variable
    upper_coor_extended = int  # Initialize new variable
    if extension_direction == "both":
        # Extend equally on both sides
        lower_coor_extended = lower_coor - extend_number  # Extend at the 5' end
        upper_coor_extended = upper_coor + extend_number  # Extend at the 3' end
    elif extension_direction == "left":
        # Extend only on the left side (5' end)
        lower_coor_extended = lower_coor - extend_number
        upper_coor_extended = upper_coor
    elif extension_direction == "right":
        lower_coor_extended = lower_coor
        # Extend only on the right side (3' end)
        upper_coor_extended = upper_coor + extend_number

    # -----------------------------------------------------------------------------
    # STEP 5: Ensure coordinates remain within valid genome boundaries
    # -----------------------------------------------------------------------------
    # Ensure lower coordinate is not negative (BLAST coordinates start at 1)
    if lower_coor_extended <= 0:
        lower_coor_extended = 1

    # Ensure upper coordinate doesn't exceed chromosome length
    ## Get the maximum length of the chromosome to avoid exceeding boundaries
    cmd = f"blastdbcmd -db {genome_fasta} -entry {element['sseqid']} -outfmt '%l'"
    chrom_max_len = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()
    chrom_max_len = int(chrom_max_len)

    ## Now, ensure the coordinate
    if upper_coor_extended > chrom_max_len:
        upper_coor_extended = chrom_max_len

    # Calculate the new sequence length after boundary adjustments
    new_subject_len = upper_coor_extended - lower_coor_extended + 1

    # -----------------------------------------------------------------------------
    # STEP 6: Check if the new length is still < limit_len
    # -----------------------------------------------------------------------------
    if new_subject_len < limit_len:
        # -----------------------------------------------------------------------------
        # STEP 7: Retrieve the extended sequence using BLAST command
        # -----------------------------------------------------------------------------
        # Construct the BLAST command to extract the sequence
        cmd = (
            f"blastdbcmd -db {genome_fasta} "
            f"-entry {element['sseqid']} "
            f"-range {lower_coor_extended}-{upper_coor_extended} "
            f"-strand {element['sstrand']} "
            "-outfmt %s"
        )

        # Execute the command to get the sequence
        seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

        # -----------------------------------------------------------------------------
        # STEP 8: Perform BLAST search with the extended sequence
        # -----------------------------------------------------------------------------
        # Create a temporary FASTA file with the extended sequence
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            temp_fasta.write(f">{index}_extended\n{seq}\n")
            temp_fasta_path = temp_fasta.name

        try:
            # Perform BLAST search with the extended sequence
            blast_results = blastn_blaster(
                query_path=temp_fasta_path,
                dict_path=genome_fasta,
                perc_identity=identity,
                word_size=word_size,
                query_coor=True
            )
            
            # -----------------------------------------------------------------------------
            # STEP 9: Filter BLAST results by minimum length
            # -----------------------------------------------------------------------------
            # Calculate sequence length for filtering (if not already present)
            if 'len' not in blast_results.columns:  # TODO: check if this is needed --> REMOVE, already checked
                blast_results['len'] = blast_results['send'] - blast_results['sstart'] + 1
            
            # Filter results by minimum length
            blast_filtered = blast_results[blast_results['len'] >= min_length].copy()
            
            # Clean up temporary file
            os.unlink(temp_fasta_path)  # TODO: does this remove the file?

            # -----------------------------------------------------------------------------
            # STEP 10: Filter elements from the blast data that overlap our row data
            # -----------------------------------------------------------------------------
            # Create a DataFrame with the current element for checking
            current_element_df = pd.DataFrame([{
                'sstart': lower_coor_extended,
                'send': upper_coor_extended,
                'sseqid': element['sseqid'],
                'sstrand': element['sstrand']
            }])

            elems_to_remove = get_interval_overlap(
                blast_filtered,
                current_element_df
            )
            final_blast_filtered = match_data_and_remove(blast_filtered, elems_to_remove)  # TODO: reset index in function?
            final_blast_filtered.reset_index(drop=True, inplace=True)

            # -----------------------------------------------------------------------------
            # STEP 10: Use next_side_extension_checker to determine further extension
            # -----------------------------------------------------------------------------
            if not final_blast_filtered.empty:
                # Check extension possibilities
                final_row_df, extension_status = next_side_extension_checker(
                    current_element_df,
                    final_blast_filtered,
                    extension_direction
                )
                
                # -----------------------------------------------------------------------------
                # STEP 11: Recursive extension based on extension_status
                # -----------------------------------------------------------------------------
                if extension_status in ["both", "left", "right"]:
                    # Further extension is possible, prepare for recursive call
                    updated_element = pd.Series({
                        'sstart': final_row_df['sstart'].iloc[0],
                        'send': final_row_df['send'].iloc[0],
                        'sseqid': element['sseqid'],
                        'sstrand': element['sstrand']
                    })

                    # Recursive call with updated coordinates, direction, and previous extension tracking
                    recursive_result = _process_single_row_extension(
                        row_data=(index, updated_element),
                        genome_fasta=genome_fasta,
                        extend_number=extend_number,
                        limit_len=limit_len,
                        extension_direction=extension_status,
                        identity=identity,
                        word_size=word_size,
                        min_length=min_length
                    )
                    
                    # Return the result from the recursive call
                    if recursive_result['modified']:
                        return recursive_result
                    else:
                        # If recursive call didn't modify (reached limit), return current state
                        result.update({
                            'modified': True,
                            'len': final_row_df['send'].iloc[0] - final_row_df['sstart'].iloc[0] + 1,
                            'sstart': int(final_row_df['sstart'].iloc[0]),
                            'send': int(final_row_df['send'].iloc[0]),
                            'sseq': seq  # Note: This might need to be updated to reflect final coordinates
                        })
                        return result
                else:
                    # extension_status is None - no further extension possible
                    # Return the coordinates from next_side_extension_checker
                    final_seq_len = final_row_df['send'].iloc[0] - final_row_df['sstart'].iloc[0] + 1
                    
                    # Get the final sequence with updated coordinates
                    final_cmd = (
                        f"blastdbcmd -db {genome_fasta} "
                        f"-entry {element['sseqid']} "
                        f"-range {int(final_row_df['sstart'].iloc[0])}-{int(final_row_df['send'].iloc[0])} "
                        f"-strand {element['sstrand']} "
                        "-outfmt %s"
                    )
                    final_seq = subprocess.check_output(final_cmd, shell=True, universal_newlines=True).strip()
                    
                    result.update({
                        'modified': True,
                        'len': final_seq_len,
                        'sstart': int(final_row_df['sstart'].iloc[0]),
                        'send': int(final_row_df['send'].iloc[0]),
                        'sseq': final_seq
                    })
                    return result
            else:
                # No BLAST results after filtering - return current extension
                result.update({
                    'modified': True,
                    'len': new_subject_len,
                    'sstart': int(lower_coor_extended),
                    'send': int(upper_coor_extended),
                    'sseq': seq
                })
                return result
                
        except Exception as e:
            # Clean up temporary file in case of error
            if os.path.exists(temp_fasta_path):
                os.unlink(temp_fasta_path)
            # Re-raise the exception
            raise e
    else:
        # Sequence length after extension meets or exceeds limit_len
        # Get the final sequence
        final_cmd = (
            f"blastdbcmd -db {genome_fasta} "
            f"-entry {element['sseqid']} "
            f"-range {lower_coor_extended}-{upper_coor_extended} "
            f"-strand {element['sstrand']} "
            "-outfmt %s"
        )
        final_seq = subprocess.check_output(final_cmd, shell=True, universal_newlines=True).strip()
        
        result.update({
            'modified': True,
            'len': new_subject_len,
            'sstart': int(lower_coor_extended),
            'send': int(upper_coor_extended),
            'sseq': final_seq
        })
        return result


def sequence_extension(
        data_input: pd.DataFrame,
        genome_fasta: str,
        extend_number: int,
        limit_len: int,
        identity: int = 60,
        word_size: int = 15,
        min_length: int = 100,
        extension_direction: str = "both",
        n_jobs: int = -1
) -> pd.DataFrame:
    """
    Adjusts subject sequences in a data frame to meet a specified length requirement.

    This function processes a DataFrame of genomic sequences, identifying those that are
    shorter than a specified threshold and extending them to meet minimum length requirements.
    The extension can be performed by adjusting the start and end coordinates on both sides,
    only on the left side, or only on the right side while ensuring they remain within 
    valid genome boundaries. After extension, it performs BLAST searches with the extended
    sequences and filters results by minimum length.

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
    extension_direction : str, optional
        Direction for extension. Options are:
        - "both": Extend on both sides (default)
        - "left": Extend only on the left side (5' end)
        - "right": Extend only on the right side (3' end)
    identity : int, optional
        Identity percentage for BLAST search. Default is 60.
    word_size : int, optional
        Word size for BLAST search. Default is 15.
    min_length : int, optional
        Minimum sequence length for filtering BLAST results. Default is 100.
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
        Additional information about BLAST results is stored but not directly added
        to the returned DataFrame.

    Raises
    ------
    subprocess.CalledProcessError
        If an error occurs when running the system command to retrieve the genome sequence.
    KeyError
        If any required column is missing in the provided data frame.
    ValueError
        If extension_direction is not one of "both", "left", or "right".

    Notes
    -----
    The function assumes that the input 'data_input' follows a specific structure with
    the necessary columns to process coordinate adjustments and sequence extraction.
    Coordinates are adjusted considering both strands ('plus' and 'minus') with careful 
    attention to avoid exceeding genome boundaries.
    
    BLAST results from the extension process are computed and filtered but are currently
    stored in the processing results for potential future use rather than being directly
    integrated into the returned DataFrame.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Process all rows in parallel using joblib
    # -----------------------------------------------------------------------------
    # Distribute the processing of each row across multiple processors
    # Each row is processed independently by the _process_single_row_extension function
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_row_extension)(
            row_data, genome_fasta, extend_number, limit_len, extension_direction,
            identity, word_size, min_length
        )
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
            data_input.loc[index, 'len'] = result['len']
            data_input.loc[index, 'sstart'] = result['sstart']
            data_input.loc[index, 'send'] = result['send']
            data_input.loc[index, 'sseq'] = result['sseq']

            # Note: BLAST results are available in result['blast_results'] 
            # but are not currently integrated into the main DataFrame
            # This could be implemented in future iterations if needed

    # Return the updated DataFrame with extended sequences
    return data_input
