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
from .aesthetics import print_message_box


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
        # Only left side is satisfied:
        ## Left side --> can be extended still
        ## Right side --> Can't be extended, limit
        if previous_extension == 'right':
            # Previously, was extended through the right only (left was banned). So there's no more extension for the right part
            result_data.loc[result_data.index[0], 'send'] = end_coordinate_with_query
            return result_data, None
        elif previous_extension == 'left':
            # If previously was extended on left, send coordinate should stay as it entered
            result_data.loc[result_data.index[0], 'send'] = extended_end
            return result_data, "left"
        else:
            # In this case, "both" will be the previous one, so send coordinate needs to be limited here.
            result_data.loc[result_data.index[0], 'send'] = end_coordinate_with_query
            return result_data, "left"

    elif not left_satisfied and right_satisfied:
        # Only right side is satisfied:
        ## Right side --> can be extended still
        ## Left side --> Can't be extended, limit
        if previous_extension == 'left':
            # Previously, was extended through the left only (right was banned). So there's no more extension for the left part
            result_data.loc[result_data.index[0], 'sstart'] = start_coordinate_with_query
            return result_data, None
        elif previous_extension == 'right':
            # If previously was extended on right, sstart coordinate should stay as it entered
            result_data.loc[result_data.index[0], 'sstart'] = extended_start
            return result_data, "right"
        else:
            # In this case, "both" will be the previous one, so start coordinate needs to be limited here.
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
        max_recursion_depth: int = 10,
        current_depth: int = 0,
        history: list[dict] | None = None,
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
    max_recursion_depth : int, optional
        Maximum allowed recursion depth to prevent infinite recursion. Default is 10.
    current_depth : int, optional
        Current recursion depth level. Used internally for tracking. Default is 0.
    """
    from .blaster import run_blastn_alignment

    # -----------------------------------------------------------------------------
    # STEP 0: Extract row index and initialize result structure
    # -----------------------------------------------------------------------------
    index, element = row_data

    # -----------------------------------------------------------------------------
    # STEP 1: Check recursion depth limit
    # -----------------------------------------------------------------------------
    if current_depth >= max_recursion_depth:
        # Maximum recursion depth reached, return current state
        result = {
            'index': index,
            'modified': False,
            'recursion_info': {
                'max_depth': current_depth,
                'status': 'max_depth_reached'
            },
            'extension_steps': history if history is not None else [{'sstart': int(element['sstart']), 'send': int(element['send'])}]
        }
        return result

    # -----------------------------------------------------------------------------
    # STEP 2: Validate extension direction parameter
    # -----------------------------------------------------------------------------
    if extension_direction not in ["both", "left", "right"]:
        raise ValueError(f"Invalid extension_direction: {extension_direction}. Must be 'both', 'left', or 'right'.")

    # -----------------------------------------------------------------------------
    # STEP 3: Extract data and calculate current sequence length
    # -----------------------------------------------------------------------------
    lower_coor = element['sstart']
    upper_coor = element['send']

    # Initialize and manage history of extension steps (store coordinates only)
    if history is None:
        history = [{'sstart': int(lower_coor), 'send': int(upper_coor)}]

    def _add_step_to_history(sstart_val: int, send_val: int):
        nonlocal history
        if not history or history[-1]['sstart'] != int(sstart_val) or history[-1]['send'] != int(send_val):
            history.append({'sstart': int(sstart_val), 'send': int(send_val)})

    # Calculate the current length of the sequence (inclusive of start and end positions)
    subject_len = upper_coor - lower_coor + 1

    # Initialize a result with default values (assuming no modification needed)
    result = {
        'index': index,
        'modified': False,
        'recursion_info': {
            'max_depth': current_depth,
            'status': 'no_extension_needed'
        },
        'extension_steps': history
    }

    # -----------------------------------------------------------------------------
    # STEP 4: Check if sequence already meets length requirement
    # -----------------------------------------------------------------------------
    if subject_len >= limit_len:
        return result

    # -----------------------------------------------------------------------------
    # STEP 5: Extend the sequence based on the specified direction
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
    # STEP 6: Ensure coordinates remain within valid genome boundaries
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
    # STEP 7: Check if the new length is still < limit_len
    # -----------------------------------------------------------------------------
    if new_subject_len < limit_len:
        # -----------------------------------------------------------------------------
        # STEP 8: Retrieve the extended sequence using BLAST command
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
        # STEP 9: Perform BLAST search with the extended sequence
        # -----------------------------------------------------------------------------
        # Create a temporary FASTA file with the extended sequence
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            temp_fasta.write(f">{index}_extended\n{seq}\n")
            temp_fasta_path = temp_fasta.name

        try:
            # Perform BLAST search with the extended sequence
            blast_results = run_blastn_alignment(
                query_path=temp_fasta_path,
                dict_path=genome_fasta,
                perc_identity=identity,
                word_size=word_size,
                query_coor=True
            )
            
            # -----------------------------------------------------------------------------
            # STEP 10: Filter BLAST results by minimum length
            # -----------------------------------------------------------------------------
            # Filter results by minimum length
            blast_filtered = blast_results[blast_results['len'] >= min_length].copy()
            
            # Clean up temporary file
            os.unlink(temp_fasta_path)

            # -----------------------------------------------------------------------------
            # STEP 11: Filter elements from the blast data that overlap our row data
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
            final_blast_filtered = match_data_and_remove(blast_filtered, elems_to_remove)
            final_blast_filtered.reset_index(drop=True, inplace=True)

            # -----------------------------------------------------------------------------
            # STEP 12: Use next_side_extension_checker to determine further extension
            # -----------------------------------------------------------------------------
            if not final_blast_filtered.empty:
                # Check extension possibilities for the next recursive call:
                ## both: both sides can be extended
                ## left: only left side will be extended
                ## right: only right side can be extended
                ## None: neither side can be extended
                final_row_df, extension_status = next_side_extension_checker(
                    current_element_df,
                    final_blast_filtered,
                    extension_direction
                )
                
                # -----------------------------------------------------------------------------
                # STEP 13: Recursive extension based on extension_status
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
                    # Record this step before recursing
                    _add_step_to_history(int(final_row_df['sstart'].iloc[0]), int(final_row_df['send'].iloc[0]))

                    recursive_result = _process_single_row_extension(
                        row_data=(index, updated_element),
                        genome_fasta=genome_fasta,
                        extend_number=extend_number,
                        limit_len=limit_len,
                        extension_direction=extension_status,
                        identity=identity,
                        word_size=word_size,
                        min_length=min_length,
                        max_recursion_depth=max_recursion_depth,
                        current_depth=current_depth + 1,
                        history=history
                    )
                    
                    # Return the result from the recursive call (it already has recursion info)
                    return recursive_result

                else:
                    # extension_status is None - no further extension possible
                    # Return the coordinates from next_side_extension_checker
                    final_sstart = int(final_row_df['sstart'].iloc[0])
                    final_send = int(final_row_df['send'].iloc[0])
                    final_seq_len = final_send - final_sstart + 1

                    # Record this final step
                    _add_step_to_history(final_sstart, final_send)
                    
                    # Get the final sequence with updated coordinates
                    final_cmd = (
                        f"blastdbcmd -db {genome_fasta} "
                        f"-entry {element['sseqid']} "
                        f"-range {final_sstart}-{final_send} "
                        f"-strand {element['sstrand']} "
                        "-outfmt %s"
                    )
                    final_seq = subprocess.check_output(final_cmd, shell=True, universal_newlines=True).strip()
                    
                    result.update({
                        'modified': True,
                        'len': final_seq_len,
                        'sstart': final_sstart,
                        'send': final_send,
                        'sseq': final_seq,
                        'recursion_info': {
                            'max_depth': current_depth,
                            'status': 'extension_completed'
                        },
                        'extension_steps': history
                    })
                    return result
            else:
                # No BLAST results after filtering - return current extension
                _add_step_to_history(int(lower_coor_extended), int(upper_coor_extended))
                result.update({
                    'modified': True,
                    'len': new_subject_len,
                    'sstart': int(lower_coor_extended),
                    'send': int(upper_coor_extended),
                    'sseq': seq,
                    'recursion_info': {
                        'max_depth': current_depth,
                        'status': 'no_blast_results'
                    },
                    'extension_steps': history
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
        # Get the final sequence with the normal coordinates (not the extended ones)
        final_cmd = (
            f"blastdbcmd -db {genome_fasta} "
            f"-entry {element['sseqid']} "
            f"-range {lower_coor}-{upper_coor} "
            f"-strand {element['sstrand']} "
            "-outfmt %s"
        )
        final_seq = subprocess.check_output(final_cmd, shell=True, universal_newlines=True).strip()

        # Record this step as the final chosen coordinates
        _add_step_to_history(int(lower_coor), int(upper_coor))
        
        result.update({
            'modified': True,
            'len': subject_len,
            'sstart': int(lower_coor),
            'send': int(upper_coor),
            'sseq': final_seq,
            'recursion_info': {
                'max_depth': current_depth,
                'status': 'length_requirement_met'
            },
            'extension_steps': history
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
        max_recursion_depth: int = 10,  # ADD: Maximum recursion depth parameter
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
    
    max_recursion_depth : int, optional
        Maximum allowed recursion depth to prevent infinite recursion. Default is 10.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Process all rows in parallel using joblib
    # -----------------------------------------------------------------------------
    # Distribute the processing of each row across multiple processors
    # Each row is processed independently by the _process_single_row_extension function
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_row_extension)(
            row_data, genome_fasta, extend_number, limit_len, extension_direction,
            identity, word_size, min_length, max_recursion_depth  # ADD: Pass recursion limit
        )
        for row_data in data_input.iterrows()
    )

    # -----------------------------------------------------------------------------
    # STEP 2: Print extension results summary
    # -----------------------------------------------------------------------------
    print("\n=== Extension Results ===")
    extended_count = 0
    recursion_stats = {}
    
    for result in results:
        recursion_info = result.get('recursion_info', {})
        depth = recursion_info.get('max_depth', 0)
        status = recursion_info.get('status', 'unknown')
        
        if result['modified']:
            extended_count += 1
            if depth == 0:
                print(f"Extending row {result['index']} (starting extension - completed)")
            else:
                print(f"Extending row {result['index']} ({depth} recursive calls)")
            
            # Track recursion statistics
            if depth not in recursion_stats:
                recursion_stats[depth] = 0
            recursion_stats[depth] += 1

    # Print summary statistics
    print(f"\nExtension Summary:")
    print(f"  - Total sequences processed: {len(results)}")
    print(f"  - Sequences extended: {extended_count}")
    print(f"  - Sequences not extended: {len(results) - extended_count}")
    
    if recursion_stats:
        print(f"\nRecursion Depth Statistics:")
        for depth in sorted(recursion_stats.keys()):
            if depth == 0:
                print(f"  - No recursion needed: {recursion_stats[depth]} sequences")
            else:
                print(f"  - {depth} recursive calls: {recursion_stats[depth]} sequences")

    # -----------------------------------------------------------------------------
    # STEP 3: Resolve conflicts across sequences using extension histories
    # -----------------------------------------------------------------------------
    # Build accepted intervals progressively in priority order (by index)
    results_by_index = {r['index']: r for r in results}
    ordered_indices = sorted(results_by_index.keys())

    accepted_df = pd.DataFrame(columns=['sseqid', 'sstart', 'send', 'sstrand'])
    chosen_map: dict[int, dict | None] = {}

    adjustments = 0
    removed = 0
    total = len(ordered_indices)

    for idx in ordered_indices:
        row = data_input.loc[idx]
        res = results_by_index[idx]
        steps = res.get('extension_steps') or [{'sstart': int(row['sstart']), 'send': int(row['send'])}]

        # Try from most-extended to least (reverse chronological history)
        chosen_step = None
        for step in reversed(steps):
            candidate = pd.DataFrame([{
                'sseqid': row['sseqid'],
                'sstart': int(step['sstart']),
                'send': int(step['send']),
                'sstrand': row['sstrand']
            }])
            if accepted_df.empty or get_interval_overlap(candidate, accepted_df).empty:
                chosen_step = {'sstart': int(step['sstart']), 'send': int(step['send'])}
                break
        if chosen_step is None:
            # No non-overlapping step exists; remove this sequence entirely
            chosen_map[idx] = None
            removed += 1
            continue

        # Count adjustments if not using the last (most-extended) step
        if steps and (chosen_step['sstart'] != int(steps[-1]['sstart']) or chosen_step['send'] != int(steps[-1]['send'])):
            adjustments += 1

        # Append to accepted set
        accepted_df = pd.concat([accepted_df, pd.DataFrame([{
            'sseqid': row['sseqid'],
            'sstart': chosen_step['sstart'],
            'send': chosen_step['send'],
            'sstrand': row['sstrand']
        }])], ignore_index=True)

        chosen_map[idx] = chosen_step

    # -----------------------------------------------------------------------------
    # STEP 4: Update DataFrame with chosen coordinates and sequences
    # -----------------------------------------------------------------------------
    removed_indices: list[int] = []
    for idx in ordered_indices:
        chosen_step = chosen_map[idx]
        if chosen_step is None:
            removed_indices.append(idx)
            continue
        row = data_input.loc[idx]
        chosen_len = int(chosen_step['send']) - int(chosen_step['sstart']) + 1
        # Fetch sequence for chosen step
        final_cmd = (
            f"blastdbcmd -db {genome_fasta} "
            f"-entry {row['sseqid']} "
            f"-range {int(chosen_step['sstart'])}-{int(chosen_step['send'])} "
            f"-strand {row['sstrand']} "
            "-outfmt %s"
        )
        final_seq = subprocess.check_output(final_cmd, shell=True, universal_newlines=True).strip()

        data_input.loc[idx, 'len'] = chosen_len
        data_input.loc[idx, 'sstart'] = int(chosen_step['sstart'])
        data_input.loc[idx, 'send'] = int(chosen_step['send'])
        data_input.loc[idx, 'sseq'] = final_seq

    # Drop removed sequences and reset index for cleanliness
    if removed_indices:
        data_input = data_input.drop(index=removed_indices).reset_index(drop=True)

    # Print conflict resolution summary
    print(f"\nConflict resolution summary:")
    print(f"  - Total sequences considered: {total}")
    print(f"  - Sequences adjusted due to conflicts: {adjustments}")
    print(f"  - Sequences removed due to irresolvable conflicts: {removed}")

    # Return the updated DataFrame with resolved sequences
    return data_input
