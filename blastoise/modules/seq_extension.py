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
from typing import Dict, Any, Tuple, List
# noinspection PyPackageRequirements
from joblib import Parallel, delayed
from .genomic_ranges import merge_overlapping_intervals, get_interval_overlap
from .strand_location import match_data_and_remove
from .blaster import run_blastn_alignment


def next_side_extension_checker(
        data_input: pd.DataFrame,
        query_data: pd.DataFrame,
        previous_extension: str = None
) -> tuple[pd.DataFrame, str|None]:
    """
    Checks the next side extension for a given data input, based on overlaps with
    query data and previous extension constraints.

    This function evaluates whether the left or right sides of the genomic range
    in the input data are completely satisfied based on the given query data.
    Depending on the conditions, it adjusts the start and end coordinates of the
    range in the input data and determines the next permissible extension for the
    input data.

    Parameters
    ----------
    data_input : pd.DataFrame
        A DataFrame with genomic range information. It is assumed to contain
        columns 'sstart' and 'send', representing the start and end coordinates
        of the range.

    query_data : pd.DataFrame
        A DataFrame with query data for comparison with `data_input`. It is
        assumed to contain columns 'qstart' and 'qend', representing the start
        and end coordinates of the query intervals.

    previous_extension : str or None, optional
        A string indicating the side of the previous extension. It can be
        'left', 'right', or None if no previous extension constraint exists.

    Returns
    -------
    tuple[pd.DataFrame, str or None]
        A tuple where:
        - The first element is a modified DataFrame (a copy of `data_input`)
          adjusted based on the given conditions.
        - The second element is a string indicating the next permissible
          extension ('left', 'right', 'both'), or None if no further
          extension is possible.
    """

    # Extract start and end coordinates from the single row in data_input
    extended_start: int = data_input.iloc[0].sstart
    extended_end: int = data_input.iloc[0].send

    # Merge overlapping intervals in query data
    query_merge: pd.DataFrame = merge_overlapping_intervals(query_data, start_col='qstart', end_col='qend')

    # Get the minimun and maximun values query data gives as information
    query_left_value: int = query_merge['qstart'].min()
    query_right_value: int = query_merge['qend'].max()

    # With the query data, get the coverage ranges
    start_coordinate_with_query: int = extended_start + query_left_value - 1
    end_coordinate_with_query: int = extended_start + query_right_value - 1

    # Check coverage on both sides
    left_satisfied: bool = start_coordinate_with_query == extended_start
    right_satisfied: bool = end_coordinate_with_query == extended_end

    # Create a copy of the input data for modification
    result_data: pd.DataFrame = data_input.copy()

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
        # Neither side is completely satisfied.
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
    # from .blaster import run_blastn_alignment  # TODO: is it needed?

    # -----------------------------------------------------------------------------
    # STEP 0: Extract row index and initialize result structure
    # -----------------------------------------------------------------------------
    index: int
    element: pd.Series
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
            }
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
    lower_coor: int = element.sstart
    upper_coor: int = element.send


    # Calculate the current length of the sequence (inclusive of start and end positions)
    subject_len: int = upper_coor - lower_coor + 1

    # Initialize a result with default values (assuming no modification needed)
    result = {
        'index': index,
        'modified': False,
        'recursion_info': {
            'max_depth': current_depth,
            'status': 'no_extension_needed'
        }
    }

    # -----------------------------------------------------------------------------
    # STEP 4: Check if sequence already meets length requirement
    # -----------------------------------------------------------------------------
    # If the max length is reached, return current status
    if subject_len >= limit_len:
        return result

    # -----------------------------------------------------------------------------
    # STEP 5: Extend the sequence based on the specified direction
    # -----------------------------------------------------------------------------
    lower_coor_extended = int
    upper_coor_extended = int
    if extension_direction == "both":
        # Extend equally on both sides
        lower_coor_extended = lower_coor - extend_number  # Extend at the 5' end
        upper_coor_extended = upper_coor + extend_number  # Extend at the 3' end
    elif extension_direction == "left":
        lower_coor_extended = lower_coor - extend_number  # Extend only on the left side (5' end)
        upper_coor_extended = upper_coor
    elif extension_direction == "right":
        lower_coor_extended = lower_coor
        upper_coor_extended = upper_coor + extend_number  # Extend only on the right side (3' end)

    # -----------------------------------------------------------------------------
    # STEP 6: Ensure coordinates remain within valid genome boundaries
    # -----------------------------------------------------------------------------
    # Ensure lower coordinate is not negative (BLAST coordinates start at 1)
    if lower_coor_extended <= 0:
        lower_coor_extended = 1

    # Ensure upper coordinate doesn't exceed chromosome length
    ## Get the maximum length of the chromosome to avoid exceeding boundaries
    cmd: str = f"blastdbcmd -db {genome_fasta} -entry {element['sseqid']} -outfmt '%l'"
    chrom_max_len: int = int(subprocess.check_output(cmd, shell=True, universal_newlines=True).strip())

    ## Now, ensure the coordinate
    if upper_coor_extended > chrom_max_len:
        upper_coor_extended = chrom_max_len

    # Calculate the new sequence length after boundary adjustments
    new_subject_len: int = upper_coor_extended - lower_coor_extended + 1

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
        seq: str = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

        # -----------------------------------------------------------------------------
        # STEP 9: Perform BLAST search with the extended sequence
        # -----------------------------------------------------------------------------
        # Create a temporary FASTA file with the extended sequence
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
            temp_fasta.write(f">{index}_extended\n{seq}\n")
            temp_fasta_path: str = temp_fasta.name

        try:
            # Perform BLAST search with the extended sequence
            blast_results: pd.DataFrame = run_blastn_alignment(
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
            blast_results = blast_results[blast_results['len'] >= min_length].copy()
            
            # Clean up temporary file
            os.unlink(temp_fasta_path)

            # -----------------------------------------------------------------------------
            # STEP 11: Filter elements from the blast data that overlap our row data
            # -----------------------------------------------------------------------------
            # Create a DataFrame with the current element for checking
            current_element = pd.DataFrame([{
                'sstart': lower_coor_extended,  # type: int
                'send': upper_coor_extended,  # type: int
                'sseqid': element.sseqid,  # type: int
                'sstrand': element.sstrand  # type: int
            }])

            # Check which elements from `blast_results` overlap with `current_element`
            elems_to_remove: pd.DataFrame = get_interval_overlap(blast_results, current_element)

            # Remove `elems_to_remove` from `blast_results`
            blast_results = match_data_and_remove(blast_results, elems_to_remove)
            blast_results.reset_index(drop=True, inplace=True)

            # -----------------------------------------------------------------------------
            # STEP 12: Use next_side_extension_checker to determine further extension
            # -----------------------------------------------------------------------------
            if not blast_results.empty:
                # Check extension possibilities for the next recursive call:
                ## both: both sides can be extended
                ## left: only left side will be extended
                ## right: only right side can be extended
                ## None: neither side can be extended
                row_extended: pd.DataFrame
                row_extended, extension_direction = next_side_extension_checker(
                    current_element,
                    blast_results,
                    extension_direction
                )
                
                # -----------------------------------------------------------------------------
                # STEP 13: Recursive extension based on extension_direction
                # -----------------------------------------------------------------------------
                if extension_direction in ["both", "left", "right"]:
                    # Further extension is possible, prepare for recursive call
                    updated_element: pd.Series = pd.Series({
                        'sstart': row_extended['sstart'].iloc[0],
                        'send': row_extended['send'].iloc[0],
                        'sseqid': element['sseqid'],
                        'sstrand': element['sstrand']
                    })

                    # Recursive call with updated coordinates, direction, and previous extension tracking
                    # Record this step before recursing
                    recursive_result: Dict[str, Any] = _process_single_row_extension(
                        row_data=(index, updated_element),
                        genome_fasta=genome_fasta,
                        extend_number=extend_number,
                        limit_len=limit_len,
                        extension_direction=extension_direction,
                        identity=identity,
                        word_size=word_size,
                        min_length=min_length,
                        max_recursion_depth=max_recursion_depth,
                        current_depth=current_depth + 1
                    )
                    
                    # Return the result from the recursive call (it already has recursion info)
                    return recursive_result

                else:
                    # extension_direction is None - no further extension possible
                    # Return the coordinates from next_side_extension_checker
                    final_sstart: int = row_extended.iloc[0].sstart
                    final_send: int = row_extended.iloc[0].send
                    final_seq_len: int = final_send - final_sstart + 1

                    # Get the final sequence with updated coordinates
                    final_cmd: str = (
                        f"blastdbcmd -db {genome_fasta} "
                        f"-entry {element['sseqid']} "
                        f"-range {final_sstart}-{final_send} "
                        f"-strand {element['sstrand']} "
                        "-outfmt %s"
                    )
                    final_seq: str = subprocess.check_output(final_cmd, shell=True, universal_newlines=True).strip()
                    
                    result.update({
                        'modified': True,
                        'len': final_seq_len,
                        'sstart': final_sstart,
                        'send': final_send,
                        'sseq': final_seq,
                        'recursion_info': {
                            'max_depth': current_depth,
                            'status': 'extension_completed'
                        }
                    })
                    return result
            else:
                # No BLAST results after filtering - return current extension
                result.update({
                    'modified': True,
                    'len': new_subject_len,
                    'sstart': int(lower_coor_extended),
                    'send': int(upper_coor_extended),
                    'sseq': seq,
                    'recursion_info': {
                        'max_depth': current_depth,
                        'status': 'no_blast_results'
                    }
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
        final_cmd: str = (
            f"blastdbcmd -db {genome_fasta} "
            f"-entry {element['sseqid']} "
            f"-range {lower_coor}-{upper_coor} "
            f"-strand {element['sstrand']} "
            "-outfmt %s"
        )
        final_seq: str = subprocess.check_output(final_cmd, shell=True, universal_newlines=True).strip()

        result.update({
            'modified': True,
            'len': subject_len,
            'sstart': int(lower_coor),
            'send': int(upper_coor),
            'sseq': final_seq,
            'recursion_info': {
                'max_depth': current_depth,
                'status': 'length_requirement_met'
            }
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
    Extends sequences within a genome based on given conditions and resolves potential
    conflicts or overlaps between them.

    Processes input data by extending sequences in the specified genome, calculating
    recursion-based extensions if applicable. It then resolves conflicts among the
    processed sequences by maintaining priority and ensuring no overlaps, finally
    returning an updated dataset with adjusted sequences.

    Parameters
    ----------
    data_input : pandas.DataFrame
        Input DataFrame having data about sequences and their coordinates.
    genome_fasta : str
        Path to the FASTA file of the genome against which the sequences will be extended.
    extend_number : int
        Number of bases by which sequences should be extended.
    limit_len : int
        Maximum allowed length for the extended sequences.
    identity : int, optional
        Minimum percentage identity for successfully extending a sequence, by default 60.
    word_size : int, optional
        Size of the word/block used in the sequence analysis, by default 15.
    min_length : int, optional
        Minimum length a sequence must reach after extension, by default 100.
    extension_direction : str, optional
        Direction of extension, either 'both', 'upstream', or 'downstream', by default "both".
    max_recursion_depth : int, optional
        Maximum depth of recursion for extending sequences, by default 10.
    n_jobs : int, optional
        Number of parallel jobs for processing inputs, by default -1.

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame where sequences are extended and conflicts resolved.
    """

    # -----------------------------------------------------------------------------
    # STEP 1: Process all rows in parallel using joblib
    # -----------------------------------------------------------------------------
    # Distribute the processing of each row across multiple processors
    # Each row is processed independently by the _process_single_row_extension function
    results: List[Dict] = Parallel(n_jobs=n_jobs)(
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
    extended_count: int = 0
    recursion_stats: Dict[int, int] = {}
    
    for result in results:
        recursion_info: Dict[str, Any] = result.get('recursion_info', {})
        depth: int = recursion_info.get('max_depth', 0)
        status: str = recursion_info.get('status', 'unknown')
        
        if result['modified']:
            extended_count += 1
            if depth == 0:
                print(f"Extending row {result['index']} (Extension completed - no recursion made)")
            else:
                print(f"Extending row {result['index']} ({depth} recursive calls - {status})")
            
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
    # STEP 3: Resolve conflicts across sequences
    # -----------------------------------------------------------------------------
    # Priority is by original index (A has priority over B if A's index < B's index)
    results_by_index: Dict[int, Dict] = {r['index']: r for r in results}
    ordered_indices: List[int] = sorted(results_by_index.keys())

    accepted_df: pd.DataFrame = pd.DataFrame(columns=['sseqid', 'sstart', 'send', 'sstrand'])

    removed: int = 0
    removed_idx: List[int] = []
    total: int = len(ordered_indices)

    for idx in ordered_indices:
        row: pd.Series = data_input.loc[idx]
        res: Dict[str, Any] = results_by_index[idx]

        # Determine the single final coordinates for this row
        final_sstart: int
        final_send: int
        final_seq_len: int
        final_seq: str|None
        if res.get('modified'):
            final_sstart = res.get('sstart')
            final_send = res.get('send')
            final_seq_len = res.get('len')
            final_seq = res.get('sseq')
        else:
            # Unmodified: keep original coordinates
            final_sstart = row.sstart
            final_send = row.send
            final_seq_len = row.len
            final_seq = None

        candidate: pd.DataFrame = pd.DataFrame([{
            'sseqid': row.sseqid,  # type: str
            'sstart': final_sstart,  # type: int
            'send': final_send,  # type: int
            'sstrand': row.sstrand,  # type: str
            'len': final_seq_len,  # type: int
            'sseq': final_seq  # type: str|None
        }])

        # If conflicts with already accepted intervals, remove B entirely
        if not accepted_df.empty and not get_interval_overlap(candidate, accepted_df).empty:
            removed += 1
            removed_idx.append(idx)
            continue

        # Accept this candidate
        accepted_df = pd.concat([accepted_df, candidate], ignore_index=True)
        
    # Now, get sequence for each element that has None in `accepted_df['sseq']`
    none_mask: pd.Series[bool] = accepted_df['sseq'].isna()
    
    if none_mask.any():
        print(f"\nRetrieving sequences for {none_mask.sum()} elements with missing sequences...")
        
        for idx in accepted_df[none_mask].index:
            row: pd.Series = accepted_df.loc[idx]
            
            # Get the sequence using blastdbcmd
            cmd: str = (
                f"blastdbcmd -db {genome_fasta} "
                f"-entry {row['sseqid']} "
                f"-range {row['sstart']}-{row['send']} "
                f"-strand {row['sstrand']} "
                "-outfmt %s"
            )
            
            seq: str = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()
            accepted_df.loc[idx, 'sseq'] = seq

    # Print conflict resolution summary
    print(f"\nConflict resolution summary:")
    print(f"  - Total sequences considered: {total}")
    print(f"  - Sequences removed due to conflicts: {removed}")

    # Return the updated DataFrame with resolved sequences
    return accepted_df
