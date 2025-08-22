"""
SIDER Functions Module

This module contains utility functions for SIDER (Short Interspersed Degenerated Retroposons) filtering operations.
These functions support the main SIDER filtering pipeline by providing functionality for:
- Command-line argument parsing
- Sequence processing with BLASTN
- Directory and database setup
- Sequence filtering based on SIDER criteria
- Processing of recaught sequences
- Saving results to disk

This module is designed to be imported by the main sider_filter.py script.
"""

import argparse
import os
import sys
import logging
import pandas as pd
import subprocess
from typing import Dict, Tuple
from joblib import Parallel, delayed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from blastoise.modules.blaster import create_blast_database
from .extra_functions import general_blastn_blaster




def process_sequence(
        row: pd.Series,
        blastn_db_path: str,
        word_size: int,
        evalue: float,
        min_subjects: int = 5
) -> Dict[str, str]:
    """
    BLAST a single nucleotide sequence against a reference database and decide
    whether it meets the SIDER acceptance threshold.

    Parameters
    ----------
    row : pd.Series
        Row containing at least the columns ``name_id`` and ``sequence``.
    blastn_db_path : str
        Path to the BLAST-formatted nucleotide database.
    word_size : int
        BLASTN word size.
    evalue : float
        BLASTN E-value cutoff.
    min_subjects : int, default 5
        Minimum number of distinct subject hits required for acceptance.

    Returns
    -------
    Dict[str, str]
        ``{"name_id": <str>, "status": "Accepted" | "Rejected"}``.
    """
    # Use the sequence directly from the CSV
    if 'seq' not in row:
        logging.getLogger('sider_filter').warning(f"No sequence found for {row['name_id']}")
        return {'name_id': row['name_id'], 'status': 'Rejected'}

    sequence = row['seq']
    name_id = row['name_id']

    # Create a temporary query file for BLASTN
    query = f"<(echo -e '>{name_id}\\n{sequence}')"

    try:
        # Run BLASTN
        blastn_df = general_blastn_blaster(
            query_path=query,
            dict_path=blastn_db_path,
            word_size=word_size,
            evalue=evalue
        )

        # Check BLASTN results - a sequence is accepted if it has hits to at least min_subjects different subjects
        if not blastn_df.empty and blastn_df['sseqid'].nunique() >= min_subjects:
            return {'name_id': name_id, 'status': 'Accepted'}
        else:
            return {'name_id': name_id, 'status': 'Rejected'}
    except Exception as e:
        logging.getLogger('sider_filter').error(f"Error processing sequence {name_id}: {str(e)}")
        return {'name_id': name_id, 'status': 'Rejected'}




def filter_sequences(
        data: pd.DataFrame,
        blastn_db_path: str,
        word_size: int,
        evalue: float,
        min_subjects: int = 5
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter sequences using SIDER criteria by applying BLASTN-based filtering.

    A sequence is accepted if it has BLASTN hits to at least `min_subjects` different subjects,
    in the expected minimum `evalue`, otherwise it is rejected. The function processes sequences
    in parallel using joblib.

    For each sequence:
    1. Creates a unique name_id from sequence coordinates
    2. Runs BLASTN against the provided database
    3. Checks if the number of unique subject hits meets minimum threshold
    4. Labels sequence as 'Accepted' or 'Rejected'

    Parameters:
    -----------
    data: pd.DataFrame
        A data frame containing sequence data with columns for sseqid, sstrand, sstart, send,
        and sseq (the sequence itself)
    blastn_db_path: str
        Path to the pre-formatted BLASTN database to search against
    word_size: int
        Word size parameter for BLASTN search (-W parameter)
    evalue: float
        E-value threshold for BLASTN hits (-e parameter)
    min_subjects: int
        Minimum number of unique subjects a sequence must hit to be accepted (default: 5)

    Returns:
    --------
    Tuple[pd.DataFrame, pd.DataFrame]
        A tuple containing:
        - First DataFrame with accepted sequences that meet the SIDER criteria
        - Second DataFrame with rejected sequences that don't meet the criteria
    """
    logger = logging.getLogger('sider_filter')
    logger.info("Preparing data for analysis")

    # Create a unique identifier for each sequence
    data['name_id'] = data.apply(
        lambda row: f"{row['chromosome']}_{row['strand']}_{row['start']}-{row['end']}",
        axis=1
    )

    logger.info("Applying SIDER filter")
    logger.info(f"Analyzing {len(data)} sequences...")

    # Process sequences in parallel
    try:
        results = Parallel(n_jobs=-1)(
            delayed(process_sequence)(row, blastn_db_path, word_size, evalue, min_subjects) 
            for _, row in data.iterrows()
        )

        # Convert results to DataFrame and merge with the original data
        results_df = pd.DataFrame(results)
        data = pd.merge(data, results_df, on='name_id', how='left')

        # Split into accepted and rejected
        accepted_data = data[data['status'] == 'Accepted'].copy()
        rejected_data = data[data['status'] == 'Rejected'].copy()

        logger.info(f"Accepted sequences: {len(accepted_data)}")
        logger.info(f"Rejected sequences: {len(rejected_data)}")

        return accepted_data, rejected_data

    except Exception as e:
        logger.error(f"Error during sequence filtering: {str(e)}")
        raise


def process_recaught_data(
    rejected_data: pd.DataFrame, 
    temp_dir: str, 
    blastn_db_path: str, 
    recaught_file_path: str,
    identity: int,
    word_size: int,
    recaught_threshold: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process rejected sequences to recapture potentially valid sequences.

    Parameters:
    -----------
    rejected_data: pd.DataFrame
        A DataFrame containing rejected sequences
    temp_dir: str
        Path to the temporary directory
    blastn_db_path: str
        Path to the BLASTN database
    recaught_file_path: str
        Path to the recaught file
    identity: Minimum identity percentage
    word_size: int
        Word size parameter for BLASTN
    recaught_threshold: int
        E-value threshold for recapturing

    Returns:
    --------
    Tuple[pd.DataFrame, pd.DataFrame]
        Updated DataFrames with accepted and rejected sequences
    """
    logger = logging.getLogger('sider_filter')
    logger.info("Processing recaught data")

    # If there are no rejected sequences, return empty DataFrames
    if rejected_data.empty:
        logger.info("No rejected sequences to process")
        return pd.DataFrame(), rejected_data

    try:
        # Reset index to ensure continuous indexing for iloc operations
        rejected_data_reset = rejected_data.reset_index(drop=True)

        # Create a FASTA file from the rejected data for BLASTN using BioPython
        fasta_file_path = os.path.join(temp_dir, "negative_database.fasta")

        # Create a list of SeqRecord objects
        records = []
        for idx, row in rejected_data_reset.iterrows():
            if 'seq' in row:
                # Create a SeqRecord with the sequence and an ID
                record = SeqRecord(
                    Seq(row['seq']),
                    id=f"Seq_{idx}_{row['chromosome']}",
                    description=""
                )
                records.append(record)
            else:
                logger.warning(f"No sequence found for rejected row {idx}")

        # Write the records to a FASTA file
        SeqIO.write(records, fasta_file_path, "fasta")

        # Create BLASTN database
        create_blast_database(path_input=fasta_file_path, path_output=blastn_db_path)

        # Run BLASTN on a recaught file
        caught_data = general_blastn_blaster(
            query_path=recaught_file_path,
            dict_path=blastn_db_path,
            perc_identity=identity,
            word_size=word_size,
            evalue=recaught_threshold
        )

        # Process recaught data
        if not caught_data.empty:
            # Filter by e-value
            caught_data = caught_data.sort_values(by=["evalue"])
            logger.info(f"Recaught data: {caught_data.shape[0]} elements")

            # Extract sequence indices from sseqid
            caught_data["index"] = caught_data["sseqid"].str.extract(r"Seq_(\d+)_").astype(int)

            # Get indices of recaught sequences
            recaught_indices = caught_data["index"].unique()

            # Validate indices are within bounds
            valid_indices = [idx for idx in recaught_indices if idx < len(rejected_data_reset)]

            if len(valid_indices) != len(recaught_indices):
                logger.warning(
                    f"Some recaught indices were out of bounds. Using {len(valid_indices)} out of {len(recaught_indices)} indices.")

            # Move recaught sequences from rejected to accept
            if valid_indices:
                recaught_rows = rejected_data_reset.iloc[valid_indices].copy()
                recaught_rows['status'] = 'Accepted (Recaught)'
                accepted_data = recaught_rows

                # Remove recaught sequences from rejected data
                remaining_rejected = rejected_data_reset.drop(rejected_data_reset.index[valid_indices])

                logger.info(f"Recaught sequences: {len(recaught_rows)}")
                logger.info(f"Updated rejected sequences: {len(remaining_rejected)}")

                return accepted_data, remaining_rejected

        logger.info("No sequences were recaptured")
        return pd.DataFrame(), rejected_data

    except Exception as e:
        logger.error(f"Error during recaught data processing: {str(e)}")
        return pd.DataFrame(), rejected_data


def save_results(
    accepted_data: pd.DataFrame, 
    rejected_data: pd.DataFrame, 
    output_dir: str, 
    temp_dir: str
) -> Tuple[str, str]:
    """
    Persist accepted/rejected sequences to disk and return the resulting file paths.

    Parameters
    ----------
    accepted_data : pd.DataFrame
        Sequences that passed filtering.
    rejected_data : pd.DataFrame
        Sequences that failed filtering.
    output_dir : str
        Destination directory for the final files.
    temp_dir : str
        Working directory for intermediate artifacts.

    Returns
    -------
    Tuple[str, str]
        (positive_file_path, negative_file_path)

    """
    logger = logging.getLogger('sider_filter')
    logger.info("Saving results")

    # Define output paths
    positive_path = os.path.join(output_dir, "siders_df.csv")
    negative_path = os.path.join(output_dir, "non_siders_df.csv")

    try:
        # Save DataFrames to CSV
        accepted_data.to_csv(positive_path, index=False)
        rejected_data.to_csv(negative_path, index=False)

        logger.info(f"Positive results saved to: {positive_path}")
        logger.info(f"Negative results saved to: {negative_path}")

        # Set file permissions
        subprocess.run(["chmod", "-R", "a+w", temp_dir], check=True)
        subprocess.run(["chmod", "a+w", positive_path], check=True)
        subprocess.run(["chmod", "a+w", negative_path], check=True)

        return positive_path, negative_path

    except Exception as e:
        logger.error(f"Error saving results: {str(e)}")
        raise
