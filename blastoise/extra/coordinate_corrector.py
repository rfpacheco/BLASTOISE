"""
Coordinate Corrector Script

This script adjusts sequence coordinates using BLASTN and BEDOPS processing.
It processes input sequences, filters them based on coordinate overlap and strand validity,
merges the results with BEDOPS, and computes updated start and end coordinates.

The script requires several input parameters and produces a CSV file with the corrected coordinates.
"""

import argparse
import os
import sys
import logging
import pandas as pd
import subprocess
from typing import Dict, List, Tuple, Any

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from modules.blaster import blastn_dic
from extra.second_functions import get_sequence, bedops_merge, general_blastn_blaster
from modules.aesthetics import print_message_box

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('coordinate_corrector')


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the coordinate corrector script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
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
        "-o", "--output", 
        type=str, 
        required=True,
        help="Path to the output CSV file where corrected coordinates will be saved."
    )

    parser.add_argument(
        "-ws", "--word_size", 
        type=int, 
        required=True,
        help="Word size parameter for BLASTN."
    )

    parser.add_argument(
        "-min", "--min_length", 
        type=int, 
        required=True,
        help="Minimum length of the sequence to be considered after correction."
    )

    return parser.parse_args()


def setup_directories(input_file_dir: str, dict_path: str) -> Tuple[str, str]:
    """
    Create necessary directories and prepare BLASTN database.

    Args:
        input_file_dir: Directory of the input file
        dict_path: Path to the genome FASTA file

    Returns:
        Tuple containing paths to the temporary directory and BLASTN database
    """
    logger.info("Setting up directories and BLASTN database")

    # Prepare a subfolder for temporary files
    temp_dir = os.path.join(input_file_dir, "tmpCoordinateCorrector")
    os.makedirs(temp_dir, exist_ok=True)

    # Prepare BLASTn dict
    dict_folder_path = os.path.join(temp_dir, "blastn_dict")
    os.makedirs(dict_folder_path, exist_ok=True)

    blastn_db_path = os.path.join(dict_folder_path, os.path.basename(dict_path))

    try:
        blastn_dic(path_input=dict_path, path_output=blastn_db_path)
        logger.info(f"BLASTN database created at {blastn_db_path}")
    except Exception as e:
        logger.error(f"Error creating BLASTN database: {str(e)}")
        raise

    return temp_dir, blastn_db_path


def correct_coordinates(
    df: pd.DataFrame, 
    blastn_db_path: str, 
    temp_dir: str, 
    word_size: int, 
    min_length: int
) -> Dict[str, List[List[Any]]]:
    """
    Adjust sequence coordinates using BLASTN and BEDOPS processing.

    This function takes a DataFrame containing sequence information, uses BLASTN to process sequences,
    filters sequences based on coordinate overlap and strand validity, merges the results with BEDOPS,
    and computes updated start and end coordinates.

    Args:
        df: DataFrame containing sequence data
        blastn_db_path: Path to the BLASTN database
        temp_dir: Path to the temporary directory
        word_size: Word size parameter for BLASTN
        min_length: Minimum length of the sequence to be considered

    Returns:
        Dictionary containing the corrected coordinates
    """
    logger.info("Starting coordinate correction process")

    # Create a dictionary to store results
    results_dict = {}

    # Process each row in the DataFrame
    for pos, (index, row) in enumerate(df.iterrows(), start=1):
        # Prepare the query data
        name_id = f"{row['sseqid']}_{row['sstrand']}_{row['sstart']}-{row['send']}"
        seq = row["sseq"]
        query = f"<(echo -e '>{name_id}\\n{seq}')"
        start_coor = row["sstart"]
        end_coor = row["send"]
        strand_seq = row["sstrand"]
        name_chr = row["sseqid"]

        logger.info(f"Analyzing row {pos}/{df.shape[0]} with name_id {name_id}")

        # Run BLASTN
        blastn_df = general_blastn_blaster(
            query_path=query, 
            dict_path=blastn_db_path, 
            word_size=word_size
        )

        # Filter the BLASTN results
        # Remove rows with coordinates that overlap with the original coordinates
        blastn_df = filter_blastn_results(
            blastn_df, 
            start_coor, 
            end_coor, 
            name_chr
        )

        # Sort and merge the filtered results
        blastn_df.sort_values(by=["qstart"], inplace=True)
        bedops_df = bedops_merge(input_df=blastn_df, path_folder=temp_dir)

        # Initialize the results for this sequence
        results_dict[name_id] = []

        # Calculate new coordinates
        for _, bedops_row in bedops_df.iterrows():
            # Calculate new start and end coordinates
            new_start = start_coor + bedops_row["qstart"] - 1
            new_end = start_coor + bedops_row["qend"] - 1

            # Check if the sequence meets the minimum length requirement
            if abs(new_end - new_start) + 1 >= min_length:
                results_dict[name_id].append([name_chr, strand_seq, new_start, new_end])

        # Log the results
        logger.info(f"Found {len(results_dict[name_id])} new sequences for {name_id}")
        for seq in results_dict[name_id]:
            logger.info(f"  {seq[0]}:{seq[1]}:{seq[2]}-{seq[3]}")

    return results_dict


def filter_blastn_results(
    blastn_df: pd.DataFrame, 
    start_coor: int, 
    end_coor: int, 
    name_chr: str
) -> pd.DataFrame:
    """
    Filter BLASTN results to remove rows with coordinates that overlap with the original coordinates.

    Args:
        blastn_df: DataFrame containing BLASTN results
        start_coor: Start coordinate of the original sequence
        end_coor: End coordinate of the original sequence
        name_chr: Chromosome name of the original sequence

    Returns:
        Filtered DataFrame
    """
    # Filter for plus strand
    # Remove rows where:
    # (sstart is within start_coor and end_coor OR send is within start_coor and end_coor)
    # AND sseqid matches name_chr AND sstrand is "plus"
    plus_filter = (
        (((blastn_df["sstart"] >= start_coor) & (blastn_df["sstart"] <= end_coor)) |
         ((blastn_df["send"] <= end_coor) & (blastn_df["send"] >= start_coor))) &
        (blastn_df["sseqid"] == name_chr) &
        (blastn_df["sstrand"] == "plus")
    )
    blastn_df = blastn_df[~plus_filter].copy()

    # Filter for minus strand
    # Remove rows where:
    # (sstart is within start_coor and end_coor OR send is within start_coor and end_coor)
    # AND sseqid matches name_chr AND sstrand is "minus"
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

    Args:
        results_dict: Dictionary containing corrected coordinates

    Returns:
        DataFrame containing corrected coordinates
    """
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

    logger.info(f"Created DataFrame with {len(result_df)} rows")

    return result_df


def add_sequences(df: pd.DataFrame, blastn_db_path: str) -> pd.DataFrame:
    """
    Add sequence data to the DataFrame.

    Args:
        df: DataFrame containing coordinate data
        blastn_db_path: Path to the BLASTN database

    Returns:
        DataFrame with added sequence data
    """
    if df.empty:
        logger.warning("DataFrame is empty, no sequences to add")
        return df

    logger.info("Adding sequences to DataFrame")

    # Add sequence data to each row
    df['sseq'] = df.apply(
        lambda x: get_sequence(
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

    Args:
        df: DataFrame to save
        output_path: Path to the output CSV file
        temp_dir: Path to the temporary directory
    """
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


def main():
    """Main function to run the coordinate correction pipeline."""
    try:
        # Parse arguments
        args = parse_arguments()

        # Expand user paths
        csv_path = os.path.expanduser(args.file)
        dict_path = os.path.expanduser(args.dict_path)
        output_path = os.path.expanduser(args.output)

        # Get the directory of the input file
        input_file_dir = os.path.dirname(csv_path)

        # Setup directories and BLASTN database
        temp_dir, blastn_db_path = setup_directories(input_file_dir, dict_path)

        # Read input CSV
        data = pd.read_csv(csv_path, sep=",", header=0)
        logger.info(f"Loaded {len(data)} sequences from {csv_path}")

        # Print message
        print_message_box("Correcting coordinates")

        # Correct coordinates
        results_dict = correct_coordinates(
            data, 
            blastn_db_path, 
            temp_dir, 
            args.word_size, 
            args.min_length
        )

        # Create output DataFrame
        result_df = create_output_dataframe(results_dict)

        # Add sequences to DataFrame
        result_df = add_sequences(result_df, blastn_db_path)

        # Save results
        save_results(result_df, output_path, temp_dir)

    except Exception as e:
        logger.error(f"Error in coordinate correction pipeline: {str(e)}")
        print_message_box(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
    print_message_box("Coordinate correction completed")
