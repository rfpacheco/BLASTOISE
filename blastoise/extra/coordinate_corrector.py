"""
Coordinate Corrector Script

This script adjusts sequence coordinates using BLASTN and PyRanges processing.
It processes input sequences, filters them based on coordinate overlap and strand validity,
merges the results with PyRanges, and computes updated start and end coordinates.

The script requires several input parameters and produces a CSV file with the corrected coordinates.
"""

import os
import sys
import logging
import pandas as pd

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from blastoise.modules.aesthetics import print_message_box
from .utils.extra_functions import setup_directories
from .utils.coordinate_corrector_functions import (
    parse_arguments,
    correct_coordinates,
    create_output_dataframe,
    add_sequences,
    save_results
)


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('coordinate_corrector')


def main():
    """
    Main function to run the coordinate correction pipeline.

    This function orchestrates the entire coordinate correction process by:
    1. Parsing command-line arguments
    2. Setting up directories and BLASTN database
    3. Reading input data
    4. Correcting coordinates
    5. Creating output DataFrame
    6. Adding sequences to DataFrame
    7. Saving results to disk

    Raises
    ------
    Exception
        If any error occurs during the pipeline execution.
    """
    try:
        # Parse arguments
        args = parse_arguments()

        # Expand user paths
        csv_path = os.path.expanduser(args.file)
        dict_path = os.path.expanduser(args.dict_path)

        # User as an output path to the folder where the csv_path is
        ## File name as the original file without the extension + corrected.csv
        out_file_name = os.path.splitext(os.path.basename(csv_path))[0] + "_corrected.csv"
        output_path = os.path.join(os.path.dirname(csv_path), out_file_name)
        output_path = str(output_path)

        # Get the directory of the input file
        input_file_dir = os.path.dirname(csv_path)

        # Setup directories and BLASTN database
        temp_dir, blastn_db_path = setup_directories(
            base_dir=input_file_dir,
            dict_path=dict_path,
            temp_dir_name="tmpCoordinateCorrector",
            logger_name="coordinate_corrector"
        )

        # Read input CSV
        data = pd.read_csv(csv_path, sep=",", header=0)
        logger.info(f"Loaded {len(data)} sequences from {csv_path}")

        # Print message
        print_message_box("Correcting coordinates")

        # Correct coordinates
        results_dict = correct_coordinates(
            data, 
            blastn_db_path,
            args.word_size, 
            args.min_length,
            args.max_length,
            args.identity,
            args.evalue
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
