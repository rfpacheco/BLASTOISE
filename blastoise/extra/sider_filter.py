"""
SIDER Filter Script

This script filters sequences using SIDER (Short Interspersed Degenerated Retroposons) criteria.
It processes input sequences, applies BLASTN-based filtering, and categorizes them as accepted or rejected.
Additionally, it attempts to recapture sequences that might have been incorrectly rejected.

Outputs:
- Two CSV files saved in the same directory as the input file, named using the input filename's basename:
  - <input_basename>--sider.csv: Sequences that meet the SIDER criteria
  - <input_basename>--non_sider.csv: Sequences that do not meet the SIDER criteria
"""

import os
import sys
import logging
import pandas as pd
import argparse
import shutil

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Create logger instance
logger = logging.getLogger(__name__)

from blastoise.modules.aesthetics import print_message_box
from .utils.extra_functions import setup_directories
from .utils.sider_filter_functions import (
    filter_sequences,
    process_recaught_data,
    save_results
)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for the SIDER filter script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Filter sequences using SIDER criteria",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-d", "--data",
        type=str, 
        required=True,
        help="Path to the input CSV file containing sequence data."
    )

    parser.add_argument(
        "-g", "--genome",
        type=str, 
        required=True,
        help="Path to the genome FASTA file for BLASTN database creation."
    )


    parser.add_argument(
        "-rf", "--recaught_file", 
        type=str, 
        required=True,
        help="Path to the FASTA file used for recapturing sequences."
    )

    parser.add_argument(
        "-rt", "--recaught_threshold", 
        type=float,
        default=1.0E-03,
        help="E-value threshold for recapturing sequences."
    )

    parser.add_argument(
        "-ws", "--word_size", 
        type=int, 
        default=15,
        help="Word size parameter for BLASTN."
    )

    parser.add_argument(
        "-e", "--evalue", 
        type=float, 
        default=1.0E-09,
        help="E-value threshold for initial BLASTN filtering."
    )

    parser.add_argument(
        "-i", "--identity", 
        type=int, 
        default=60,
        help="Minimum percentage identity for sequence recapturing."
    )

    parser.add_argument(
        "-ms", "--min_subjects", 
        type=int, 
        default=5,
        help="Minimum number of unique subjects required for a sequence to be accepted."
    )

    parser.add_argument(
        '-j', '--jobs',
        type=int,
        default=-1,
        help='Number of jobs for parallel processing. -1 means using all processors.'
    )

    return parser.parse_args()


def main():
    """Main function to run the SIDER filter pipeline."""
    temp_dir = None
    try:
        # Parse arguments
        args: argparse.Namespace = parse_arguments()

        # Expand user paths
        csv_path: str = os.path.expanduser(args.data)
        dict_path: str = os.path.expanduser(args.genome)
        # Use the directory of the input file as the output location
        output_dir: str = os.path.dirname(csv_path)
        recaught_file_path: str = os.path.expanduser(args.recaught_file)

        # Setup directories and BLASTN database
        temp_dir: str; blastn_db_path: str
        temp_dir, blastn_db_path = setup_directories(
            base_dir=output_dir,
            dict_path=dict_path,
            temp_dir_name="tmpSiderFilter",
        )

        # Read input CSV
        data: pd.DataFrame = pd.read_csv(csv_path, sep=",", header=0)
        print(f"Loaded {len(data)} sequences from {csv_path}")

        # Filter sequences
        print_message_box("Applying SIDER filter")
        accepted_data: pd.DataFrame; rejected_data: pd.DataFrame
        accepted_data, rejected_data = filter_sequences(
            data, 
            blastn_db_path, 
            args.word_size, 
            args.evalue,
            args.min_subjects,
            args.jobs
        )

        # Process recaught data
        print_message_box("Processing recaught data")
        recaught_accepted: pd.DataFrame; updated_rejected: pd.DataFrame
        recaught_accepted, updated_rejected = process_recaught_data(
            rejected_data,
            temp_dir,
            blastn_db_path,
            recaught_file_path,
            args.identity,
            args.word_size,
            args.recaught_threshold
        )

        # Combine original accepted and recaught accepted data
        if not recaught_accepted.empty:
            final_accepted = pd.concat([accepted_data, recaught_accepted], ignore_index=True)
            logger.info(f"Final accepted sequences: {len(final_accepted)}")
        else:
            final_accepted = accepted_data

        # Save results
        positive_path, negative_path = save_results(
            final_accepted,
            updated_rejected,
            output_dir,
            temp_dir,
            csv_path,
        )

        # Calculate the number of recaught elements
        num_recaught = len(recaught_accepted) if not recaught_accepted.empty else 0

        # Print counts to STDOUT
        print(f"Number of elements in accepted data: {len(final_accepted)}")
        print(f"Number of elements in rejected data: {len(updated_rejected)}")
        print(f"Number of elements recaught: {num_recaught}")

        # Print completion message
        message_to_print = (
            f"SIDER filtering completed. Results saved to:\n"
            f"- Positive: {positive_path}\n"
            f"- Negative: {negative_path}"
        )
        print_message_box(message_to_print)

    except Exception as e:
        logger.error(f"Error in SIDER filter pipeline: {str(e)}")
        print_message_box(f"Error: {str(e)}")
        sys.exit(1)
    finally:
        # Clean up temporary directory to avoid residue files
        try:
            if temp_dir and os.path.isdir(temp_dir):
                shutil.rmtree(temp_dir)
                logger.info(f"Temporary directory removed: {temp_dir}")
        except Exception as cleanup_err:
            logger.warning(f"Failed to remove temporary directory {temp_dir}: {cleanup_err}")


if __name__ == "__main__":
    main()
    print_message_box("SIDER filtering completed")
