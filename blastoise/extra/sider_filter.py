"""
SIDER Filter Script

This script filters sequences using SIDER (Short Interspersed Degenerated Retroposons) criteria.
It processes input sequences, applies BLASTN-based filtering, and categorizes them as accepted or rejected.
Additionally, it attempts to recapture sequences that might have been incorrectly rejected.

The script requires several input parameters and produces two output CSV files in the same directory as the input file:
- siders_df.csv: Contains sequences that meet the SIDER criteria
- non_siders_df.csv: Contains sequences that do not meet the SIDER criteria
"""

import os
import sys
import logging
import pandas as pd

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from modules.aesthetics import print_message_box
from extra.sider_filter_functions import (
    parse_arguments,
    setup_directories,
    filter_sequences,
    process_recaught_data,
    save_results
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('sider_filter')


def main():
    """Main function to run the SIDER filter pipeline."""
    try:
        # Parse arguments
        args = parse_arguments()

        # Expand user paths
        csv_path = os.path.expanduser(args.file)
        dict_path = os.path.expanduser(args.dict_path)
        # Use the directory of the input file as the output location
        output_dir = os.path.dirname(csv_path)
        recaught_file_path = os.path.expanduser(args.recaught_file)

        # Setup directories and BLASTN database
        temp_dir, blastn_db_path = setup_directories(output_dir, dict_path)

        # Read input CSV
        data = pd.read_csv(csv_path, sep=",", header=0)
        logger.info(f"Loaded {len(data)} sequences from {csv_path}")

        # Filter sequences
        print_message_box("Applying SIDER filter")
        accepted_data, rejected_data = filter_sequences(
            data, 
            blastn_db_path, 
            args.word_size, 
            args.evalue,
            args.min_subjects
        )

        # Process recaught data
        print_message_box("Processing recaught data")
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
            temp_dir
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


if __name__ == "__main__":
    main()
    print_message_box("SIDER filtering completed")
