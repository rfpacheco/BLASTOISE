import argparse
import os
import sys
import pandas as pd
import subprocess
from typing import Dict, List, Optional, Union, Any

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from modules.blaster import blastn_dic
from extra.main_functions import coordinates_corrector
from modules.aesthetics import print_message_box


def parse_arguments() -> argparse.Namespace:
    """
    Parses command-line arguments for the coordinate corrector script.

    This function sets up an argument parser with options for input and
    output file paths as well as parameters controlling the behavior of the
    coordinate correction process.

    Returns:
    --------
        argparse.Namespace: An object containing parsed command-line arguments.

    Arguments:
    ----------
        -f, --file (str): Path to the input CSV file (required).
        -d, --dict_path (str): Path to the genome fasta file (required).
        -o, --output (str): Path to the output CSV file (required).
        -ws, --word_size (int): word_size value for BLASTN (required).
        -min, --min_length (int): Minimum length of the sequence to be considered
            (required).
    """
    parser = argparse.ArgumentParser(description="Correct coordinates for sequences in a CSV file")
    parser.add_argument("-f", "--file", type=str, required=True,
                        help="Path to the input CSV file.")
    parser.add_argument("-d", "--dict_path", type=str, required=True,
                        help="Path to the genome fasta file.")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to the output CSV file.")
    parser.add_argument("-ws", "--word_size", type=int, required=True,
                        help="word_size value of BLASTN")
    parser.add_argument("-min", "--min_length", type=int, required=True,
                        help="Minimum length of the sequence to be considered.")
    return parser.parse_args()


if __name__ == "__main__":
    # Load all arguments
    args = parse_arguments()
    csv_path = os.path.expanduser(args.file)
    dict_path = os.path.expanduser(args.dict_path)
    output_path = os.path.expanduser(args.output)
    word_size = args.word_size
    min_length = args.min_length

    # Prepare subfolder
    csv_parent_path = os.path.dirname(csv_path)
    folder_path = os.path.join(csv_parent_path, "tmpCoordinateCorrector")  # tmp folder to place all data
    os.makedirs(folder_path, exist_ok=True)

    # Prepare BLASTn dict
    dict_folder_path = os.path.join(folder_path, "blastn_dict")
    os.makedirs(dict_folder_path, exist_ok=True)
    dict_file_path_out = os.path.join(dict_folder_path, os.path.basename(dict_path))
    blastn_dic(path_input=dict_path, path_output=dict_file_path_out)

    # Read input CSV
    data = pd.read_csv(csv_path, sep=",", header=0)
    print_message_box("Correcting coordinates")
    
    # Perform coordinate correction
    json_path = coordinates_corrector(
        df=data,
        dict_path=dict_file_path_out,
        folder_path=folder_path,
        word_size=word_size,
        min_length=min_length
    )
    
    # Convert JSON to CSV
    with open(json_path, 'r') as f:
        import json
        json_data = json.load(f)
    
    # Extract data from JSON and create a DataFrame
    rows = []
    for key, value in json_data.items():
        for element in value:
            rows.append({
                'sseqid': element[0],
                'sstrand': element[1],
                'sstart': element[2],
                'send': element[3]
            })
    
    # Create DataFrame and save to CSV
    result_df = pd.DataFrame(rows)
    if not result_df.empty:
        # Get sequences for each row
        from extra.second_functions import get_sequence
        result_df['sseq'] = result_df.apply(
            lambda x: get_sequence(
                start_coor=x['sstart'],
                end_coor=x['send'],
                strand=x['sstrand'],
                chromosome=x['sseqid'],
                path_genome=dict_file_path_out
            ), axis=1
        )
        
        # Save to CSV
        result_df.to_csv(output_path, index=False)
        print_message_box(f"Coordinate correction completed. Results saved to {output_path}")
    else:
        print_message_box("No sequences found after coordinate correction.")
    
    # Add right to groups and users
    try:
        subprocess.run(["chmod", "-R", "a+w", folder_path], check=True)
        subprocess.run(["chmod", "a+w", output_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while changing permissions: {e}")

print_message_box("Coordinate correction completed")