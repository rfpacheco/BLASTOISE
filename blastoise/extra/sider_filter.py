import argparse
import os
import sys
import pandas as pd
import subprocess
import json

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from modules.blaster import blastn_dic
from extra.main_functions import json_sider_filter, sider_json_to_csv
from modules.aesthetics import print_message_box


def parse_arguments() -> argparse.Namespace:
    """
    Parses command-line arguments for the SIDER filter script.

    This function sets up an argument parser with options for input and
    output file paths as well as parameters controlling the behavior of the
    SIDER filtering process.

    Arguments:
    ----------
        -f, --file (str): Path to the input CSV file (required).
        -d, --dict_path (str): Path to the genome fasta file (required).
        -o, --output (str): Path to the output directory for positive and negative CSV files (required).
        -rf, --recaught_file (str): Path to the recaught file (required).
        -rt, --recaught_threshold (float): Recaught threshold (default: 1.0E-03).
        -ws, --word_size (int): word_size value for BLASTN (required).
        -e, --evalue (float): E-value threshold for BLASTN (default: 1.0E-09).
        -i, --identity (int): Minimum identity of the sequence to be considered (required).
    """
    parser = argparse.ArgumentParser(description="Filter sequences using SIDER criteria")
    parser.add_argument("-f", "--file", type=str, required=True,
                        help="Path to the input CSV file.")
    parser.add_argument("-d", "--dict_path", type=str, required=True,
                        help="Path to the genome fasta file.")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to the output directory for positive and negative CSV files.")
    parser.add_argument("-rf", "--recaught_file", type=str, required=True,
                        help="Path to the recaught file.")
    parser.add_argument("-rt", "--recaught_threshold", type=float, default=1.0E-03,
                        help="Recaught threshold.")
    parser.add_argument("-ws", "--word_size", type=int, required=True,
                        help="word_size value of BLASTN")
    parser.add_argument("-e", "--evalue", type=float, default=1.0E-09,
                        help="E-value threshold for BLASTN.")
    parser.add_argument("-i", "--identity", type=int, required=True,
                        help="Minimum identity of the sequence to be considered.")
    return parser.parse_args()


if __name__ == "__main__":
    # Load all arguments
    args = parse_arguments()
    csv_path = os.path.expanduser(args.file)
    dict_path = os.path.expanduser(args.dict_path)
    output_dir = os.path.expanduser(args.output)
    recaught_file_path = os.path.expanduser(args.recaught_file)
    recaught_threshold = args.recaught_threshold
    word_size = args.word_size
    evalue = args.evalue
    identity = args.identity

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Prepare a subfolder for temporary files
    folder_path = os.path.join(output_dir, "tmpSiderFilter")
    os.makedirs(folder_path, exist_ok=True)

    # Prepare BLASTn dict
    dict_folder_path = os.path.join(folder_path, "blastn_dict")
    os.makedirs(dict_folder_path, exist_ok=True)
    dict_file_path_out = os.path.join(dict_folder_path, os.path.basename(dict_path))
    blastn_dic(path_input=dict_path, path_output=dict_file_path_out)

    # Read input CSV and convert to JSON format for filtering
    data = pd.read_csv(csv_path, sep=",", header=0)
    
    # Create a JSON structure from the CSV data
    json_data = {}
    for index, row in data.iterrows():
        name_id = f"{row['sseqid']}_{row['sstrand']}_{row['sstart']}-{row['send']}"
        if name_id not in json_data:
            json_data[name_id] = []
        json_data[name_id].append([row['sseqid'], row['sstrand'], row['sstart'], row['send']])
    
    # Save the JSON data to a temporary file
    json_path = os.path.join(folder_path, "input_data.json")
    with open(json_path, 'w') as f:
        json.dump(json_data, f, indent=4)
    
    # Apply SIDER filter
    print_message_box("Applying SIDER filter")
    filtered_json_path = json_sider_filter(
        json_file=json_path,
        folder_path=folder_path,
        dict_path=dict_file_path_out,
        word_size=word_size,
        evalue=evalue
    )

    # Convert filtered JSON to CSV files (positive and negative)
    print_message_box("Converting filtered data to CSV")
    sider_json_to_csv(
        json_file=filtered_json_path,
        folder_path=folder_path,
        dict_path=dict_file_path_out,
        recaught_file=recaught_file_path,
        recaught_threshold=recaught_threshold,
        word_size=word_size,
        perc_identity=identity
    )
    
    # The sider_json_to_csv function already saves the output files to:
    # - os.path.join(os.path.dirname(folder_path), "siders_df.csv")
    # - os.path.join(os.path.dirname(folder_path), "non_siders_df.csv")
    
    positive_path = os.path.join(output_dir, "siders_df.csv")
    negative_path = os.path.join(output_dir, "non_siders_df.csv")

    message_to_print = (
        f"SIDER filtering completed. Results saved to:\n"
        f"- Positive: {positive_path}\n"
        f"- Negative: {negative_path}"
    )
    print_message_box(message_to_print)
    
    # Add right to groups and users
    try:
        subprocess.run(["chmod", "-R", "a+w", folder_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while changing permissions: {e}")

print_message_box("SIDER filtering completed")