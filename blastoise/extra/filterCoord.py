import argparse
import os
import sys
import pandas as pd
import subprocess

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


from modules.blaster import blastn_dic
from extra.main_functions import coordinates_corrector, json_sider_filter, sider_json_to_csv
from modules.aesthetics import print_message_box


# ======================================================================
def parse_arguments():
    """
    Parses command-line arguments necessary for running a bioinformatics pipeline.

    This function sets up an argument parser with various options for input and
    output file paths as well as parameters controlling the behavior of the
    pipeline. It is specifically designed to handle sequence analysis tasks and
    accepts arguments related to BLASTN configuration and sequencing thresholds.

    Returns:
        argparse.Namespace: An object containing parsed command-line arguments.

    Raises:
        SystemExit: If the command-line arguments are invalid or miss required
        options, the parser will terminate the program and print the help
        message.

    Arguments:
        -f, --file (str): Path to the input CSV file (required).
        -d, --dict_path (str): Path to the genome fasta file (required).
        -rf, --recaught_file (str): Path to the recaught file (required).
        -rt, --recaught_threshold (float): Recaught threshold (default: 1.0E-03).
        -ws, --word_size (int): word_size value for BLASTN (required).
        -min, --min_length (int): Minimum length of the sequence to be considered
        (required).
        -e, --evalue (float): E-value threshold for BLASTN (default: 1.0E-09).
        -i, --identity (int): Minimum identity of the sequence to be considered
        (required).
    """
    parser = argparse.ArgumentParser(description="")
    # noinspection DuplicatedCode
    parser.add_argument("-f", "--file", type=str, required=True, help="Path to the input CSV file.")
    parser.add_argument("-d", "--dict_path", type=str, required=True, help="Path to the genome fasta file.")
    parser.add_argument("-rf", "--recaught_file", type=str, required=True, help="Path to the recaught file.")
    parser.add_argument("-rt", "--recaught_threshold", type=float, default=1.0E-03, help="")
    parser.add_argument("-ws", "--word_size", type=int, required=True, help="word_size value of BLASTN")
    parser.add_argument("-min", "--min_length", type=int, required=True, help="Minimum length of the sequence to be considered.")
    parser.add_argument("-e", "--evalue", type=float, default=1.0E-09, help="E-value threshold for BLASTN.")
    parser.add_argument("-i", "--identity", type=int, required=True, help="Minimum identity of the sequence to be considered.")
    return parser.parse_args()

# ======================================================================
if __name__ == "__main__":
    # Load all arguments
    args = parse_arguments()
    csv_path = os.path.expanduser(args.file)
    dict_path = os.path.expanduser(args.dict_path)
    recaught_file_path = os.path.expanduser(args.recaught_file)
    recaught_threshold = args.recaught_threshold
    word_size = args.word_size
    min_length = args.min_length
    evalue = args.evalue
    identity = args.identity

    # Prepare subfolder
    csv_parent_path = os.path.dirname(csv_path)
    folder_path = os.path.join(csv_parent_path, "tmpFilterCoord")  # tmp folder to place all data before deleting

    # Prepare BLASTn dict
    dict_folder_path = os.path.join(folder_path, "blastn_dict")
    os.makedirs(dict_folder_path, exist_ok=True)
    dict_file_path_out = os.path.join(dict_folder_path, os.path.basename(dict_path))
    blastn_dic(path_input=dict_path, path_output=dict_file_path_out)

    # # First SIDER filter
    data = pd.read_csv(csv_path, sep=",", header=0)
    print_message_box("Correcting coordinates 1/3")
    json_path = coordinates_corrector(
        df=data,
        dict_path=dict_file_path_out,
        folder_path=folder_path,
        word_size=word_size,
        min_length=min_length
    )

    ## sider filter to new coordinates
    print_message_box("Applying filter 2/3")
    filtered_json_path = json_sider_filter(
        json_file=json_path,
        folder_path=folder_path,
        dict_path=dict_file_path_out,
        word_size=word_size,
        evalue=evalue
    )

    ## from json to csv
    print_message_box("Re-caught data and output files 3/3")
    sider_json_to_csv(
        json_file=filtered_json_path,
        folder_path=folder_path,
        dict_path=dict_file_path_out,
        recaught_file=recaught_file_path,
        recaught_threshold = recaught_threshold,
        word_size=word_size,
        perc_identity=identity
    )

    # Add right to groups and users
    try:
        subprocess.run(["chmod", "-R", "a+w", folder_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while changing permissions: {e}")

print_message_box("THE END")