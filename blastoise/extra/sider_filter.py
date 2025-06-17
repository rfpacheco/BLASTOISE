import argparse
import os
import sys
import pandas as pd
import subprocess
from joblib import Parallel, delayed

# Add the parent directory of 'blastoise' to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from modules.blaster import blastn_dic
from extra.second_functions import general_blastn_blaster
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

    # Read input CSV
    data = pd.read_csv(csv_path, sep=",", header=0)

    # Group data by sequence identifiers
    print_message_box("Preparing data for analysis")
    data['name_id'] = data.apply(lambda row: f"{row['sseqid']}_{row['sstrand']}_{row['sstart']}-{row['send']}", axis=1)

    # Apply SIDER filter directly on the DataFrame
    print_message_box("Applying SIDER filter")

    # Function to process a single sequence
    def process_sequence(row):
        # Use the sequence directly from the CSV
        if 'sseq' not in row:
            return {'name_id': row['name_id'], 'status': 'Rejected'}

        sequence = row['sseq']

        # Create a temporary query file for BLASTN
        name_id = row['name_id']
        query = f"<(echo -e '>{name_id}\\n{sequence}')"

        # Run BLASTN
        blastn_df = general_blastn_blaster(
            query_path=query,
            dict_path=dict_file_path_out,
            word_size=word_size,
            evalue=evalue
        )

        # Check BLASTN results
        if not blastn_df.empty and blastn_df["sseqid"].nunique() >= 5:
            return {'name_id': name_id, 'status': 'Accepted'}
        else:
            return {'name_id': name_id, 'status': 'Rejected'}

    # Process sequences in parallel
    print(f"Analyzing {len(data)} sequences...")
    results = Parallel(n_jobs=-1)(
        delayed(process_sequence)(row) for _, row in data.iterrows()
    )

    # Convert results to DataFrame
    results_df = pd.DataFrame(results)

    # Merge results with original data
    data = pd.merge(data, results_df, on='name_id', how='left')

    # Split into accepted and rejected
    accepted_data = data[data['status'] == 'Accepted'].copy()
    rejected_data = data[data['status'] == 'Rejected'].copy()

    print(f"Accepted sequences: {len(accepted_data)}")
    print(f"Rejected sequences: {len(rejected_data)}")

    # Process recaught data
    print_message_box("Processing recaught data")

    # Create a FASTA file from the rejected data for BLASTN
    if not rejected_data.empty:
        # The sequence data should already be in the CSV
        # No need to call get_data_sequence

        # Create a FASTA file for BLASTN
        fasta_file_path = os.path.join(folder_path, "negative_database.fasta")
        with open(fasta_file_path, 'w') as f:
            for idx, row in rejected_data.iterrows():
                f.write(f">Seq_{idx}_{row['sseqid']}\n{row['sseq']}\n")

        # Create BLASTN database
        blastn_dic(path_input=fasta_file_path, path_output=dict_file_path_out)

        # Run BLASTN on recaught file
        caught_data = general_blastn_blaster(
            query_path=recaught_file_path,
            dict_path=dict_file_path_out,
            perc_identity=identity,
            word_size=word_size
        )

        # Process recaught data
        if not caught_data.empty:
            # Filter by e-value
            caught_data = caught_data[caught_data["evalue"] <= recaught_threshold].sort_values(by=["evalue"])
            print(f"\nRecaught data: {caught_data.shape[0]} elements")

            # Extract sequence indices from sseqid
            caught_data["index"] = caught_data["sseqid"].str.extract(r"Seq_(\d+)_").astype(int)

            # Get indices of recaught sequences
            recaught_indices = caught_data["index"].unique()

            # Move recaught sequences from rejected to accepted
            recaught_rows = rejected_data.iloc[recaught_indices].copy()
            if not recaught_rows.empty:
                recaught_rows['status'] = 'Accepted (Recaught)'
                accepted_data = pd.concat([accepted_data, recaught_rows], ignore_index=True)
                rejected_data = rejected_data.drop(recaught_indices)

                print(f"Recaught sequences: {len(recaught_rows)}")
                print(f"Updated accepted sequences: {len(accepted_data)}")
                print(f"Updated rejected sequences: {len(rejected_data)}")

    # Save the final datasets
    positive_path = os.path.join(output_dir, "siders_df.csv")
    negative_path = os.path.join(output_dir, "non_siders_df.csv")

    accepted_data.to_csv(positive_path, index=False)
    rejected_data.to_csv(negative_path, index=False)

    message_to_print = (
        f"SIDER filtering completed. Results saved to:\n"
        f"- Positive: {positive_path}\n"
        f"- Negative: {negative_path}"
    )
    print_message_box(message_to_print)

    # Add right to groups and users
    try:
        subprocess.run(["chmod", "-R", "a+w", folder_path], check=True)
        subprocess.run(["chmod", "a+w", positive_path], check=True)
        subprocess.run(["chmod", "a+w", negative_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while changing permissions: {e}")

print_message_box("SIDER filtering completed")
