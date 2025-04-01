import argparse
import pandas as pd
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from modules.blaster import blastn_dic
from extra.second_functions import get_sequence, simple_fasta_creator, csv_to_fasta_creator, general_blastn_blaster

# ======================================================================
# Defining parse arguments
# ======================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter elements based on chromosome count and e-value.")
    parser.add_argument("-f", "--file", type=str, required=True, help="Path to the input CSV file.")
    parser.add_argument("-d", "--dict_path", type=str, required=True, help="Path to the genome fasta file.")
    parser.add_argument("-ws", "--word_size", type=str, required=True, help="word_size value of BLASTN")
    parser.add_argument("-rf", "--recaught_file", type=str, required=True, help="Path to the recaught file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the folder to save the CSV files.")
    return parser.parse_args()

# ======================================================================
# Defining needed functions
# ======================================================================
# Define function to extract fasta sequences
def fasta_extractor(pathfile, outfile, extract_list):
    with open(outfile, "w") as out_file:
        # Remember "enumerate" starts in "1"
        for count, fasta in enumerate(SeqIO.parse(open(pathfile), "fasta"), start=0):  # from Bio import SeqIO
            # name, sequence = fasta.id, str(fasta.seq)
            if count in extract_list:
                SeqIO.write(fasta, out_file, "fasta")

# ======================================================================
# Main function
# ======================================================================
# noinspection DuplicatedCode
def main():
    # Prepare files, paths & tmp folder.
    args = parse_arguments()
    file_path = os.path.expanduser(args.file)
    dict_path = os.path.expanduser(args.dict_path)
    recaught_file_path = os.path.expanduser(args.recaught_file)
    output_path = os.path.expanduser(args.output)
    word_size = args.word_size

    folder_path = os.path.join(output_path, "tmpTrueSiderFilter")
    os.makedirs(folder_path, exist_ok=True)

    data = pd.read_csv(file_path, sep=",", header=0)  # read the file
    data["sseq"] = None # New None column to fill later

    # Prepare pandas Series False arrays with the row length of the data file
    matches = pd.Series([False] * data.shape[0])  # For matches
    not_matches = pd.Series([False] * data.shape[0])  # For not matches
    accepted = 0
    rejected = 0

    # Create the blast dic
    genome_folder_path = os.path.join(folder_path, "blastn_dic")
    os.makedirs(genome_folder_path, exist_ok=True)
    genome_file_path_out = os.path.join(genome_folder_path, os.path.basename(dict_path))
    blastn_dic(dict_path, genome_file_path_out)

    for index, row in data.iterrows():
        print("=" * 50)
        index = hash(index)
        print(f"Analyzing row {index + 1} of {data.shape[0]}")

        # Get seq
        sequence = get_sequence(
            start_coor=row["sstart"],
            end_coor=row["send"],
            strand="plus",  # TODO: make an arg
            chromosome=row["sseqid"],
            path_genome=genome_file_path_out
        )
        # Add sequence to the col "sseq" pandas series & to the original data
        row["sseq"] = sequence
        data.loc[index, "sseq"] = sequence

        fasta_path = os.path.join(folder_path, "mySequence.fasta")
        print(row)
        simple_fasta_creator(sequence, index, fasta_path)
        blastn_data = general_blastn_blaster(
            query_path=fasta_path,
            dict_path=genome_file_path_out,
            evalue=1.0E-09,
            word_size=word_size
        )

        # TODO: improve reusability
        if blastn_data.empty or blastn_data.shape[0] <= 1:
            not_matches[index] = True
            rejected += 1
            print("\t\tREJECTED")
        else:
            if blastn_data["sseqid"].nunique() >= 5:
                matches[index] = True
                accepted += 1
                print("\t\tACCEPTED")
            else:
                not_matches[index] = True
                rejected += 1
                print("\t\tREJECTED")
        print(f"\t\t\t\t\tAccepted: {accepted} - Rejected: {rejected}")

    # TODO: improve reusability
    print(f"The total number of matches is: {matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of matches is: {round(matches.sum() / data.shape[0] * 100, 2)}%")
    print("~" * 50)
    print(f"The total number of not matches is: {not_matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of not matches is: {round(not_matches.sum() / data.shape[0] * 100, 2)}%")

    yes_data = data[matches]
    no_data = data[~matches]

    # if no_data has lines
    if not no_data.empty:
        # Recaught elements in 'no_data'
        no_data_recaught_folder_path = os.path.join(folder_path, 'recaught_in_negatives')
        os.makedirs(no_data_recaught_folder_path, exist_ok=True)

        # Create fasta about the negative files
        no_data_fasta_path = os.path.join(no_data_recaught_folder_path, 'no_data.fasta')
        csv_to_fasta_creator(no_data, no_data_fasta_path)

        # Make a BLASTn dict with that:
        blastn_dic(no_data_fasta_path, no_data_fasta_path)

        # Search for recaught data
        caught_data = general_blastn_blaster(
            query_path=recaught_file_path,
            dict_path=no_data_fasta_path,
            word_size=60
        )
        # TODO: improve reusability
        if not caught_data.empty:
            # Remove ones with an evalue <= 10**-3
            caught_data = caught_data[caught_data['evalue'] <= 1.0E-03].sort_values(by=['evalue'])
            print("")
            print("*" * 50)
            print(f"\nRecaught data: {caught_data.shape[0]} elements")

            # Create a column with the number in "sseqid"
            caught_data['index'] = caught_data['sseqid'].str.extract(r'_(\d+)_')
            caught_data['index'] = pd.to_numeric(caught_data['index'])

            # Get a list with the index column
            index_list = caught_data['index'].sort_values().unique().tolist()

            # Extract sequences from the 'no_data'
            no_data_recaught = no_data[no_data.index.isin(index_list)]

            # Join yes_data and no_data_recaught
            final_yes_data = pd.concat([yes_data, no_data_recaught], axis=0, ignore_index=True)
            final_yes_data.sort_values(by=['sseqid', 'sstart'], inplace=True)

            # Remove no_data_recaught from no data
            final_no_data = pd.concat([no_data, no_data_recaught]).drop_duplicates(keep=False)

            # Print results:
            print(f"\n\t - Accepted data + recaught: {final_yes_data.shape[0]} elements")
            print(f"\t - Rejected data - recaught: {final_no_data.shape[0]} elements")

        else:
            final_yes_data = yes_data
            final_no_data = no_data
            print("\n\t - No recaught data")

    else:
        final_yes_data = yes_data
        final_no_data = no_data

    # Save both data:
    final_yes_data_path = os.path.join(output_path, 'passed_sider_test.csv')
    final_no_data_path = os.path.join(output_path, 'not_passed_sider_test.csv')
    final_yes_data.to_csv(final_yes_data_path, index=False, header=True)
    final_no_data.to_csv(final_no_data_path, index=False, header=True)

    # Remove the folder folder_path
    if os.path.exists(folder_path):
        os.system(f'rm -rf {folder_path}')


if __name__ == "__main__":
    main()