import argparse  # type: ignore

import pandas as pd
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter elements based on chromosome count and e-value.")
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('-d', '--dict_path', type=str, required=True, help='Path to the genome fasta file.')
    parser.add_argument('--word_size', type=str, required=True, help='word_size value of BLASTN')
    parser.add_argument('--recaught_file', type=str, required=True, help='Path to the recaught file.')
    return parser.parse_args()

def blastn_dic(path_input, path_output):
    # "parse_seqids" is used to keep the sequence ID in the output.
    cmd = f"makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}"
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def blastn_blaster(query_path, dict_path, evalue, word_size):
    cmd = "blastn -word_size " + str(word_size) + " -query " \
          + query_path + " -db " \
          + dict_path \
          + " -evalue " + str(evalue) \
          + " -outfmt 10"
    blast_data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return blast_data

def recaught_blast(query_path, dict_path, perc_identity, word_size):
    cmd = "blastn -word_size " + str(word_size) + " -query " \
        + query_path + " -db " \
        + dict_path \
        + " -perc_identity " + str(perc_identity) \
        + " -outfmt '10 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen'"
    recaught_df = subprocess.check_output(cmd, shell=True, universal_newlines=True)  # Important the E value
    recaught_df = pd.DataFrame([x.split(",") for x in recaught_df.split("\n") if x])
    if not recaught_df.empty:
        recaught_df.columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
        recaught_df[['pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']] = recaught_df[['pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']].apply(pd.to_numeric)
    else:
        recaught_df = pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
    return recaught_df

def fasta_creator(sequence, fasta_index, fasta_output_path):
    rec = SeqRecord(Seq(sequence),
                    id="Seq_" + str(fasta_index),
                    description="Leishmania"
                    )
    SeqIO.write(rec, fasta_output_path, "fasta")

def csv_to_fasta_creator(csv_data, fasta_output_path):
    matrix = []
    for csv_index, sequence in csv_data.iterrows():
        rec = SeqRecord(Seq(sequence['sseq']),
                        id=f"Seq_{csv_index}_{sequence['sseqid']}",
                        description="Leishmania"
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")

# Define function to extract fasta sequences
def fasta_extractor(pathfile, outfile, extract_list):
    with open(outfile, "w") as out_file:
        # Remember "enumerate" starts in "1"
        for count, fasta in enumerate(SeqIO.parse(open(pathfile), "fasta"), start=0):  # from Bio import SeqIO
            # name, sequence = fasta.id, str(fasta.seq)
            if count in extract_list:
                SeqIO.write(fasta, out_file, "fasta")

if __name__ == "__main__":
    args = parse_arguments()
    file_1_path = os.path.expanduser(args.file)
    file_1_parent = os.path.dirname(file_1_path)
    folder_path = os.path.join(file_1_parent, 'true_sider_filter')
    os.makedirs(folder_path, exist_ok=True)

    data = pd.read_csv(args.file, sep=",", header=0)
    matches = pd.Series([False] * data.shape[0])
    not_matches = pd.Series([False] * data.shape[0])
    accepted = 0
    rejected = 0

    # Create the blast dic
    genome_file_path = os.path.expanduser(args.dict_path)
    genome_folder_path = os.path.join(folder_path, 'blastn_dic')
    os.makedirs(genome_folder_path, exist_ok=True)
    genome_file_path_out = os.path.join(genome_folder_path, os.path.basename(genome_file_path))
    blastn_dic(genome_file_path, genome_file_path_out)

    for index, row in data.iterrows():
        print("="*50)
        print(f"Analyzing row {hash(index) + 1} of {data.shape[0]}")

        fasta_path = os.path.join(folder_path, "mySequence.fasta")
        print(row)
        fasta_creator(row["sseq"], index, fasta_path)
        blastn_data = blastn_blaster(fasta_path, genome_file_path_out, 1.0E-09, args.word_size)

        if blastn_data.count("\n") <= 1:
            not_matches[index] = True
            rejected += 1
            print("\t\tREJECTED")
        else:
            blastn_data = blastn_data.strip().split("\n")
            blast_data_df = pd.DataFrame([x.split(",") for x in blastn_data if x])
            if blast_data_df[1].nunique() >= 5:
                matches[index] = True
                accepted += 1
                print("\t\tACCEPTED")
            else:
                not_matches[index] = True
                rejected += 1
                print("\t\tREJECTED")
        print(f"\t\t\t\t\tAccepted: {accepted} - Rejected: {rejected}")

    print(f"The total number of matches is: {matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of matches is: {round(matches.sum() / data.shape[0] * 100, 2)}%")
    print("~"*50)
    print(f"The total number of not matches is: {not_matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of not matches is: {round(not_matches.sum() / data.shape[0] * 100, 2)}%")

    yes_data = data[matches]
    # yes_data.to_csv(os.path.join(args.output_dir, "positives_testing_elements.csv"), index=False, header=True)
    no_data = data[~matches]
    # no_data.to_csv(os.path.join(args.output_dir, "negatives_testing_elements.csv"), index=False, header=True)


    # Recaught elements in 'no_data'
    no_data_recaught_folder_path = os.path.join(folder_path, 'recaught_in_negatives')
    os.makedirs(no_data_recaught_folder_path, exist_ok=True)

    # Create fasta about the negative files
    no_data_fasta_path = os.path.join(no_data_recaught_folder_path, 'no_data.fasta')
    csv_to_fasta_creator(no_data, no_data_fasta_path)

    # Make a BLASTn dict with that:
    blastn_dic(no_data_fasta_path, no_data_fasta_path)

    # Search for recaught data
    caught_data = recaught_blast(args.recaught_file, no_data_fasta_path, 60, args.word_size)
    if not caught_data.empty:
        # Remove ones with an evalue <= 10**-3
        caught_data = caught_data[caught_data['evalue'] <= 1.0**-3].sort_values(by=['evalue'])
        print("")
        print("*"*50)
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

    # Save both data:
    final_yes_data_path = os.path.join(file_1_parent, 'final_yes_data.csv')
    final_no_data_path = os.path.join(file_1_parent, 'final_no_data.csv')
    final_yes_data.to_csv(final_yes_data_path, index=False, header=True)
    final_no_data.to_csv(final_no_data_path, index=False, header=True)

    # Remove the folder folder_path
    if os.path.exists(folder_path):
        os.system(f'rm -rf {folder_path}')


