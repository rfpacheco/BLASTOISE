import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from blastoise.modules.blaster import blastn_dic

# ======================================================================================================================
def simple_fasta_creator(sequence, fasta_index, fasta_output_path):
    rec = SeqRecord(Seq(sequence),
                    id="Seq_" + str(fasta_index),
                    description=""
                    )
    SeqIO.write(rec, fasta_output_path, "fasta")

# ======================================================================================================================
def csv_to_fasta_creator(csv_data, fasta_output_path):
    matrix = []
    for csv_index, sequence in csv_data.iterrows():
        rec = SeqRecord(Seq(sequence['sseq']),
                        id=f"Seq_{csv_index}_{sequence['sseqid']}",
                        description=""
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")

# ======================================================================================================================
# noinspection DuplicatedCode
def blastn_blaster(query_path, dict_path, evalue, word_size):
    cmd = "blastn -word_size " + str(word_size) + " -query " \
          + query_path + " -db " \
          + dict_path \
          + " -evalue " + str(evalue) \
          + " -outfmt 10"
    blast_data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return blast_data

# ======================================================================================================================
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

# ======================================================================================================================
# TODO: make it simpler or more modular
def sider_filter(df, dict_path, folder_path, word_size, recaught_file):
    matches = pd.Series([False] * df.shape[0])
    not_matches = pd.Series([True] * df.shape[0])
    accepted = 0
    rejected = 0

    for index, row in df.iterrows():
        print("="*50)
        print(f"Analyzing row {index + 1} of {df.shape[0]}")

        # TODO: change to tmp .fasta files
        fasta_path = os.path.join(folder_path, "mySequence.fasta")
        print(row)
        simple_fasta_creator(row["sseq"], index, fasta_path)
        blastn_data = blastn_blaster(
            query_path=fasta_path,
            dict_path=dict_path,
            evalue=1.0E-09,  # TODO: implement as argument
            word_size=word_size
        )

        # noinspection DuplicatedCode
        if blastn_data.count("\n") <= 1:  # Only match with itself
            not_matches[index] = True
            rejected += 1
            print("\t\tREJECTED")
        else:
            blastn_data = blastn_data.strip().split("\n")
            blast_data_df = pd.DataFrame([x.split(",") for x in blastn_data if x])
            if blast_data_df[1].nunique() >= 5:  # 5 is from out SIDER test  # TODO: improve as parameter
                matches[index] = True
                accepted += 1
                print("\t\tACCEPTED")
            else:
                not_matches[index] = True
                rejected += 1
                print("\t\tREJECTED")
        print(f"\t\t\t\t\tAccepted: {accepted} - Rejected: {rejected}")

    # noinspection DuplicatedCode
    print(f"The total number of matches is: {matches.sum()} out of {df.shape[0]}")
    print(f"The percentage of matches is: {round(matches.sum() / df.shape[0] * 100, 2)}%")
    print("~"*50)
    print(f"The total number of not matches is: {not_matches.sum()} out of {df.shape[0]}")
    print(f"The percentage of not matches is: {round(not_matches.sum() / df.shape[0] * 100, 2)}%")

    yes_data = df[matches]
    no_data = df[~matches]

    # Recaught elements in 'no_data'
    no_data_recaught_folder_path = os.path.join(folder_path, 'recaught_in_negatives')
    os.makedirs(no_data_recaught_folder_path, exist_ok=True)

    # Create fasta about the negative files
    no_data_fasta_path = os.path.join(no_data_recaught_folder_path, 'no_data.fasta')
    csv_to_fasta_creator(no_data, no_data_fasta_path)

    # Make a BLASTn dict with that:
    blastn_dic(no_data_fasta_path, no_data_fasta_path)

    # Search for recaught data
    # TODO: perc identity could be an argument
    caught_data = recaught_blast(recaught_file, no_data_fasta_path, 60, word_size)
    # noinspection DuplicatedCode
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
    return final_yes_data, final_no_data





