from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd
import subprocess

# ======================================================================
def simple_fasta_creator(sequence, fasta_index, fasta_output_path):
    rec = SeqRecord(Seq(sequence),
                    id="Seq_" + str(fasta_index),
                    description=""
                    )
    SeqIO.write(rec, fasta_output_path, "fasta")

# ======================================================================
def csv_to_fasta_creator(csv_data, fasta_output_path):
    matrix = []
    for csv_index, sequence in csv_data.iterrows():
        rec = SeqRecord(Seq(sequence['sseq']),
                        id=f"Seq_{csv_index}_{sequence['sseqid']}",
                        description=""
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")

# ======================================================================
# noinspection DuplicatedCode
def get_sequence(start_coor, end_coor, strand, chromosome, path_genome):
    cmd = f'blastdbcmd -db {path_genome} -entry {chromosome} -range {start_coor}-{end_coor} -strand {strand} -outfmt %s'
    sequence = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')
    sequence = sequence.stdout.strip()
    return sequence

# ======================================================================
# noinspection DuplicatedCode
def get_sequence_json_to_csv(start_coor, end_coor, strand, chromosome, path_genome):
    cmd = (
        f'blastdbcmd -db {path_genome} '
        f'-entry {chromosome} '
        f'-range {start_coor}-{end_coor} '
        f'-strand {strand} '
        f'-outfmt %s'
    )
    sequence = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')
    sequence = sequence.stdout.strip()
    return sequence

# ======================================================================================================================
# noinspection DuplicatedCode
def bedops_merge(input_df, path_folder):
    # Create a temporary bed file
    path_bedops_file = os.path.join(path_folder, "tmp.bed")
    data_bedops = input_df[['qstart', 'qend']].copy()  # in qstart and qend I don't have the "minus" coordinates problem
    data_bedops.insert(0, 'new_column', 'test')  # Add a new column with every row with the same value 'test'
    data_bedops.to_csv(path_bedops_file, sep="\t", header=False, index=False)

    # Call and process the bedops merge command
    cmd = f"bedops --merge {path_bedops_file}"
    data = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True,
                          executable='/usr/bin/bash')
    data = data.stdout  # Get the output
    data = pd.DataFrame([x.split("\t") for x in data.split("\n") if x], columns=['sseqid', 'qstart', 'qend'])
    data[['qstart', 'qend']] = data[['qstart', 'qend']].apply(pd.to_numeric)  # Convert to numeric

    return data

# ======================================================================
# BLAST FUNCTIONS
# ======================================================================
# noinspection DuplicatedCode
def blastn_blaster(query_path, dict_path, evalue, word_size):
    cmd = "blastn -word_size " + str(word_size) + " -query " \
          + query_path + " -db " \
          + dict_path \
          + " -evalue " + str(evalue) \
          + " -outfmt 10"
    blast_data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return blast_data

# noinspection DuplicatedCode
def simple_blastn_blaster(query_path, dict_path):  # TODO: change 11 to argument
    cmd = "blastn -word_size 11" \
        + " -query " + query_path \
        + " -db " + dict_path \
        + " -outfmt '10 qseqid sseqid sstrand pident qstart qend sstart send evalue bitscore length qlen qcovs slen mismatch gapopen gaps'"
    data = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')  # Important the E value
    data = data.stdout
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    if not data.empty:  # If the dataframe is not empty
        data.columns = ["qseqid", "sseqid", "sstrand", "pident", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "qlen", "qcovs", "slen", "mismatch", "gapopen", "gaps"]
        data[['pident',  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'length', 'qlen', 'qcovs', 'slen', 'mismatch', 'gapopen', 'gaps']] = data[['pident',  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'length', 'qlen', 'qcovs', 'slen', 'mismatch', 'gapopen', 'gaps']].apply(pd.to_numeric)
    else:  # If the dataframe is empty
        data = pd.DataFrame(columns=["qseqid", "sseqid", "sstrand", "pident", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "qlen", "qcovs", "slen", "mismatch", "gapopen", "gaps"])  # Create an empty dataframe
    return data


# noinspection DuplicatedCode
def json_blastn_blaster(query, path_genome, evalue):
    cmd = (
        f'blastn -word_size 15 '  # TODO: could be a parameter
        f'-query {query} '
        f'-db {path_genome} '
        f'-evalue {evalue} '
        f'-outfmt 10'
    )
    data = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')
    data = data.stdout
    data_df = pd.DataFrame(
        [x.split(',') for x in data.split('\n') if x],
        columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    )
    if not data_df.empty:
        data_df[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']] = data_df[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']].apply(pd.to_numeric)  # Convert to numeric
    else:  # If empty, return an empty dataframe
        return pd.DataFrame()

    return data_df


# noinspection DuplicatedCode
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
