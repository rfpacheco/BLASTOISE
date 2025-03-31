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
def general_blastn_blaster(query_path, dict_path, word_size, perc_identity=None, evalue=None):
    """
        Executes a BLASTN command-line query, processes its output, and returns
        the results as a Pandas DataFrame.

        This function builds and executes a BLASTN command using the provided query
        file and database path, while allowing customization of parameters like
        word size, percentage identity, and E-value. The results from the BLASTN
        query are parsed into a Pandas DataFrame for easier subsequent analysis.

        Parameters:
            query_path (str): Path to the query file for BLAST analysis.
            dict_path (str): Path to the database against which to run the BLAST
                analysis.
            word_size (int): Word size for the BLAST search, specifying the minimum
                length for exact matches.
            perc_identity (Optional[float]): Percentage identity threshold for
                matches (if provided, it is included in the BLAST query).
            evalue (Optional[float]): E-value threshold for matches (if provided,
                it is included in the BLAST query).

        Returns:
            pd.DataFrame: A Pandas DataFrame containing the BLAST output with the
            following columns:
                - qseqid: Query sequence ID.
                - sseqid: Subject sequence ID.
                - sstrand: Strand of the subject sequence.
                - qstart: Start position in the query sequence.
                - qend: End position in the query sequence.
                - sstart: Start position in the subject sequence.
                - send: End position in the subject sequence.
                - evalue: Expect value for the match.
                - bitscore: Score of the match in bits.
                - length: Alignment length denoting the length of the match.
                - qlen: Length of the query sequence.
                - slen: Length of the subject sequence.

        Raises:
            TypeError: If parameters do not match the expected type or are missing.
            ValueError: If the BLAST query execution or processing fails.

        Notes:
            This function uses the BLASTN command-line tool, and 'blastn' must be
            available and properly configured in the system path or environment. It
            runs the command through a subprocess with Bash as the shell.
    """
    cmd = f"blastn -word_size {word_size} -query {query_path} -db {dict_path}"

    if perc_identity is not None:
        cmd += f" -perc_identity {perc_identity}"

    if evalue is not None:
        cmd += f" -evalue {evalue}"

    # Use the more detailed output format from simple_blastn_blaster
    cmd += " -outfmt '10 qseqid sseqid sstrand qstart qend sstart send evalue bitscore length qlen slen'"

    data = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True,
                          executable='/usr/bin/bash')
    data = data.stdout

    # Parse the output into a DataFrame
    data_df = pd.DataFrame([x.split(",") for x in data.split("\n") if x])

    if not data_df.empty:
        data_df.columns = ["qseqid", "sseqid", "sstrand", "qstart", "qend", "sstart", "send",
                           "evalue", "bitscore", "length", "qlen", "slen"]

        # Convert numeric columns to numeric type
        numeric_cols = ["qstart", "qend", "sstart", "send", "evalue", "bitscore",
                        "length", "qlen", "slen"]
        data_df[numeric_cols] = data_df[numeric_cols].apply(pd.to_numeric)
    else:
        # Create an empty DataFrame with the correct columns
        data_df = pd.DataFrame(columns=["qseqid", "sseqid", "sstrand", "qstart", "qend", "sstart", "send",
                                        "evalue", "bitscore", "length", "qlen", "slen"]
        )

    return data_df