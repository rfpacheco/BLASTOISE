from typing import Optional, NoReturn
import pandas as pd
import os
import subprocess
import argparse


# Now let's get the sequences
def get_data_sequence(data: pd.DataFrame, strand: str, genome_fasta: str) -> pd.DataFrame:
    """
    :param data: A pandas DataFrame containing BLAST search results with columns 'sseqid', 'sstart', and 'send'.
    :param strand: The DNA strand ('plus' or 'minus') indicating the sequence orientation.
    :param genome_fasta: The path to the genome fasta file used for reference.
    :return: A pandas DataFrame with BLAST coordinates and the corresponding sequences.
    """
    sequences = []
    for index, row in data.iterrows():
        print(f"\n\t - Getting sequence {index + 1}/{data.shape[0]}")
        chrom = row["sseqid"]
        start = row["sstart"]
        end = row["send"]
        blast_cmd = f"blastdbcmd -db {genome_fasta} -entry {chrom} -range {start}-{end} -strand {strand} -outfmt %s"

        sequence = subprocess.check_output(blast_cmd, shell=True, universal_newlines=True).replace('\n', '')

        sequences.append({
            "sseqid": chrom,
            "sstart": start,
            "send": end,
            "sstrand": strand,
            "sseq": sequence
        })

    sequences_df = pd.DataFrame(sequences)
    return sequences_df


# Prepare dict creation
def blastn_dic(path_input: str, path_output: str) -> Optional[None]:
    """
    :param path_input: Path to the input file containing nucleotide sequences to be used in creating the BLAST database.
    :param path_output: Path where the resulting BLAST database should be stored.
    :return: None if the operation completes successfully.
    """
    # "parse_seqids" is used to keep the sequence ID in the output.
    blast_cmd = f"makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}"
    subprocess.run(blast_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main(df_path: str, genome_fasta_path: str, example_strand: str) -> NoReturn:
    """
    :param df_path: Path to the CSV file containing genomic 0.1.data.
    :param genome_fasta_path: Path to the genome FASTA file used for sequence extraction.
    :param example_strand: Strand information for sequence extraction.
    :return: The function does not return any value. It processes the genomic 0.1.data, merges strands using BEDOPS, extracts sequences, and writes the merged sequences to a CSV file.
    """
    df_path = os.path.expanduser(df_path)
    df = pd.read_csv(df_path, sep=',', header=0)
    df = df[['sseqid', 'sstart', 'send']].copy()  # take only needed files

    df_parent_path = os.path.dirname(df_path)
    bedops_folder_path = os.path.join(df_parent_path, "join_strands")
    os.makedirs(bedops_folder_path, exist_ok=True)

    bedops_file_path = os.path.join(bedops_folder_path, 'to_join_strands.bed')
    df.sort_values(by=['sseqid', 'sstart'], inplace=True)  # Important for BEDOPS
    df.to_csv(bedops_file_path, sep='\t', index=False, header=False)

    # Bedops merge
    cmd = f"bedops --merge {bedops_file_path}"
    df_merged = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    df_merged = pd.DataFrame([x.split("\t") for x in df_merged.split("\n") if x],
                             columns=["sseqid", "sstart", "send"])
    df_merged[["sstart", "send"]] = df_merged[["sstart", "send"]].apply(pd.to_numeric)
    print(f"From a total of {df.shape[0]} there are {df_merged.shape[0]} sequences after BEDOPS merge")

    genome_fasta_path = os.path.expanduser(genome_fasta_path)
    genome_fasta_name = os.path.basename(genome_fasta_path)
    blast_dic_folder_path = os.path.join(bedops_folder_path, "blast_dic")
    blast_dic_path = os.path.join(blast_dic_folder_path, genome_fasta_name)
    os.makedirs(blast_dic_folder_path, exist_ok=True)
    blastn_dic(genome_fasta_path, blast_dic_path)

    print("\nGetting data sequences:")
    print("="*20)
    df_merged_seqs = get_data_sequence(df_merged, example_strand, blast_dic_path)

    # Get the parent path of 'df_path'
    parent_path = os.path.dirname(df_path)
    merged_file_name = os.path.basename(df_path)
    merged_file_name = merged_file_name.replace(".csv", "")
    merged_file_name = f"{merged_file_name}_merged.csv"
    merged_path = os.path.join(parent_path, merged_file_name)
    df_merged_seqs.to_csv(merged_path, sep=',', header=True, index=False)
    print(f"Saved data in: {merged_path}")

    if os.path.exists(bedops_folder_path):
        os.system(f"rm -rf {bedops_folder_path}")

    # Add right to groups and users
    subprocess.run(["chmod", "-R", "a+w", merged_path], check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some genomic 0.1.data.")
    parser.add_argument("-f", "--file", type=str, help="Path to the CSV file containing the sequences data.")
    parser.add_argument("-d", "--dict_path", type=str, help="Path to the genome FASTA file.")
    parser.add_argument("--strand", type=str, help="Strand to use (plus or minus).")

    args = parser.parse_args()
    main(args.file, args.dict_path, args.strand)

