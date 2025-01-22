# Needed modules
import argparse
import os
import subprocess
import pandas as pd
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#########
# NOTES #
#########
# This code is for my SIDER elements which are PLUS strand ALLWAYS


# =======================================================================================
# =======================================================================================
# add the main data frame to an argparse value
def parse_arguments():
    parser = argparse.ArgumentParser(description='Find more adjusted coordinates.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('-d', '--dict_path', type=str, required=True, help='Path to the genome fasta file.')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save the output CSV files.')
    return parser.parse_args()


# Define all functions needed
def blastn_dic(path_input, path_output):
    # "parse_seqids" is used to keep the sequence ID in the output.
    cmd = f"makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}"
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    

def blastn_blaster(query_path, dict_path):
    cmd = "blastn -word_size 15" \
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


def bedops_merge(input_df, path_folder):
    # Create a temporary bed file
    path_bedops_file = os.path.join(path_folder, "tmp.bed")
    data_bedops = input_df[['qstart', 'qend']].copy()  # in qstart and qend I don't have the "minus" coordinates problem
    data_bedops.insert(0, 'new_column', 'test')  # Add a new column with every row with the same value 'test'
    data_bedops.to_csv(path_bedops_file, sep="\t", header=False, index=False)

    # Call and process the bedops merge command
    cmd = f"bedops --merge {path_bedops_file}"
    data = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')
    data = data.stdout   # Get the output
    data = pd.DataFrame([x.split("\t") for x in data.split("\n") if x], columns=['sseqid', 'qstart', 'qend'])
    data[['qstart', 'qend']] = data[['qstart', 'qend']].apply(pd.to_numeric)  # Convert to numeric
    
    return data


if __name__ == '__main__':
    # 1) Load data
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)  # Create the output directory if it does not exist
    data = pd.read_csv(args.file, sep=',', header=0)
    print(f'Loading data from {args.file}')

    # 2) Create the BLASTn dictionary
    path_dict_out = os.path.join(args.output_dir, "blastn_dict")  # Create folder to store the dictionary
    blastn_dic(path_input=args.dict_path, path_output=path_dict_out)
    print(f'Creating BLASTn dictionary from {args.dict_path}\n')

    # 2) Create dict where the data will be stored
    main_dict = {}

    # 3) Create main loop
    for index, row in data.iterrows():
        # 3.1) Prepare the query data 
        name_id = f"{row['sseqid']}_{row['sstrand']}_{row['sstart']}-{row['send']}"  # Create the name_id with the original info of the seq.
        seq = row['sseq']
        query = f"<(echo -e '>{name_id}\n{seq}')"  # Create the query with the name_id and the seq in a bash tmp file
        start_coor = row['sstart']  # Get the start coordinate for later
        end_coor = row['send']  # Get the end coordinate for later
        strand_seq = row['sstrand']  # Get the strand of the sequence for later
        name_chr = row['sseqid']  # Get the name of the chromosome for later
        print(f'Analyzing row {index + 1}/{data.shape[0]} with name_id {name_id}')

        # 3.2) Call the blastn function 
        blastn_df = blastn_blaster(query_path=query, dict_path=path_dict_out)  # Call the blastn function

        # 3.3) Filter the blastn_df
        # remove the row with the same sstart and send values thand start_coor and end_coor in plus way
        blastn_df = blastn_df[~(
            ((blastn_df["sstart"] >= start_coor) & (blastn_df["sstart"] <= end_coor)) |  # (sstart is within the start and end coordinates OR 
            ((blastn_df["send"] <= end_coor) & (blastn_df["send"] >= start_coor)) &      # send is within the start and end coordinates) AND
            (blastn_df["sseqid"] == name_chr) &                                          # sseqid matches name_chr AND
            (blastn_df["sstrand"] == "plus")                                             # sstrand is "plus"
        )].copy()
        
        # Same but for the minus one, where the coor are inverted
        blastn_df = blastn_df[~(
            ((blastn_df["sstart"] <= end_coor) & (blastn_df["sstart"] >= start_coor)) |
            ((blastn_df["send"] >= start_coor) & (blastn_df["send"] <= end_coor)) &
            (blastn_df["sseqid"] == name_chr) &
            (blastn_df["sstrand"] == "minus"))].copy()
        
        
        # 3.4) Call bedops merge function
        blastn_df.sort_values(by=['qstart'], inplace=True)  # Sort the blastn_df by qstart. IMPORTANT for the bedops_merge function
        bedops_df = bedops_merge(input_df=blastn_df, path_folder=args.output_dir)

        # 3.5) Prepare the dict
        main_dict[name_id] = []

        # 3.5) Get the new coordinates
        for _, row in bedops_df.iterrows():
            # In the next steps it's important to add "-1" because if the bedops starts in "1" it means the start_coord will stay the same, so instead of adding +1 it should add +0. And we get that by adding -1
            new_start = start_coor + row['qstart'] - 1  # Get the new start coordinate by adding the start_coor and the qstart
            new_end = start_coor + row['qend'] - 1  # Get the new end coordinate by adding the start_coor and the qend. Important to be the "start" and not the "end"
            if abs(new_end - new_start) + 1 > 100:
                main_dict[name_id].append([name_chr, strand_seq, new_start, new_end])  # Append the new coordinates to the main_dict
            else:  # If the length is less than 100
                continue  # Continue to the next iteration
        
        # 3.6) Print the results (just for output info)
        print(f'\tFINISHED ==> {len(main_dict[name_id])} new sequences:')
        for seq in main_dict[name_id]:
            print(f'\t\t{seq[0]}:{seq[1]}:{seq[2]}-{seq[3]}')
        print('')  # Print a new line

    # 4) Save the dict as a JSON file
    with open(os.path.join(args.output_dir, "main_dict.json"), "w") as file:
        json.dump(main_dict, file, indent=4)  # Save the main_dict as a JSON file