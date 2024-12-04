# Define needed modules
import argparse
import json
import os
import pandas as pd
import subprocess

# ======================================================================
# Defining parse arguments
# ======================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Creating new data bases from JSON python dict')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the JSON file that contains the python tictionary.')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory output path to save all the data')
    parser.add_argument('-db', '--database', type=str, required=True, help='Path to the genome fasta file.')
    parser.add_argument('-ndb', '--neg_database', type=str, required=True, help='Path to the old negative elements that did not pass the filter')
    return parser.parse_args()
# ======================================================================
# Defining neeeded functions
# ======================================================================
def blastn_dic(path_input, path_output):
    cmd = f'makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}'
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def get_sequence(start_coor, end_coor, strand, chromosome, path_genome):
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

# ======================================================================
# Main function
# ======================================================================
if __name__ == '__main__':
    args = parse_arguments()  # Parse arguments

    # Read JSON file containing the dictionary
    with open(args.file, 'r') as file:
        data_dict = json.load(file)  # With this I have a dictionary python data type
        # Structured of the dictionary
            # element[0] = chromosome
            # element[1] = strand
            # element[2] = start coordinate
            # element[3] = end coordinate
            # element[4] = Accepted or Rejected

    # Read the negative database
    neg_database = pd.read_csv(args.neg_database, sep=',', header=0)
    # In the negative database
        # element['sseqid'] = chromosome
        # element['sstart'] = start coordinate
        # element['send'] = end coordinate
        # element['sstrand'] = strand
        # element['sseq'] = sequence

    # Let's start by creating the positive database
    positive_database = [element[0:4] for value in data_dict.values() for element in value if element[4] == 'Accepted']  # Save from chromosome to end coordinate
    negative_database = [element[0:4] for value in data_dict.values() for element in value if element[4] != 'Accepted']

    # Transform the database to a pandas dataframe
    positive_database = pd.DataFrame(positive_database, columns=['sseqid', 'sstrand', 'sstart', 'send'])
    negative_database = pd.DataFrame(negative_database, columns=['sseqid', 'sstrand', 'sstart', 'send'])

    # Reorder the columns to 'sseqid', 'sstart' , 'send', 'sstrand'
    positive_database = positive_database[['sseqid', 'sstart', 'send', 'sstrand']]
    negative_database = negative_database[['sseqid', 'sstart', 'send', 'sstrand']]

    # Take the sequence from the genome for both databases
    ## First let's create a blastn dictionary
    path_genome = args.database
    path_blast_dict_folder = os.path.join(args.output_dir, 'blastn_dict')
    os.makedirs(path_blast_dict_folder, exist_ok=True)
    path_blast_dict_file = os.path.join(path_blast_dict_folder, os.path.basename(path_genome))
    blastn_dic(path_genome, path_blast_dict_file)

    ## Now let's get the sequence for the positive database
    positive_database['sseq'] = positive_database.apply(
        lambda x: get_sequence(start_coor=x['sstart'], end_coor=x['send'], strand=x['sstrand'], chromosome=x['sseqid'], path_genome=path_blast_dict_file), axis=1
    )

    # And do the same for the negative database
    negative_database['sseq'] = negative_database.apply(
        lambda x: get_sequence(start_coor=x['sstart'], end_coor=x['send'], strand=x['sstrand'], chromosome=x['sseqid'], path_genome=path_blast_dict_file), axis=1
    )

    # Add the negative database to the old one
    negative_database = pd.concat([negative_database, neg_database], ignore_index=True, axis=0)

    # Make data types conversion just in case
    positive_database[['sstart', 'send']] = positive_database[['sstart', 'send']].apply(pd.to_numeric)
    negative_database[['sstart', 'send']] = negative_database[['sstart', 'send']].apply(pd.to_numeric)

    # Reorder data 
    positive_database.sort_values(by=['sseqid', 'sstart'], inplace=True)
    negative_database.sort_values(by=['sseqid', 'sstart'], inplace=True)

    # Save the databases
    positive_database.to_csv(os.path.join(args.output_dir, 'positive_database.csv'), sep=',', index=False, header=True)
    negative_database.to_csv(os.path.join(args.output_dir, 'negative_database.csv'), sep=',', index=False, header=True)