# Importing needed modules
import argparse
import json
import subprocess
import os
import pandas as pd

 # ======================================================================
 # Defining parse arguments
 # ======================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter elements based on chromosome count and e-value.")
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the input JSON file that contains the dictionary.')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save the output data.')
    parser.add_argument('-db', '--database', type=str, required=True, help='Path to the genome fasta file.')
    return parser.parse_args()

# ======================================================================
# Defining neeeded functions
# ======================================================================
def blastn_dic(path_input, path_output):
    # "parse_seqids" is used to keep the sequence ID in the output.
    cmd = f"makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}"
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def get_sequence(start_coor, end_coor, strand, chromosome, path_genome):
    cmd = f'blastdbcmd -db {path_genome} -entry {chromosome} -range {start_coor}-{end_coor} -strand {strand} -outfmt %s'
    sequence = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')
    sequence = sequence.stdout.strip()
    return sequence

def blastn_blaster(query, path_genome, evalue):
    cmd = (
        f'blastn -word_size 15 '
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

    return data_df  # I only need to return the data, no need to transform it into a pandas dataframe since I want to check the number of lines with the line separator '\n'
    
# ======================================================================
# Main function
# ======================================================================
if __name__ == '__main__':
    args = parse_arguments()  # Parse arguments
    # Read JSON file containing the dictionary
    with open(args.file, 'r') as file:
        json_data = json.load(file)  # Loading json file into a python dict
        # In the dict for each element:
            # element[0] = chromosome
            # element[1] = strand
            # element[2] = start coordinate
            # element[3] = end coordinate
    print(f'Analyzing {len(json_data)} elements.')
    total_elements_counter = 0
    for key, value in json_data.items():
        for _, element in enumerate(value, start=0):
            total_elements_counter += 1
    print(f'\tTotal elements to analyze: {total_elements_counter}')
            
    # Create a BLASTn dictionary
    path_blast_folder = os.path.join(args.output_dir, 'blastn_dict')
    path_blast_dict_file = os.path.join(path_blast_folder, os.path.basename(args.database))
    os.makedirs(path_blast_folder, exist_ok=True)
    blastn_dic(path_input=args.database, path_output=path_blast_dict_file)

    # Loop through the dictionary to get the sequences
    accepted_elements_not_fragmented = 0
    accepted_elements_fragmented = 0
    rejected_elements_fragmented = 0
    for index, (key, value) in enumerate(json_data.items(), start=0):
        print('')
        print(f'Analyzing element {index + 1}/{len(json_data)} ==> {key}')
        if len(value) == 1:
            json_data[key][0].append('Accepted')  # If only one element, then it is accepted, since I don't want to check this ones, only the fragmented ones.
            print(f'\tAccepted ==> Not fragmented element')
            accepted_elements_not_fragmented += 1
            continue  # Skip to the next iteration
        else:  # If more than one element
            for i, element in enumerate(value, start=0):
                sequence = get_sequence(start_coor=element[2], 
                                        end_coor=element[3], 
                                        strand=element[1], 
                                        chromosome=element[0], 
                                        path_genome=path_blast_dict_file)
                # Make a BLASTn with this sequence with the filter:
                ## Prepare data
                name_id = f'{key}_{i}'
                query = f"<(echo -e '>{name_id}\n{sequence}')"  # create bash tmp file
                evalue = 1.0E-09

                # Run BLASTn
                blastn_df = blastn_blaster(query=query, path_genome=path_blast_dict_file, evalue=evalue)  # pandas dataframe 

                # Check BLASTn lines
                if not blastn_df.empty:
                    if blastn_df['sseqid'].nunique() >= 5:
                        json_data[key][i].append('Accepted')
                        print(f'\tAccepted ==> Fragmented element {i + 1} of {len(value)}')
                        accepted_elements_fragmented += 1
                    else:  # If not accepted
                        json_data[key][i].append('Rejected')
                        print(f'\tRejected ==> Fragmented element {i + 1} of {len(value)}')
                        rejected_elements_fragmented += 1
                else: # If empty, then it is rejected
                    json_data[key][i].append('Rejected')
                    print(f'\tRejected ==> Fragmented element {i + 1} of {len(value)}')
                    rejected_elements_fragmented += 1
    
    # Save the data
    with open(os.path.join(args.output_dir, 'filtered_data.json'), 'w') as file:
        json.dump(json_data, file, indent=4)

    # Print summary
    total_elements_accepted = accepted_elements_not_fragmented + accepted_elements_fragmented  # fragmented or not fragmented
    print('')
    print('='*50)
    print(f'- Total elements analyzed: {total_elements_counter} from {len(json_data)} original elements.')
    print(f'\t- Total accepted fragmented elements: {accepted_elements_fragmented}/{total_elements_counter}')
    print(f'\t\t- {accepted_elements_fragmented/total_elements_counter * 100:.2f} % of the total fragmented elements elements.')
    print(f'\t- Total rejected fragmented elements: {rejected_elements_fragmented}/{total_elements_counter}')
    print(f'\t\t- {rejected_elements_fragmented/total_elements_counter * 100:.2f} % of the total fragmented elements elements.')
    print(f'\t- Total accepted not fragmented elements: {accepted_elements_not_fragmented}/{total_elements_counter}')
    print(f'\t\t- {accepted_elements_not_fragmented/total_elements_counter * 100:.2f} % of the total not fragmented elements elements.')
    print('')
    print(f'- Total accepted elements: {total_elements_accepted} from the original {len(json_data)} elements.')