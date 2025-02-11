import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json

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

# noinspection DuplicatedCode
def simple_blastn_blaster(query_path, dict_path):
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


# noinspection DuplicatedCode
def json_blastn_blaster(query, path_genome, evalue):
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

    return data_df

# ======================================================================================================================
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

# ======================================================================================================================
# noinspection DuplicatedCode
def coordinates_corrector(df, dict_path, folder_path):
    main_dict = {}
    for index, row in df.iterrows():
        # 3.1) Prepare the query data
        name_id = f"{row['sseqid']}_{row['sstrand']}_{row['sstart']}-{row['send']}"  # Create the name_id with the original info of the seq.
        seq = row['sseq']
        query = f"<(echo -e '>{name_id}\n{seq}')"  # Create the query with the name_id and the seq in a bash tmp file
        start_coor = row['sstart']  # Get the start coordinate for later
        end_coor = row['send']  # Get the end coordinate for later
        strand_seq = row['sstrand']  # Get the strand of the sequence for later
        name_chr = row['sseqid']  # Get the name of the chromosome for later
        print(f'Analyzing row {index + 1}/{df.shape[0]} with name_id {name_id}')

        # 3.2) Call the blastn function
        blastn_df = simple_blastn_blaster(query_path=query, dict_path=dict_path)  # Call the blastn function  # TODO: simplify blastn function

        # 3.3) Filter the blastn_df
        # remove the row with the same sstart and send values than start_coor and end_coor in plus way
        blastn_df = blastn_df[~(
                ((blastn_df["sstart"] >= start_coor) & (
                        blastn_df["sstart"] <= end_coor)) |  # (sstart is within the start and end coordinates OR
                ((blastn_df["send"] <= end_coor) & (
                        blastn_df["send"] >= start_coor)) &  # send is within the start and end coordinates) AND
                (blastn_df["sseqid"] == name_chr) &  # sseqid matches name_chr AND
                (blastn_df["sstrand"] == "plus")  # sstrand is "plus"
        )].copy()

        # Same but for the minus one, where the coor are inverted
        blastn_df = blastn_df[~(
                ((blastn_df["sstart"] <= end_coor) & (blastn_df["sstart"] >= start_coor)) |
                ((blastn_df["send"] >= start_coor) & (blastn_df["send"] <= end_coor)) &
                (blastn_df["sseqid"] == name_chr) &
                (blastn_df["sstrand"] == "minus"))].copy()

        # 3.4) Call bedops merge function
        blastn_df.sort_values(by=['qstart'],
                              inplace=True)  # Sort the blastn_df by sstart. IMPORTANT for the bedops_merge function
        bedops_df = bedops_merge(input_df=blastn_df, path_folder=folder_path)

        # 3.5) Prepare the dict
        main_dict[name_id] = []

        # 3.5) Get the new coordinates
        for _, row2 in bedops_df.iterrows():
            # In the next steps it's important to add "-1" because if the bedops starts in "1" it means the start_coord will stay the same, so instead of adding +1 it should add +0. And we get that by adding -1
            new_start = start_coor + row2[
                'qstart'] - 1  # Get the new start coordinate by adding the start_coor and the qstart
            new_end = start_coor + row2[
                'qend'] - 1  # Get the new end coordinate by adding the start_coor and the qend. Important to be the "start" and not the "end"
            if abs(new_end - new_start) + 1 > 100:
                main_dict[name_id].append(
                    [name_chr, strand_seq, new_start, new_end])  # Append the new coordinates to the main_dict
            else:  # If the length is less than 100
                continue  # Continue to the next iteration

        # 3.6) Print the results (just for output info)
        print(f'\tFINISHED ==> {len(main_dict[name_id])} new sequences:')
        for seq in main_dict[name_id]:
            print(f'\t\t{seq[0]}:{seq[1]}:{seq[2]}-{seq[3]}')
        print('')  # Print a new line

    # 4) Save the dict as a JSON file
    with open(os.path.join(folder_path, "main_dict.json"), "w") as file:
        json.dump(main_dict, file, indent=4)  # Save the main_dict as a JSON file

    return os.path.join(folder_path, "main_dict.json")  # return path

# ======================================================================================================================
# noinspection DuplicatedCode
def get_sequence(start_coor, end_coor, strand, chromosome, path_genome):
    cmd = f'blastdbcmd -db {path_genome} -entry {chromosome} -range {start_coor}-{end_coor} -strand {strand} -outfmt %s'
    sequence = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')
    sequence = sequence.stdout.strip()
    return sequence

# ======================================================================================================================
# noinspection DuplicatedCode
def json_sider_filter(json_file, folder_path, dict_path):
    with open(json_file, 'r') as file:
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

    # Loop through the dictionary to get the sequences
    accepted_elements_not_fragmented = 0
    accepted_elements_fragmented = 0
    rejected_elements_fragmented = 0
    for index, (key, value) in enumerate(json_data.items(), start=0):
        print('')
        print(f'Analyzing element {index + 1}/{len(json_data)} ==> {key}')
        if len(value) == 1:
            json_data[key][0].append(
                'Accepted')  # If only one element, then it is accepted, since I don't want to check this ones, only the fragmented ones.
            print(f'\tAccepted ==> Not fragmented element')
            accepted_elements_not_fragmented += 1
            continue  # Skip to the next iteration
        else:  # If more than one element
            for i, element in enumerate(value, start=0):
                sequence = get_sequence(start_coor=element[2],
                                        end_coor=element[3],
                                        strand=element[1],
                                        chromosome=element[0],
                                        path_genome=dict_path)
                # Make a BLASTn with this sequence with the filter:
                ## Prepare data
                name_id = f'{key}_{i}'
                query = f"<(echo -e '>{name_id}\n{sequence}')"  # create bash tmp file
                evalue = 1.0E-09

                # Run BLASTn
                blastn_df = json_blastn_blaster(query=query,
                                                path_genome=dict_path,
                                                evalue=evalue)
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
                else:  # If empty, then it is rejected
                    json_data[key][i].append('Rejected')
                    print(f'\tRejected ==> Fragmented element {i + 1} of {len(value)}')
                    rejected_elements_fragmented += 1

    # Save the data
    with open(os.path.join(folder_path, 'filtered_data.json'), 'w') as file:
        json.dump(json_data, file, indent=4)

    # Print summary
    total_elements_accepted = accepted_elements_not_fragmented + accepted_elements_fragmented  # fragmented or not fragmented
    print('')
    print('=' * 50)
    print(f'- Total elements analyzed: {total_elements_counter} from {len(json_data)} original elements.')
    print(f'\t- Total accepted fragmented elements: {accepted_elements_fragmented}/{total_elements_counter}')
    print(
        f'\t\t- {accepted_elements_fragmented / total_elements_counter * 100:.2f} % of the total fragmented elements elements.')
    print(f'\t- Total rejected fragmented elements: {rejected_elements_fragmented}/{total_elements_counter}')
    print(
        f'\t\t- {rejected_elements_fragmented / total_elements_counter * 100:.2f} % of the total fragmented elements elements.')
    print(f'\t- Total accepted not fragmented elements: {accepted_elements_not_fragmented}/{total_elements_counter}')
    print(
        f'\t\t- {accepted_elements_not_fragmented / total_elements_counter * 100:.2f} % of the total not fragmented elements elements.')
    print('')
    print(f'- Total accepted elements: {total_elements_accepted} from the original {len(json_data)} elements.')
    return os.path.join(folder_path, 'filtered_data.json')

# ======================================================================================================================
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
def sider_json_to_csv(json_file, folder_path, dict_path, neg_db_df):
    print("1. Loading JSON file...")
    with open(json_file, 'r') as file:
        data_dict = json.load(file)  # With this I have a dictionary python data type
        # Structured of the dictionary
        # element[0] = chromosome
        # element[1] = strand
        # element[2] = start coordinate
        # element[3] = end coordinate
        # element[4] = Accepted or Rejected
    print("\t DONE")

    # Read the negative database
    print("2. Loading negative database...")
    # neg_database = pd.read_csv(neg_db_path, sep=',', header=0)  # TODO: remove if not needed
    neg_database = neg_db_df
    # In the negative database
    # element['sseqid'] = chromosome
    # element['sstart'] = start coordinate
    # element['send'] = end coordinate
    # element['sstrand'] = strand
    # element['sseq'] = sequence
    print("\t DONE")

    # Let's start by creating the positive database
    print("3. Getting accepted or rejected elements...")
    positive_database = [element[0:4] for value in data_dict.values() for element in value if
                         element[4] == 'Accepted']  # Save from chromosome to end coordinate
    negative_database = [element[0:4] for value in data_dict.values() for element in value if element[4] != 'Accepted']
    print("\t DONE")

    # Transform the database to a pandas dataframe
    print("4. Transforming to data frame...")
    positive_database = pd.DataFrame(positive_database, columns=['sseqid', 'sstrand', 'sstart', 'send'])
    negative_database = pd.DataFrame(negative_database, columns=['sseqid', 'sstrand', 'sstart', 'send'])
    print("\t DONE")

    # Reorder the columns to 'sseqid', 'sstart' , 'send', 'sstrand'
    print("5. Getting correct order...")
    positive_database = positive_database[['sseqid', 'sstart', 'send', 'sstrand']]
    negative_database = negative_database[['sseqid', 'sstart', 'send', 'sstrand']]
    print("\t DONE")

    ## Now let's get the sequence for the positive database
    print("6. Getting positive sequences...")
    positive_database['sseq'] = positive_database.apply(
        lambda x: get_sequence_json_to_csv(start_coor=x['sstart'], end_coor=x['send'], strand=x['sstrand'],
                                           chromosome=x['sseqid'],
                                           path_genome=dict_path), axis=1
    )
    print("\t DONE")

    # And do the same for the negative database
    print("7. Getting negative sequences...")
    negative_database['sseq'] = negative_database.apply(
        lambda x: get_sequence_json_to_csv(start_coor=x['sstart'], end_coor=x['send'], strand=x['sstrand'],
                                           chromosome=x['sseqid'],
                                           path_genome=dict_path), axis=1
    )
    print("\t DONE")

    # Add the negative database to the old one
    print("8. Merging negative databases...")
    negative_database = pd.concat([negative_database, neg_database], ignore_index=True, axis=0)
    print("\t DONE")

    # Make data types conversion just in case
    print("9. Transforming data...")
    positive_database[['sstart', 'send']] = positive_database[['sstart', 'send']].apply(pd.to_numeric)
    negative_database[['sstart', 'send']] = negative_database[['sstart', 'send']].apply(pd.to_numeric)
    print("\t DONE")

    # Reorder data
    print("10. Reordering data...")
    positive_database.sort_values(by=['sseqid', 'sstart'], inplace=True)
    negative_database.sort_values(by=['sseqid', 'sstart'], inplace=True)
    print("\t DONE")

    # Save the databases
    positive_database.to_csv(os.path.join(folder_path, 'positive_database.csv'), sep=',', index=False, header=True)
    negative_database.to_csv(os.path.join(folder_path, 'negative_database.csv'), sep=',', index=False, header=True)
    print("Program finished!")
