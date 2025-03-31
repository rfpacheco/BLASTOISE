import os
import pandas as pd
import subprocess
import json

# noinspection PyUnresolvedReferences
from modules.blaster import blastn_dic
# noinspection PyUnresolvedReferences
from extra.second_functions import get_sequence, bedops_merge, csv_to_fasta_creator, general_blastn_blaster


def coordinates_corrector(df, dict_path, folder_path, word_size, min_length):
    """
        Adjusts sequence coordinates using BLASTN and BEDOPS processing and saves the updated results.

        This function takes a DataFrame containing sequence information, utilizes BLASTN to process sequences,
        filters sequences based on coordinate overlap and strand validity, merges the results with BEDOPS,
        and computes updated start and end coordinates. The updated results are stored in a dictionary and saved
        as a JSON file.

        Arguments:
        df : DataFrame
            The input pandas DataFrame containing the sequence data with required columns: 'sseqid', 'sstrand',
            'sstart', 'send', and 'sseq'.
        dict_path : str
            The file path to the dictionary to be used in the BLASTN process.
        folder_path : str
            The folder path where intermediate and output files will be stored.
        word_size : int
            The word size parameter for the BLASTN alignment process.
        min_length : int
            The minimum required length of the sequences after adjusting their coordinates.

        Returns:
        str
            The file path of the saved JSON file containing the updated coordinates.

        Raises:
        TypeError
            If the types of input arguments do not match the specified types.
        ValueError
            If some of the required columns are missing from the DataFrame.
    """
    main_dict = {}
    for index, row in df.iterrows():
        # Prepare the query data
        name_id = f"{row['sseqid']}_{row['sstrand']}_{row['sstart']}-{row['send']}"  # Create the name_id with the original info of the seq.
        seq = row["sseq"]
        query = f"<(echo -e '>{name_id}\n{seq}')"  # Create the query with the name_id and the seq in a bash tmp file
        start_coor = row["sstart"]  # Get the start coordinate for later
        end_coor = row["send"]  # Get the end coordinate for later
        strand_seq = row["sstrand"]  # Get the strand of the sequence for later
        name_chr = row["sseqid"]  # Get the name of the chromosome for later
        print(f"Analyzing row {index + 1}/{df.shape[0]} with name_id {name_id}")

        # Call the blastn function
        blastn_df = general_blastn_blaster(query_path=query, dict_path=dict_path, word_size=word_size)

        # Filter the blastn_df
        # remove the row with the same `sstart` and `send` values than `start_coor` and `end_coor` in plus way
        blastn_df = blastn_df[~(
                ((blastn_df["sstart"] >= start_coor) & (blastn_df["sstart"] <= end_coor)) |  # `sstart` is within the start and end coordinates OR
                ((blastn_df["send"] <= end_coor) & (blastn_df["send"] >= start_coor)) &  # send is within the start and end coordinates) AND
                (blastn_df["sseqid"] == name_chr) &  # sseqid matches name_chr AND
                (blastn_df["sstrand"] == "plus")  # sstrand is "plus"
        )].copy()

        # Same but for the minus one, where the coor are inverted
        blastn_df = blastn_df[~(
                ((blastn_df["sstart"] <= end_coor) & (blastn_df["sstart"] >= start_coor)) |
                ((blastn_df["send"] >= start_coor) & (blastn_df["send"] <= end_coor)) &
                (blastn_df["sseqid"] == name_chr) &
                (blastn_df["sstrand"] == "minus"))].copy()

        # Call bedops merge function
        blastn_df.sort_values(by=["qstart"],
                              inplace=True)  # Sort the blastn_df by sstart. IMPORTANT for the bedops_merge function
        bedops_df = bedops_merge(input_df=blastn_df, path_folder=folder_path)

        # 3.5) Prepare the dict
        main_dict[name_id] = []

        # 3.5) Get the new coordinates
        for _, row2 in bedops_df.iterrows():
            # In the next steps it's important to add "-1" because if the bedops starts in "1" it means the start_coord will stay the same, so instead of adding +1 it should add +0. And we get that by adding -1
            new_start = start_coor + row2[
                "qstart"] - 1  # Get the new start coordinate by adding the start_coor and the qstart
            new_end = start_coor + row2[
                "qend"] - 1  # Get the new end coordinate by adding the start_coor and the qend. Important to be the "start" and not the "end"
            if abs(new_end - new_start) + 1 >= min_length:
                main_dict[name_id].append(
                    [name_chr, strand_seq, new_start, new_end])  # Append the new coordinates to the main_dict
            else:  # If the length is less than `min_length`
                continue  # Continue to the next iteration

        # Print the results (just for output info)
        print(f"\tFINISHED ==> {len(main_dict[name_id])} new sequences:")
        for seq in main_dict[name_id]:
            print(f"\t\t{seq[0]}:{seq[1]}:{seq[2]}-{seq[3]}")
        print("")  # Print a new line

    # 4) Save the dict as a JSON file
    with open(os.path.join(folder_path, "main_dict.json"), "w") as file:
        # noinspection PyTypeChecker
        json.dump(main_dict, file, indent=4)  # Save the main_dict as a JSON file

    return os.path.join(folder_path, "main_dict.json")  # return path

# ======================================================================
# noinspection DuplicatedCode
def json_sider_filter(json_file, folder_path, dict_path, word_size, evalue):
    """
        Filters and analyzes elements in a JSON file based on specific criteria, identifying whether they
        are accepted or rejected, and outputs the results to a new JSON file.

        Parameters:
        json_file: str
            The path to the JSON file containing the original data.
        folder_path: str
            The directory where the output JSON file will be saved.
        dict_path: str
            Path to the genome dictionary required for sequence retrieval and BLASTn processing.
        word_size: int
            The word size parameter for the BLASTn search.
        evalue: float
            The e-value threshold for filtering BLASTn hits.

        Returns:
        str
            The file path to the filtered JSON data stored in 'filtered_data.json'.

        Raises:
        FileNotFoundError
            If the JSON file or genome dictionary is not found.
        ValueError
            If the JSON file contains invalid or unexpected data.
    """
    with open(json_file, "r") as file:
        json_data = json.load(file)  # Loading json file into a python dict
        # In the dict for each element:
        # element[0] = chromosome
        # element[1] = strand
        # element[2] = start coordinate
        # element[3] = end coordinate
    print(f"Analyzing {len(json_data)} elements.")
    total_elements_counter = 0
    for key, value in json_data.items():
        for _, element in enumerate(value, start=0):
            total_elements_counter += 1
    print(f"\tTotal elements to analyze: {total_elements_counter}")

    # Loop through the dictionary to get the sequences
    accepted_elements_not_fragmented = 0
    accepted_elements_fragmented = 0
    rejected_elements_fragmented = 0
    for index, (key, value) in enumerate(json_data.items(), start=0):
        print("")
        print(f"Analyzing element {index + 1}/{len(json_data)} ==> {key}")
        if len(value) == 1:
            json_data[key][0].append('Accepted')  # If only one element, then it is accepted, since I don't want to check these, only the fragmented ones.
            print(f"\tAccepted ==> Not fragmented element")
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
                name_id = f"{key}_{i}"
                query = f"<(echo -e '>{name_id}\n{sequence}')"  # create bash tmp file

                # Run BLASTn
                blastn_df = general_blastn_blaster(query_path=query,
                                                   dict_path=dict_path,
                                                   word_size=word_size,
                                                   evalue=evalue)  # First evalue filter
                # Check BLASTn lines
                if not blastn_df.empty:
                    if blastn_df["sseqid"].nunique() >= 5:
                        json_data[key][i].append('Accepted')
                        print(f"\tAccepted ==> Fragmented element {i + 1} of {len(value)}")
                        accepted_elements_fragmented += 1
                    else:  # If not accepted
                        json_data[key][i].append('Rejected')
                        print(f"\tRejected ==> Fragmented element {i + 1} of {len(value)}")
                        rejected_elements_fragmented += 1
                else:  # If empty, then it is rejected
                    json_data[key][i].append('Rejected')
                    print(f"\tRejected ==> Fragmented element {i + 1} of {len(value)}")
                    rejected_elements_fragmented += 1

    # Save the data
    with open(os.path.join(folder_path, "filtered_data.json"), 'w') as file:
        # noinspection PyTypeChecker
        json.dump(json_data, file, indent=4)

    # Print summary
    total_elements_accepted = accepted_elements_not_fragmented + accepted_elements_fragmented  # fragmented or not fragmented
    print("")
    print("=" * 50)
    print(f"- Total elements analyzed: {total_elements_counter} (with fragmentation) from {len(json_data)} original elements.")
    print(f"\t- Total accepted fragmented elements: {accepted_elements_fragmented}/{total_elements_counter}")
    print(f"\t\t- {accepted_elements_fragmented / total_elements_counter * 100:.2f} % of the total fragmented elements elements.")
    print(f"\t- Total rejected fragmented elements: {rejected_elements_fragmented}/{total_elements_counter}")
    print(f"\t\t- {rejected_elements_fragmented / total_elements_counter * 100:.2f} % of the total fragmented elements elements.")
    print(f"\t- Total accepted not fragmented elements: {accepted_elements_not_fragmented}/{total_elements_counter}")
    print(f"\t\t- {accepted_elements_not_fragmented / total_elements_counter * 100:.2f} % of the total not fragmented elements elements.")
    print("")
    print(f"- Total accepted elements: {total_elements_accepted} from {total_elements_counter} elements.")

    return os.path.join(folder_path, "filtered_data.json")

# ======================================================================
# noinspection DuplicatedCode
def sider_json_to_csv(json_file, folder_path, dict_path, recaught_file, recaught_threshold, word_size, perc_identity):
    """
        Transform data from a JSON file into two CSV files (positive and negative databases),
        process accepted and rejected elements, extract sequences, and handle recaptured data.

        This function reads a structured JSON file containing chromosomal data, processes the
        contents to extract positive and negative databases, and optionally includes recaptured
        data based on BLASTn results. The positive and negative databases are saved as separate
        CSV files. Additionally, appropriate file permissions are modified.

        Attributes:
            json_file (str): Path to the input JSON file.
            folder_path (str): Path to the folder to store intermediate outputs, like the negative
                database in FASTA format.
            dict_path (str): Path to the genome dictionary for sequence retrieval.
            recaught_file (str): Path to the file containing data used for recapturing sequences.
            recaught_threshold (float): E-value threshold for filtering recaptured data.
            word_size (int): Word size parameter used in BLASTn analysis for recapturing.
            perc_identity (float): Percentage identity threshold for BLASTn alignment during
                recapturing.

        Parameters:
            json_file: The path to the JSON file that contains the data to be processed.
            folder_path: The destination folder for intermediate generated files.
            dict_path: The path to the genome dictionary file needed for sequence generation.
            recaught_file: The input file to capture sequences that may have been rejected or
                missed.
            recaught_threshold: Specifies the threshold for e-value in filtering recaught
                results.
            word_size: Controls the word size used during recapturing using BLASTn.
            perc_identity: Minimum identity percentage threshold to consider matches during
                recapturing.

        Raises:
            FileNotFoundError: If any input file paths do not exist.
            ValueError: If the input data or parameters do not meet expected criteria.
            SubprocessError: If external subprocess calls (e.g., for setting file permissions)
                fail.

        Returns:
            None

        Notes:
            The function performs multiple tasks including:
             - Loading the JSON input file into a pandas DataFrame for processing.
             - Sorting data columns for consistency and clarity.
             - Retrieving sequences for positive and negative databases using genomic
               coordinates.
             - Generating a BLASTn database and running recapture operations if needed.
             - Saving final processed data (positive and negative) into separate CSV files.
             - Modifying file permissions to ensure appropriate access.

            External utility functions and subprocesses (e.g., BLASTn) are required to run this
            function successfully.
    """
    print("1. Loading JSON file...")
    with open(json_file, "r") as file:
        data_dict = json.load(file)  # With this I have a dictionary python data type
        # Structured of the dictionary
        # element[0] = chromosome
        # element[1] = strand
        # element[2] = start coordinate
        # element[3] = end coordinate
        # element[4] = Accepted or Rejected
    print("\t DONE")

    # Let's start by creating the positive database
    print("2. Getting accepted or rejected elements...")
    positive_database = [element[0:4] for value in data_dict.values() for element in value if
                         element[4] == "Accepted"]  # Save from chromosome to end coordinate
    negative_database = [element[0:4] for value in data_dict.values() for element in value if element[4] != "Accepted"]
    print("\t DONE")

    # Transform the database to a pandas dataframe
    print("3. Transforming to data frame...")
    positive_database = pd.DataFrame(positive_database, columns=["sseqid", "sstrand", "sstart", "send"])
    negative_database = pd.DataFrame(negative_database, columns=["sseqid", "sstrand", "sstart", "send"])
    print("\t DONE")

    # Reorder the columns to 'sseqid', 'sstart' , 'send', 'sstrand'
    print("4. Getting correct order...")
    positive_database = positive_database[["sseqid", "sstart", "send", "sstrand"]]
    negative_database = negative_database[["sseqid", "sstart", "send", "sstrand"]]
    print("\t DONE")

    ## Now let's get the sequence for the positive database
    print("5. Getting positive sequences...")
    positive_database['sseq'] = positive_database.apply(
        lambda x: get_sequence(start_coor=x['sstart'], end_coor=x['send'], strand=x['sstrand'],
                               chromosome=x['sseqid'],
                               path_genome=dict_path), axis=1
    )
    print("\t DONE")

    # And do the same for the negative database
    print("6. Getting negative sequences...")
    negative_database["sseq"] = negative_database.apply(
        lambda x: get_sequence(start_coor=x["sstart"], end_coor=x["send"], strand=x["sstrand"],
                               chromosome=x["sseqid"],
                               path_genome=dict_path), axis=1
    )
    print("\t DONE")

    # Make data types conversion just in case
    print("7. Transforming data...")
    positive_database[["sstart", "send"]] = positive_database[["sstart", "send"]].apply(pd.to_numeric)
    negative_database[["sstart", "send"]] = negative_database[["sstart", "send"]].apply(pd.to_numeric)
    print("\t DONE")

    # Recaught data
    print("8. Recaught lost data data...")
    ## Create a fasta file from the pandas.DataFrame from the negative DB to make a blastN DB
    fasta_file_path = os.path.join(folder_path, "negative_database.fasta")
    csv_to_fasta_creator(negative_database, fasta_file_path)
    ## Create the BALSTn dict
    blastn_dic(path_input=fasta_file_path, path_output=dict_path)
    ## Recaught the data
    caught_data = general_blastn_blaster(
        query_path=recaught_file,
        dict_path=dict_path,
        perc_identity=perc_identity,
        word_size=word_size
    )

    ## If caught data is with data
    if not caught_data.empty:
        # Remove ones with an evalue <= 10**-3
        caught_data = caught_data[caught_data["evalue"] <= recaught_threshold].sort_values(by=["evalue"])
        print("")
        print("*" * 50)
        print(f"\nRecaught data: {caught_data.shape[0]} elements")

        # Create a column with the number in "sseqid"
        caught_data["index"] = caught_data["sseqid"].str.extract(r"_(\d+)_")
        caught_data["index"] = pd.to_numeric(caught_data["index"])

        # Get a list with the index column
        index_list = caught_data["index"].sort_values().unique().tolist()

        # Extract sequences from the 'negative_database'
        recaught_in_neg = negative_database[negative_database.index.isin(index_list)]

        # Join yes_data and recaught_in_neg
        accepted_elements = pd.concat([positive_database, recaught_in_neg], axis=0, ignore_index=True)
        accepted_elements.sort_values(by=['sseqid', 'sstart'], inplace=True)

        # Remove recaught_in_neg from no data
        rejected_elements = pd.concat([negative_database, recaught_in_neg]).drop_duplicates(keep=False)

        # Print results:
        print(f"\n\t - Accepted data + recaught: {accepted_elements.shape[0]} elements")
        print(f"\t - Rejected data - recaught: {rejected_elements.shape[0]} elements")

    else:  # if it's empty
        accepted_elements = positive_database
        rejected_elements = negative_database
        print("\n\t - No recaught data")


    # Reorder data
    print("9. Reordering data...")
    accepted_elements.sort_values(by=["sstrand", "sseqid", "sstart"], inplace=True)
    rejected_elements.sort_values(by=["sstrand", "sseqid", "sstart"], inplace=True)
    print("\t DONE")

    # Save the databases
    positive_path = os.path.join(os.path.dirname(folder_path), "siders_df.csv")
    negative_path = os.path.join(os.path.dirname(folder_path), "non_siders_df.csv")
    accepted_elements.to_csv(positive_path, sep=',', index=False, header=True)
    rejected_elements.to_csv(negative_path, sep=',', index=False, header=True)
    print("")
    print(f"Positive DataBase saved at: {positive_path}")
    print(f"Negative DataBase saved at: {negative_path}")

    # Add right to groups and users
    subprocess.run(["chmod", "-R", "a+w", positive_path], check=True)
    subprocess.run(["chmod", "-R", "a+w", negative_path], check=True)




