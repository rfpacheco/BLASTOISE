import os
import pandas as pd
import subprocess
import time
import shutil
import logging
# from pathlib import Path
from datetime import datetime

from modules.aesthetics import boxymcboxface  # Some aesthetics function
from modules.identifiers import genome_specific_chromosome_main
from modules.filters import global_filters_main
# from modules.files_manager import columns_to_numeric  # TODO: check if remove
from modules.compare import compare_main

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_dic(path_input, path_output):
    """
    Creation af a BLAST database of our whole genome. It uses the BLAST :sup:`R` command line, see BLAST
    `Command Line Application User Manual`_ for more information.


    The generation of the proper database will be placed in the directory where ``path_input`` is.
    It is recommended to use a dedicated folder to this FASTA file so the database is written next to it.

    :param path_input: path to a FASTA file.
    :type path_input: string

    :param path_output: path to the output folder where the BLAST database will be created.
    :type path_output: string

    :return: a BLAST database.
    :rtype: Multiples files (**.nhr**, **.nin**, **.nog**, **.nsd**, **.nsi** and **.nsq** extensions)
    """

    # Remember is "path.input.dic_path" for "argparse".
    try:
        # "parse_seqids" is used to keep the sequence ID in the output.
        cmd = f"makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        logging.error(f"Error: Blast Dictionary couldn't be created: {e}", exc_info=True)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_blaster(query_path, dict_path, perc_identity, word_size=15):
    """
    Executes a BLASTN search on the provided query sequence against a nucleotide database.

    Parameters:
        query_path (str): Path to the query file containing the sequence to be searched.
        dict_path (str): Path to the nucleotide database against which the query sequence is to be searched.
        perc_identity (float): Percent identity threshold for reporting matches.
        word_size (int, optional): Word size for the BLASTN search algorithm. Default is 15.

    Returns:
        pandas.DataFrame: DataFrame containing the BLASTN search results with columns:
            qseqid: Query sequence ID.
            sseqid: Subject sequence ID.
            pident: Percentage of identical matches.
            length: Alignment length.
            qstart: Start of alignment in query.
            qend: End of alignment in query.
            sstart: Start of alignment in subject.
            send: End of alignment in subject.
            evalue: Expectation value.
            bitscore: Bit score.
            qlen: Length of query sequence.
            slen: Length of subject sequence.
            sstrand: Strand of the subject sequence.
            sseq: Aligned part of the subject sequence.
    """

    cmd = "blastn -word_size " + str(word_size) + " -query " \
        + query_path + " -db " \
        + dict_path \
        + " -perc_identity " + str(perc_identity) \
        + " -outfmt '10 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen sstrand sseq'"
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)  # Important the E value
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    data.columns = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "sstrand", "sseq"]
    return data
    

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def repetitive_blaster(data_input, genome_fasta, folder_path, numbering, start_time, identity_1, tic_start, word_size, min_length, extend_number, coincidence_data=None):
    """
    Processes and filters genomic data by grouping, analyzing chromosomes, and comparing results with previous runs.

    Parameters:
    data_input : pandas.DataFrame
        The initial dataset containing genomic data.
    genome_fasta : str
        Path to the reference genome file in FASTA format.
    folder_path : str
        Directory path where results will be stored.
    numbering : int
        Run identifier number.
    start_time : str
        Starting timestamp of the run.
    identity_1 : float
        Identity threshold for genomic sequence analysis.
    tic_start : float
        Starting time as recorded by time.perf_counter() for performance logging.
    word_size : int
        Word size for sequence alignment.
    min_length : int
        Minimum length threshold for sequences.
    extend_number : int
        Number of bases to extend during analysis.
    coincidence_data : pandas.DataFrame, optional
        Data from a previous run for comparison.
    """

    # Call the aesthetics function RUN identifier.
    boxymcboxface("RUN " + str(numbering))
    tic_main = time.perf_counter()  # Start the timer

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    # First let's order the data by "sseqid", "sstrand", "sstart".
    data_ordered = data_input.sort_values(by=["sseqid", "sstrand", "sstart"])

    # Now let's group the data by "sseqid". We'll have a pandas groupby object.
    data_grouped = data_ordered.groupby("sseqid")
    toc = time.perf_counter()
    print("")
    print(f"1. Initial data:\n",
          f"\t- Data row length: {data_input.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Now let's call  `genome_specific_chromosome_main` for each chromosome_ID in the data using the groupby object.
    terminal_width = shutil.get_terminal_size().columns  # Get the terminal width

    print("")
    print(f"2. Individual searching and cleaning:")
    whole_group = pd.DataFrame()  # This will be the final data frame for each chromosome
    for _, (chromosome, group) in enumerate(data_grouped):
        tic = time.perf_counter()
        now_time = datetime.now()
        formatted_now_time = now_time.strftime("%Y %B %d at %H:%M")
        print("")
        print(f"{' ' * 7}{'-' * 74}")
        print(f"\t- {chromosome}:")
        start_time_text = f"Program started: {start_time}"
        end_time_text = f"Program time now: {formatted_now_time}"
        run_text = f"RUN {numbering}"
        print(f"{run_text:>{terminal_width}}")
        print(f"{start_time_text:>{terminal_width}}")
        print(f"{end_time_text:>{terminal_width}}")

        data = genome_specific_chromosome_main(data_input=group,
                                               chromosome_ID=chromosome,
                                               main_folder_path=folder_path,
                                               genome_fasta=genome_fasta,
                                               identity_1=identity_1,
                                               run_phase=numbering,
                                               coincidence_data=coincidence_data,
                                               word_size=word_size,
                                               min_length=min_length,
                                               extend_number=extend_number)
        toc = time.perf_counter()
        print("")
        print(f"\t\t- Data row length: {len(data)}\n",  # Not .shape[0] in case the data is empty
              f"\t\t- Execution time: {toc - tic:0.2f} seconds")
        whole_group = pd.concat([whole_group, data])
    print(f"{' ' * 7}{'-' * 74}")
    print("")
    print(f"\t- Blast data row length: {whole_group.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds\n")
    # -----------------------------------------------------------------------------   
    tic = time.perf_counter()
    print("")
    print(f"4. Global filtering:")
    whole_group_filtered = global_filters_main(data_input=whole_group,
                                               genome_fasta=genome_fasta,
                                               writing_path=folder_path,
                                               min_length=min_length)
    toc = time.perf_counter()
    print("")
    print(f"\t- Data row length: {whole_group_filtered.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")

    # -----------------------------------------------------------------------------
    ## Save the RUN
    runs_folder = os.path.join(folder_path, "RUNS")  # Creates the folder for the RUNS
    os.makedirs(runs_folder, exist_ok=True)  # Creates the folder for the RUNS
    run_saver_path = os.path.join(runs_folder, "run_" + str(numbering) + ".csv")  # Path to save the RUN
    whole_group_filtered.to_csv(run_saver_path, sep=",", header=True, index=False)  # Saves the RUN
    # -----------------------------------------------------------------------------
    # Compare part
    # Prepare folders and path
    comparison_folder = os.path.join(folder_path, "comparison")
    os.makedirs(comparison_folder, exist_ok=True)

    print("")
    print(f"5. Comparison vs Old Run:")
    if coincidence_data is not None:
        print("")
        print(f"\t- Last Run data:\n",
              f"\t\t- Coincidence data row length: {coincidence_data.shape[0]}\n",
              f"\t\t- New data row length: {data_input.shape[0]}")
        # This part is important to campere with the last run whole data "whoel_group_fildered", and not only the "new_data" subset.
        data_input = pd.concat([coincidence_data, data_input], ignore_index=True)
        data_input.sort_values(by=["sseqid", "sstrand", "sstart"], inplace=True)  # Sort the data frame by the start coordinate
        print(f"\t\t- Total data row length: {data_input.shape[0]}")
    else:  # when coincidence_data == None
        print(f"\t- Last Run data:\n",
              f"\t\t- First run row length: {data_input.shape[0]}")

    tic = time.perf_counter()
    print("")
    print(f"\t- Results in this RUN:")
    coincidence_data, new_data, old_data_exclusive = compare_main(last_df=whole_group_filtered,
                                                                  old_df=data_input,
                                                                  folder_path=comparison_folder,
                                                                  genome_fasta=genome_fasta)
    toc = time.perf_counter()


    print("")
    print(f"\t\t- Coincidence data row length: {coincidence_data.shape[0]}\n",
          f"\t\t- New data row length: {new_data.shape[0]}\n",
          f"\t\t- Previous data row length: {old_data_exclusive.shape[0]}")

    old_data_exclusive_less_than_100 = None

    if not old_data_exclusive.empty and (old_data_exclusive["length"] < 100).sum() > 0:  # If there are sequences less than 100 bp. The sum of TRUE (for < 100) has to be > 0
        old_data_exclusive_less_than_100 = old_data_exclusive[old_data_exclusive["length"] < 100]
        old_data_exclusive = old_data_exclusive[old_data_exclusive["length"] > 100]
    else:
        pass

    if old_data_exclusive_less_than_100 is not None: # If old_data_exclusive_less_than_100 exists
        new_data_and_old = pd.concat([new_data, old_data_exclusive_less_than_100], ignore_index=True)
        new_data_and_old.sort_values(by=["sseqid", "sstrand", "sstart"], inplace=True)
        print('\t' * 3 + f"- Less than 100 bp: {old_data_exclusive_less_than_100.shape[0]}")
        print('\t' * 3 + f"- New data + less than 100: {new_data_and_old.shape[0]}")
        print('\t' * 3 + f"- Previous data: {old_data_exclusive.shape[0]}")
    else:
        new_data_and_old = new_data

    # Join coincidence_data with old_data_exclusive
    if not coincidence_data.empty and not old_data_exclusive.empty:
        coincidence_data = pd.concat([coincidence_data, old_data_exclusive], ignore_index=True)
        coincidence_data.sort_values(by=["sseqid", "sstrand", "sstart"], inplace=True)
        print(f"\t\t- Coincidence data + Previous data: {coincidence_data.shape[0]}")
    else:
        pass
    print(f"\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Stopping part
    stopping_folder = os.path.join(folder_path, "stopping")
    os.makedirs(stopping_folder, exist_ok=True)

    if new_data.shape[0] == 0:
        coincidence_data = coincidence_data[['sseqid', 'length', 'sstart', 'send', 'sstrand', 'sseq']].copy()  # #Take only needed columns:
        
        # Make it so 'sstart' is always < than 'send'
        coincidence_data.loc[
            coincidence_data[coincidence_data['sstart'] > coincidence_data['send']].index,
            ['sstart', 'send']
        ] = coincidence_data.loc[
            coincidence_data[coincidence_data['sstart'] > coincidence_data['send']].index,
            ['send', 'sstart']
        ].values

        coincidence_data.to_csv(os.path.join(folder_path, "blastoise_df.csv"), index=False, header=True, sep=",")  # Save the data frame to a CSV file
        print("")
        print(f"6. Stopping:")
        print(f"\t- No new data found.")
        print(f"\t- BLASTOISE final data row length: {coincidence_data.shape[0]}")

        return

    else:
        # -----------------------------------------------------------------------------
        toc_main = time.perf_counter()
        print("")
        print(f"RUN {numbering} finished:\n",
              f"\t- Execution time: {toc_main - tic_main:0.2f} seconds")
        # -----------------------------------------------------------------------------
        numbering += 1  # Increase the numbering
        repetitive_blaster(data_input=new_data_and_old,
                           genome_fasta=genome_fasta,
                           folder_path=folder_path,
                           numbering=numbering,
                           start_time=start_time,
                           identity_1=identity_1,
                           tic_start=tic_start,
                           word_size=word_size,
                           min_length=min_length,
                           extend_number=extend_number,
                           coincidence_data=coincidence_data)
                        