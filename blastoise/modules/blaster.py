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
from modules.compare import compare_main
from modules.files_manager import end_always_greater_than_start
from extra.csv_to_gff import csv_to_gff

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def blastn_dic(path_input, path_output):
    """
    Executes the BLAST `makeblastdb` command to create a nucleotide database using the provided input
    file and outputs the result to the specified output path.

    This function is a wrapper around the BLAST `makeblastdb` tool and is designed to create a
    BLAST dictionary based on the input nucleotide file. The sequence IDs are preserved using
    the `parse_seqids` flag. If an error occurs during the execution of the command, it will
    be logged.

    Arguments:
        path_input (str): The path to the input nucleotide file.
        path_output (str): The path where the BLAST database files will be saved.

    Raises:
        Exception: If an error occurs during the creation of the BLAST dictionary, it logs the
        error message.
    """
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
            qstart: Start of alignment in a query.
            qend: End of alignment in a query.
            sstart: Start of alignment in a subject.
            send: End of alignment in a subject.
            evalue: Expectation value.
            bitscore: Bit score.
            qlen: Length of a query sequence.
            slen: Length of a subject sequence.
            sstrand: Strand of the subject sequence.
            sseq: Aligned part of the subject sequence.
    """

    cmd = (f"blastn "
           f"-word_size {word_size} "
           f"-query {query_path} "
           f"-db {dict_path} "
           f"-perc_identity {perc_identity} "
           f"-outfmt '10 qseqid sseqid sstart send sstrand evalue sseq'")
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])
    data.columns = ['qseqid', 'sseqid', 'sstart', 'send', 'sstrand', 'evalue', 'sseq']
    # Get 'sstart', 'send' to 'int' type
    data[['sstart', 'send']] = data[['sstart', 'send']].astype(int)
    # get 'evalue' as a 'float' type
    data['evalue'] = data['evalue'].astype(float)
    # Make sure 'send' > 'sstart'
    data = end_always_greater_than_start(data)
    # Create 'len' colum
    data['len'] = data['send'] - data['sstart'] + 1
    # Place it between 'sent' and 'sstrand' column
    data = data[['qseqid', 'sseqid', 'sstart', 'send', 'sstrand', 'evalue', 'sseq', 'len']]
    return data
    

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def repetitive_blaster(data_input, genome_fasta, folder_path, numbering, start_time, identity_1, tic_start, word_size, min_length, extend_number, limit_len, coincidence_data=None):
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
    boxymcboxface('RUN ' + str(numbering))
    tic_main = time.perf_counter()  # Start the timer

    # -----------------------------------------------------------------------------
    tic = time.perf_counter()
    # First let's order the data by "sseqid", "sstrand", "sstart".
    data_ordered = data_input.sort_values(by=['sseqid', 'sstrand', 'sstart'])
    toc = time.perf_counter()
    print("")
    print('1. Initial data:\n',
          f"\t- Data row length: {data_input.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")
    terminal_width = shutil.get_terminal_size().columns  # Get the terminal width

    print("")
    print('2. Individual searching and cleaning:')
    tic = time.perf_counter()
    now_time = datetime.now()
    formatted_now_time = now_time.strftime('%Y %B %d at %H:%M')
    print("")
    print(f"{' ' * 7}{'-' * 74}")
    start_time_text = f"Program started: {start_time}"
    end_time_text = f"Program time now: {formatted_now_time}"
    run_text = f"RUN {numbering}"
    print(f"{run_text:>{terminal_width}}")
    print(f"{start_time_text:>{terminal_width}}")
    print(f"{end_time_text:>{terminal_width}}")

    whole_group = genome_specific_chromosome_main(data_input=data_ordered,
                                                  main_folder_path=folder_path,
                                                  genome_fasta=genome_fasta,
                                                  identity_1=identity_1,
                                                  run_phase=numbering,
                                                  coincidence_data=coincidence_data,
                                                  word_size=word_size,
                                                  min_length=min_length,
                                                  limit_len=limit_len,
                                                  extend_number=extend_number)

    toc = time.perf_counter()
    print(f"\t- Data row length: {whole_group.shape[0]}\n",
          f"\t- Execution time: {toc - tic:0.2f} seconds")

      # -----------------------------------------------------------------------------
    # Compare part
    # Prepare folders and path
    comparison_folder = os.path.join(folder_path, "comparison")
    os.makedirs(comparison_folder, exist_ok=True)

    print("")
    print(f"3. Comparison VS Previous Run:")
    if coincidence_data is not None:
        print("")
        print(f"\t- Previous Run data:\n",
              f"\t\t- Coincidence data row length: {coincidence_data.shape[0]}\n",
              f"\t\t- New data row length: {data_input.shape[0]}")
        # This part is important to compare with the last run whole data "whoel_group_fildered", and not only the "new_data" subset.
        data_input = pd.concat([coincidence_data, data_input], ignore_index=True)
        data_input.sort_values(by=["sseqid", "sstrand", "sstart"], inplace=True)  # Sort the data frame by the start coordinate
        print(f"\t\t- Total data row length: {data_input.shape[0]}")
    else:  # when coincidence_data == None
        print(f"\t- Previous Run data:\n",
              f"\t\t- First Run row length: {data_input.shape[0]}")

    tic = time.perf_counter()
    print("")
    print(f"\t- Results in this RUN:")
    coincidence_data, new_data, old_data_exclusive = compare_main(whole_group,
                                                                  data_input,
                                                                  genome_fasta)
    toc = time.perf_counter()


    print("")
    print(f"\t\t- Coincidence data from run 'n' and 'n-1': {coincidence_data.shape[0]}\n",
          f"\t\t- New data detected only from run 'n': {new_data.shape[0]}\n",
          f"\t\t- Previous data detecte only from 'n-1': {old_data_exclusive.shape[0]}")

    old_data_exclusive_less_than_100 = None

    old_data_exclusive['len'] = abs(old_data_exclusive['send'] - old_data_exclusive['sstart']) + 1
    # noinspection PyUnresolvedReferences
    # If `old_data_exclusive` has lines and there's some 'len' < `lim_length`
    if not old_data_exclusive.empty and (old_data_exclusive['len'] < min_length).sum() > 0:  # If there are sequences less than 100 bp. The sum of TRUE (for < 100) has to be > 0
        old_data_exclusive_less_than_100 = old_data_exclusive[old_data_exclusive['len'] < min_length]
        old_data_exclusive = old_data_exclusive[old_data_exclusive['len'] >= min_length]
        # Now drop the 'len' column for both, to not break the following concat
        old_data_exclusive_less_than_100.drop(columns=['len'], inplace=True)
        old_data_exclusive.drop(columns=['len'], inplace=True)
        print(f"\t\t\t- Sequences < {min_length} bp: {old_data_exclusive_less_than_100.shape[0]}")
        print(f"\t\t\t- Sequences >= {min_length} bp: {old_data_exclusive.shape[0]}")
    else:
        pass

    if old_data_exclusive_less_than_100 is not None: # If `old_data_exclusive_less_than_100` exists
        new_data_and_old = pd.concat([new_data, old_data_exclusive_less_than_100], ignore_index=True)
        new_data_and_old.sort_values(by=["sseqid", "sstrand", "sstart"], inplace=True)
        print('\t' * 3 + f"- New data + less than {min_length}: {new_data_and_old.shape[0]}")
    else:
        new_data_and_old = new_data

    # if `coincidence_data` and `old_data_exclusive` has lines --> join them;
    # This way we have the `coincidence_data` containing joined data from 'n' and 'n-1' run and
    # the exclusive data from 'n-1' that does not appear in 'n'
    if not coincidence_data.empty and not old_data_exclusive.empty:
        coincidence_data = pd.concat([coincidence_data, old_data_exclusive], ignore_index=True)
        coincidence_data.sort_values(by=['sseqid', 'sstrand', 'sstart'], inplace=True)
        print(f"\t\t- Coincidence data + Previous data: {coincidence_data.shape[0]}")

    else:
        pass
    print(f"\t\t- Execution time: {toc - tic:0.2f} seconds")
    # -----------------------------------------------------------------------------
    # Stopping part
    stopping_folder = os.path.join(folder_path, "stopping")
    os.makedirs(stopping_folder, exist_ok=True)

    if new_data.shape[0] == 0:
        coincidence_data = coincidence_data[['sseqid', 'sstart', 'send', 'sstrand']].copy()  # #Take only the necessary columns:
        
        # Make it so 'sstart' is always < than 'send'
        coincidence_data = end_always_greater_than_start(coincidence_data)


        coincidence_data.to_csv(os.path.join(folder_path, "blastoise_df.csv"), index=False, header=True, sep=",")  # Save the data frame to a CSV file
        csv_to_gff(os.path.join(folder_path, "blastoise_df.csv"))
        print("")
        print(f"4. Stopping:")
        print(f"\t- No new data found.")
        print(f"\t- BLASTOISE final data row length: {coincidence_data.shape[0]}")

        return

    else:
        # -----------------------------------------------------------------------------
        toc_main = time.perf_counter()
        print("")
        print(f"RUN {numbering} finished:\n",
              f"\t- Execution time: {toc_main - tic_main:0.2f} seconds")

        # Now save the information for this data
        save_run_file = pd.concat([new_data_and_old, coincidence_data], ignore_index=True)
        save_run_file.sort_values(by=['sseqid', 'sstart'], inplace=True)
        save_run_file['len'] = save_run_file['send'] - save_run_file['sstart'] + 1
        runs_folder = os.path.join(folder_path, "RUNS")  # Creates the folder for the RUNS
        os.makedirs(runs_folder, exist_ok=True)  # Creates the folder for the RUNS
        run_saver_path = os.path.join(runs_folder, "run_" + str(numbering) + ".csv")  # Path to save the RUN
        save_run_file.to_csv(run_saver_path, sep=",", header=True, index=False)  # Saves the RUN
        csv_to_gff(run_saver_path)
        print(f"\t- Coincidence data + 'n' exclusive data + 'n-1' exclusive data: {save_run_file.shape[0]}")
        print(f"\t- Data file saved at {run_saver_path}")

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
                           limit_len=limit_len,
                           coincidence_data=coincidence_data)
                        