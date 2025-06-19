import os
import pandas as pd
import subprocess
import logging
from typing import Optional, Tuple
import pyranges as pr

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add the parent directory of 'blastoise' to sys.path if needed
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from modules.blaster import blastn_dic

# ======================================================================
def csv_to_fasta_creator(csv_data: pd.DataFrame, fasta_output_path: str) -> None:
    """
    Convert CSV data into a FASTA format file.

    This function reads a CSV-like data structure, extracts specific sequence information,
    and writes it to a FASTA file format. It iterates through the rows of the provided
    dataframe, converts each sequence into a SeqRecord object, and appends it to a list.
    Finally, the list of SeqRecords is saved to the specified output file in FASTA format.

    Parameters:
    -----------
    csv_data : pandas.DataFrame
        A dataframe containing the sequence data. It should include at least two
        columns, 'sseq' for the sequence string and 'sseqid' for sequence identifier.
    fasta_output_path : str 
        The path to the output FASTA file where the converted sequences will be saved.
    """
    matrix = []
    for csv_index, sequence in csv_data.iterrows():
        rec = SeqRecord(Seq(sequence.loc[0, 'sseq']),
                        id=f"Seq_{csv_index}_{sequence.loc[0, 'sseqid']}",
                        description=""
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, 'fasta')

# ======================================================================
def fetch_dna_sequence(start_coor: int, end_coor: int, strand: str, chromosome: str, path_genome: str) -> str:
    """
    Retrieve a DNA sequence from a genome database using BLAST+ commands.

    This function fetches a genomic sequence from a specified chromosome
    and coordinates the range based on the strand orientation. The function
    builds a `blastdbcmd` command to interact with the genome database and
    extract the required sequence. The output is captured and returned as
    a string with whitespace stripped.

    Parameters:
    -----------
    start_coor: int
        The start coordinate of the desired sequence range.
    end_coor: int 
        The end coordinate of the desired sequence range.
    strand: str
        The strand orientation ('plus' or 'minus').
    chromosome: str
        The chromosome identifier in the genome database.
    path_genome: str
        The file path to the BLAST genome database.

    Returns:
    --------
        str: The retrieved genomic sequence as a string.

    Raises:
    -------
        None
    """
    cmd = f'blastdbcmd -db {path_genome} -entry {chromosome} -range {start_coor}-{end_coor} -strand {strand} -outfmt %s'
    sequence = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable="/usr/bin/bash")
    sequence = sequence.stdout.strip()
    return sequence

# ======================================================================================================================
def pyranges_merge(input_df: pd.DataFrame, path_folder: str = None) -> pd.DataFrame:
    """
    Process input data by merging overlapping genomic intervals using PyRanges.
    This function replaces the previous bedops_merge function, providing the same
    functionality but using PyRanges instead of BEDOPS tools.

    Parameters:
    -----------
    input_df: pd.DataFrame
        Input DataFrame containing at least 'qstart' and 'qend' columns, which represent genomic interval start and end 
        positions.
    path_folder: str, optional
        Path to a folder. This parameter is kept for backward compatibility but is not used.

    Returns:
    --------
        pd.DataFrame: A DataFrame containing the merged interval data with columns
            'sseqid', 'qstart', and 'qend'.
    """
    # Create a copy of the input dataframe with the required columns
    data_pyranges = input_df[["qstart", "qend"]].copy()

    # Add a dummy chromosome column (required by PyRanges)
    data_pyranges.insert(0, "Chromosome", "test")

    # Rename columns to match PyRanges expectations
    data_pyranges = data_pyranges.rename(columns={
        "qstart": "Start",
        "qend": "End"
    })

    # Convert to PyRanges
    pr_data = pr.PyRanges(data_pyranges)

    # Merge overlapping intervals
    merged = pr_data.merge()

    # Convert back to DataFrame and rename columns
    result = merged.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "qstart",
        "End": "qend"
    })

    return result

# ======================================================================
# BLAST FUNCTIONS
# ======================================================================
def general_blastn_blaster(query_path: str, dict_path: str, word_size: int, 
                          perc_identity: Optional[float] = None, evalue: Optional[float] = None) -> pd.DataFrame:
    """
    Executes a BLASTN command with the specified parameters and processes the output into a
    structured DataFrame.

    Parameters:
    -----------
    query_path: str
        The path to the file containing the query sequences.
    dict_path: str
        The path to the dictionary database for the BLASTN search.
    word_size: int 
        The word size for the BLASTN algorithm, specifying the size of
        initial matches during the search.
    perc_identity: Optional[float]
        The minimum percentage identity required for the matches
        (default is None and does not apply this filter if not provided).
    evalue: Optional[float]
        The e-value threshold for reporting matches. If not provided,
        all matches are returned.

    Returns:
    --------
        pd.DataFrame: A DataFrame containing the BLASTN search results. The DataFrame includes
            the following columns (all numeric columns are converted to numeric types):
            - "qseqid": Query sequence ID
            - "sseqid": Subject sequence ID 
            - "sstrand": Strand of the subject sequence
            - "qstart": Start position on the query sequence for the alignment
            - "qend": End position on the query sequence for the alignment
            - "sstart": Start position on the subject sequence for the alignment
            - "send": End position on the subject sequence for the alignment
            - "evalue": E-value of the alignment
            - "bitscore": a bit score of the alignment
            - "length": Length of the aligned region
            - "qlen": Length of the query sequence
            - "slen": Length of the subject sequence
            If no results are found, an empty DataFrame with these columns is returned.

    Raises:
    -------
        subprocess.SubprocessError: If the BLASTN subprocess encounters an issue during execution.
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

# ======================================================================
def setup_directories(base_dir: str, dict_path: str, temp_dir_name: str, logger_name: str) -> Tuple[str, str]:
    """
    Create the necessary directories and prepare a BLASTN database.

    Parameters
    ----------
    base_dir : str
        Base directory where temporary files will be created.
    dict_path : str
        Path to the genome FASTA file used to create the BLASTN database.
    temp_dir_name : str
        Name of the temporary directory to create.
    logger_name : str
        Name of the logger to use for logging messages.

    Returns
    -------
    Tuple[str, str]
        A tuple containing:
        - temp_dir (str): Path to the created temporary directory for intermediate files
        - blastn_db_path (str): Path to the generated BLASTN database created from the input genome FASTA
    """
    logger = logging.getLogger(logger_name)
    logger.info("Setting up directories and BLASTN database")

    # Ensure base directory exists
    os.makedirs(base_dir, exist_ok=True)

    # Prepare a subfolder for temporary files
    temp_dir = os.path.join(base_dir, temp_dir_name)
    os.makedirs(temp_dir, exist_ok=True)

    # Prepare BLASTn dict
    dict_folder_path = os.path.join(temp_dir, "blastn_dict")
    os.makedirs(dict_folder_path, exist_ok=True)

    blastn_db_path = os.path.join(dict_folder_path, os.path.basename(dict_path))

    try:
        blastn_dic(path_input=dict_path, path_output=blastn_db_path)
        logger.info(f"BLASTN database created at {blastn_db_path}")
    except Exception as e:
        logger.error(f"Error creating BLASTN database: {str(e)}")
        raise

    return temp_dir, blastn_db_path
