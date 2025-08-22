import os
import pandas as pd
import subprocess
from typing import Optional, Tuple

# Add the parent directory of 'blastoise' to sys.path if needed
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from blastoise.modules.blaster import create_blast_database


# ======================================================================
# BLAST FUNCTIONS
# ======================================================================
def general_blastn_blaster(
        query_path: str,
        dict_path: str,
        word_size: int,
        perc_identity: Optional[int] = None,
        evalue: Optional[float] = None
) -> pd.DataFrame:
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
    cmd += " -outfmt '10 qseqid sseqid sstrand qstart qend sstart send pident evalue bitscore length qlen slen'"

    data = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True,
                          executable='/usr/bin/bash')
    data = data.stdout

    # Parse the output into a DataFrame
    data_df = pd.DataFrame([x.split(",") for x in data.split("\n") if x])

    if not data_df.empty:
        data_df.columns = ["qseqid", "sseqid", "sstrand", "qstart", "qend", "sstart", "send",
                           "pident", "evalue", "bitscore", "length", "qlen", "slen"]

        # Convert numeric columns to numeric type
        numeric_cols = ["qstart", "qend", "sstart", "send", "pident", "evalue", "bitscore",
                        "length", "qlen", "slen"]
        data_df[numeric_cols] = data_df[numeric_cols].apply(pd.to_numeric)
    else:
        # Create an empty DataFrame with the correct columns
        data_df = pd.DataFrame(columns=["qseqid", "sseqid", "sstrand", "qstart", "qend", "sstart", "send",
                                        "pident", "evalue", "bitscore", "length", "qlen", "slen"]
        )

    return data_df

# ======================================================================
def setup_directories(
        base_dir: str,
        dict_path: str,
        temp_dir_name: str,
) -> Tuple[str, str]:
    """
    Sets up the required directories for BLASTn database preparation and returns
    the paths to the temporary directory and the BLASTn database.

    The function ensures the base directory and temporary directory exist, creates
    a subdirectory for the BLASTn dictionary, and sets up a BLASTn database file
    based on the given input dictionary path.

    Parameters
    ----------
    base_dir : str
        Path to the base directory where directories will be set up.
    dict_path : str
        Path to the input dictionary file that will be used to create the BLASTn
        database.
    temp_dir_name : str
        Name of the subdirectory under base_dir for temporary files.

    Returns
    -------
    Tuple[str, str]
        A tuple containing the path to the temporary directory and the path to the
        created BLASTn database.

    Raises
    ------
    Exception
        Raised if there is an error creating the BLASTn database.
    """
    # Ensure base directory exists
    os.makedirs(base_dir, exist_ok=True)

    # Prepare a subfolder for temporary files
    temp_dir: str = os.path.join(base_dir, temp_dir_name)
    os.makedirs(temp_dir, exist_ok=True)

    # Prepare BLASTn dict
    dict_folder_path: str = os.path.join(temp_dir, "blastn_dict")
    os.makedirs(dict_folder_path, exist_ok=True)

    blastn_db_path: str = os.path.join(dict_folder_path, os.path.basename(dict_path))

    try:
        create_blast_database(path_input=dict_path, path_output=blastn_db_path)
        print(f"BLASTN database created at {blastn_db_path}")
    except Exception as e:
        print(f"Error creating BLASTN database: {str(e)}")
        raise

    return temp_dir, blastn_db_path
