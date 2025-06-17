import pandas as pd
import subprocess
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fasta_creator(
        data_input: pd.DataFrame,
        fasta_output_path: str,
        id_names: Optional[pd.DataFrame] = None
) -> None:
    """
    Creates a FASTA file from a DataFrame containing sequence and metadata.

    This function processes sequence data from a specified DataFrame, formats the information as FASTA records, and 
    writes those records to a file in FASTA format. Optionally, it can include original identifiers provided via a 
    secondary DataFrame, appending both original and extended identifiers to each sequence record.

    Parameters
    ----------
    data_input : pd.DataFrame
        Input data containing sequence information. Assumes the input DataFrame includes a 'sseqid', 'sstart', 'send', 
        'sstrand', and 'sseq' column for the necessary metadata.
    fasta_output_path : str
        File path where the resulting FASTA file will be written.
    id_names : pd.DataFrame, optional
        DataFrame containing original sequence identifiers and metadata, used to further annotate FASTA records. 
        Default is None.
    
    Returns
    -------
    None
    """
    if id_names is None:
        id_names = data_input

    matrix = []
    for index, (_, sequence) in enumerate(data_input.iterrows()):
        # index += 1 # To start the index in 1
        data_input_string = f"{sequence['sseqid']}-{sequence['sstart']}-{sequence['send']}-{sequence['sstrand']}"

        if id_names is not None:
            # If there are the original coordinates and the extended ones, record them both
            original = f"{id_names.iloc[index]['sseqid']}-" \
                       f"{id_names.iloc[index]['sstart']}-" \
                       f"{id_names.iloc[index]['send']}-" \
                       f"{id_names.iloc[index]['sstrand']}"
            data_input_string = f"{data_input_string}_{original}" # The first is the extended, the second is the original

        rec = SeqRecord(Seq(sequence.loc["sseq"]),  # In the 5 position is the seq
                        id="Seq-" + str(index) + "_" + data_input_string,
                        description=""
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")


def get_data_sequence(data: pd.DataFrame, genome_fasta: str) -> pd.DataFrame:
    """
    Fetches nucleotide sequences from a specified genome database using BLAST commands.

    This function retrieves specified sequence ranges from a genome database in a fasta format.
    The input strand direction is specified, and BLAST commands are executed to obtain the sequences based on the 
    provided data.

    Parameters
    ----------
    data : pandas.DataFrame
        A pandas DataFrame containing the sequence details. It must
        contain the columns 'sseqid', 'sstart', and 'send'.
    genome_fasta : str
        The file path to the genome database in fasta format to query against.

    Returns
    -------
    pandas.DataFrame
        A pandas DataFrame containing the retrieved sequences. The DataFrame includes
        columns for 'sseqid', 'sstart', 'send', 'sstrand', and 'sseq'.
    """
    sequences = []
    for _, row in data.iterrows():
        sseqid = row['sseqid']
        start = row['sstart']
        end = row['send']
        strand = row['sstrand']
        cmd = f"blastdbcmd -db {genome_fasta} -entry {sseqid} -range {start}-{end} -strand {strand} -outfmt %s"

        sequence = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True,
                                  executable="/usr/bin/bash").stdout.strip()

        sequences.append({
            'sseqid': sseqid,
            'sstart': start,
            'send': end,
            'sstrand': strand,
            'sseq': sequence
        })

    sequences_df = pd.DataFrame(sequences)

    return sequences_df


def columns_to_numeric(data_input: pd.DataFrame, columns_to_convert: list[str] | None = None) -> pd.DataFrame:
    """
    Convert specified columns of a DataFrame to numeric datatype.

    This function takes a DataFrame and a list of column names and converts the specified columns to a numeric datatype 
    using `pd.to_numeric`. If no columns are specified, it defaults to converting the columns: 'length', 'sstart', and 
    'send'. Any non-convertible values are coerced into NaN.

    Parameters
    ----------
    data_input : pandas.DataFrame
        The pandas DataFrame containing the data to be transformed.
    columns_to_convert : list[str], optional
        A list of column names to be converted to numeric data type.
        Defaults to ['length', 'sstart', 'send'].

    Returns
    -------
    pandas.DataFrame
        The modified DataFrame with specified columns converted to numeric.
    """
    if columns_to_convert is None:
        columns_to_convert = ['sstart', 'send']

    for column in columns_to_convert:
        data_input[column] = pd.to_numeric(data_input[column], errors='coerce')
    return data_input


def end_always_greater_than_start(data_input: pd.DataFrame) -> pd.DataFrame:
    """
    Ensures that the 'send' value is always greater than 'sstart' value by swapping them if needed.
    
    Parameters:
    -----------
    data_input: pd.DataFrame
        A DataFrame containing 'sstart' and 'send' columns

    Returns:
    --------
    pd.DataFrame:
        A DataFrame with 'sstart' and 'send' values properly ordered
    """
    data_input.loc[
        data_input[data_input['sstart'] > data_input['send']].index,
        ['sstart', 'send']
    ] = data_input.loc[
        data_input[data_input['sstart'] > data_input['send']].index,
        ['send', 'sstart']
    ].values

    return data_input