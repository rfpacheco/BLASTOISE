import subprocess
import pandas as pd
from typing import Dict, Any, Tuple
from joblib import Parallel, delayed


def _process_single_row_extension(
        row_data: Tuple[Any, pd.Series],
        genome_fasta: str,
        extend_number: int,
        limit_len: int
) -> Dict[str, Any]:
    """
    Process a single row for sequence extension.

    This helper function processes a single row from the DataFrame, evaluates the sequence length,
    and extends the sequence if it's shorter than the specified limit. The extension is done by
    adjusting the start and end coordinates, ensuring they don't go beyond valid boundaries,
    and then retrieving the extended sequence using a BLAST command.

    Parameters
    ----------
    row_data : Tuple[Any, pd.Series]
        A tuple containing the index and the row data from the DataFrame.
    genome_fasta : str
        The path to the genome FASTA file used for extracting sequence data.
    extend_number : int
        The minimum length required for the sequence. If a sequence is shorter than this
        value, it's extended to this length.
    limit_len : int
        A specific length limit that if met, the process of extension is skipped.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the processed row data with keys:
        'index', 'length', 'sstart', 'send', 'sseq', and a flag 'modified'.
    """
    index, element = row_data
    lower_coor = element['sstart']
    upper_coor = element['send']

    subject_len = upper_coor - lower_coor + 1  # Calculate the nucleotide size of the element

    # Default return values (no modification)
    result = {
        'index': index,
        'modified': False
    }

    if subject_len < limit_len:  # If the sequence is less than limit_len, there's room to expand it
        lower_coor = lower_coor - extend_number  # Extends the lower sequences by `extend_number`
        upper_coor = upper_coor + extend_number  # Extends the upper sequence by `extend_number`

        # Take into account that the numbers cannot be negative neither go against the limit length of the sequence
        if lower_coor <= 0:  # Because in BLASTn, 0 does not exist, it starts in 1
            lower_coor = 1

        # `upper_coor` does not need to change, because BLASTn won't extract longer values that the maximum sequence length.
        # but let's calculate, to take into account the sizes.
        cmd = f"blastdbcmd -db {genome_fasta} -entry {element['sseqid']} -outfmt '%l'"  # Command to extract the max len
        chrom_max_len = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()  # Extract the max len
        chrom_max_len = int(chrom_max_len)  # Convert to int
        if upper_coor > chrom_max_len:
            upper_coor = chrom_max_len

        new_subject_len = upper_coor - lower_coor + 1  # Calculate the nucleotide size of the element after extension.
        if new_subject_len < limit_len:  # If the sequence is still less than `limit_len` after the adjustment, continue the extension
            # Get the sequence, but extended
            cmd = (
                f"blastdbcmd -db {genome_fasta} "
                f"-entry {element['sseqid']} "
                f"-range {lower_coor}-{upper_coor} "
                f"-strand {element['sstrand']} "
                "-outfmt %s"
            )

        seq = subprocess.check_output(cmd, shell=True, universal_newlines=True).strip()

        # Return the modified values
        result.update({
            'modified': True,
            'length': new_subject_len,
            'sstart': int(lower_coor),
            'send': int(upper_coor),
            'sseq': seq
        })

    return result


def sequence_extension(
        data_input: pd.DataFrame,
        genome_fasta: str,
        extend_number: int,
        limit_len: int,
        n_jobs: int = -1
) -> pd.DataFrame:
    """
    Adjusts subject sequences in a data frame to meet a specified length requirement.

    Iterates through a data frame containing sequence alignment details, evaluates the length
    of each subject sequence, and adjusts the coordinates if the sequence is shorter than
    the specified limit. If the sequence is shorter, it extends the sequence by adjusting
    start and end positions equally on both sides while avoiding coordinates outside genome
    bounds. Also updates relevant columns in the data frame with the new sequence details.
    Uses parallel processing to improve performance when dealing with multiple chromosomes.

    Parameters
    ----------
    data_input : pd.DataFrame
        A data frame containing sequence alignment details. Expected to have the following 
        columns: 'sstrand', 'sstart', 'send', 'sseqid', and alignment metrics such as
        'pident', 'length', 'qstart', etc.
    genome_fasta : str
        The path to the genome FASTA file used for extracting sequence data.
    extend_number : int 
        The minimum length required for the sequence. If a sequence is shorter than this
        value, it's extended to this length.
    limit_len : int
        A specific length limit that if met, the process of extension is skipped.
    n_jobs : int, optional
        The number of jobs to run in parallel. -1 means using all processors.
        Default is -1.

    Returns
    -------
    pd.DataFrame
        A modified data frame where sequences having lengths below the specified minimum
        are extended. Updated metrics and extracted sequence are reflected in the modified
        data frame.

    Raises
    ------
    subprocess.CalledProcessError
        If an error occurs when running the system command to retrieve the genome sequence.
    KeyError
        If any required column is missing in the provided data frame.

    Notes:
    ------
    The function assumes that the input 'data_input' follows a specific structure with
    the necessary columns to process coordinate adjustments and sequence extraction.
    Coordinates are adjusted considering both strands ('plus' and '-') with careful attention
    to avoid exceeding genome boundaries.
    """
    # Process rows in parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(_process_single_row_extension)(row_data, genome_fasta, extend_number, limit_len)
        for row_data in data_input.iterrows()
    )

    # Update the DataFrame with the results
    for result in results:
        if result['modified']:
            index = result['index']
            data_input.loc[index, 'length'] = result['length']
            data_input.loc[index, 'sstart'] = result['sstart']
            data_input.loc[index, 'send'] = result['send']
            data_input.loc[index, 'sseq'] = result['sseq']

    return data_input
