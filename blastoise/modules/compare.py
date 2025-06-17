import pandas as pd

from modules.bedops import bedops_coincidence
from modules.files_manager import columns_to_numeric


def compare_main(  # TODO: replace BEDOPS with PyRanges
        newest_data: pd.DataFrame,
        contrast_data: pd.DataFrame,
        genome_fasta: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    # TODO: replace BEDOS in docs to PyRanges
    This function performs comparative analysis between two data sets representing genomic
    features on different strands ('plus' and 'minus') of DNA. It invokes the BEDOPS tool
    for overlap analysis and processes the results to identify shared, unique, and exclusive
    data points. The concatenated results for both strands are returned separately.

    Parameters
    ----------
    newest_data : pd.DataFrame
        A DataFrame of the newest data to analyze, containing genomic feature information 
        such as strand and coordinates.
    contrast_data : pd.DataFrame
        A DataFrame of contrast data to compare against the newest data, with similar 
        genomic feature information.
    genome_fasta : str
        The genome FASTA file used for reference during the analysis.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        A tuple containing three DataFrames:
        - coincidence_data: Genomic features common between `newest_data` and
          `contrast_data`.
        - new_data: Genomic features exclusive to the `newest_data`.
        - old_data_exclusive: Genomic features exclusive to the `contrast_data`.
    """
    # -----------------------------------------------------------------------------
    # Select data
    newest_data_plus = newest_data[newest_data["sstrand"] == "plus"].copy()
    newest_data_minus = newest_data[newest_data["sstrand"] == "minus"].copy()
    contrast_data_plus = contrast_data[contrast_data["sstrand"] == "plus"].copy()
    contrast_data_minus = contrast_data[contrast_data["sstrand"] == "minus"].copy()
        # -----------------------------------------------------------------------------

    # Call BEDOPS on plus
    print("")
    print("\t\t- '+' strand:")
    coincidence_plus, new_data_plus, old_data_exclusive_plus= bedops_coincidence(
        newest_data_plus, contrast_data_plus, "plus", genome_fasta
    )
    if not new_data_plus.empty:  # If the data frame has lines
        new_data_plus = columns_to_numeric(new_data_plus)  # Convert columns to numeric
    if not coincidence_plus.empty:  # If the data frame has lines
        coincidence_plus = columns_to_numeric(coincidence_plus)  # Convert columns to numeric
    if not old_data_exclusive_plus.empty:  # If the data frame has lines
        old_data_exclusive_plus = columns_to_numeric(old_data_exclusive_plus)  # Convert columns to numeric
    # -----------------------------------------------------------------------------
    # Call BEDOPS on minus.
    # And now call BEDOPS on minus
    print("")
    print("\t\t- '-' strand:")
    coincidence_minus, new_data_minus, old_data_exclusive_minus = bedops_coincidence(
        newest_data_minus, contrast_data_minus, "minus", genome_fasta
    )
    # Restore the coordinates
    if not new_data_minus.empty:  # If the data frame is not empty
        new_data_minus = columns_to_numeric(new_data_minus)
    if not coincidence_minus.empty:
        coincidence_minus = columns_to_numeric(coincidence_minus)
    if not old_data_exclusive_minus.empty:
        old_data_exclusive_minus = columns_to_numeric(old_data_exclusive_minus)
    # -----------------------------------------------------------------------------
    # Concatenate both Data Frames
    if not new_data_plus.empty or not new_data_minus.empty:
        new_data = pd.concat([new_data_plus, new_data_minus], ignore_index=True)
    else:  # If both Data Frames are empty
        new_data = pd.DataFrame()

    if not coincidence_plus.empty or not coincidence_minus.empty:
        coincidence_data = pd.concat([coincidence_plus, coincidence_minus], ignore_index=True)  # joins both Data Frames
    else:  # If both Data Frames are empty
        coincidence_data = pd.DataFrame()

    if not old_data_exclusive_plus.empty or not old_data_exclusive_minus.empty:
        old_data_exclusive = pd.concat([old_data_exclusive_plus, old_data_exclusive_minus], ignore_index=True)
    else:  # If both Data Frames are empty
        old_data_exclusive = pd.DataFrame()
    # -----------------------------------------------------------------------------
    return coincidence_data, new_data, old_data_exclusive