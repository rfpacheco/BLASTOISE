import pandas as pd

from modules.bedops import bedops_coincidence
from modules.files_manager import df_columns_restore, columns_to_numeric

def compare_main(newest_data, contrast_data, genome_fasta):
    """
    Processes and compares genomic dataframes for visualization and analysis based on
    BEDOPS tool outputs for both plus and minus DNA strands.

    Parameters:
    newest_data: pandas.DataFrame
        The latest dataframe containing genomic data to compare.
    contrast_data: pandas.DataFrame
        The older dataframe containing genomic data to compare against.
    genome_fasta: str
        Path to the genome FASTA file for reference in BEDOPS processing.

    Returns:
    tuple
        A tuple containing:
        - coincidence_data (pandas.DataFrame): Dataframe with regions common between newest_data and contrast_data. It will contain the data already present in contrast_data. However, this data may be modified by new coordinates with newest_data
        - new_data (pandas.DataFrame): Dataframe with regions exclusive to the latest dataframe.
        - old_data_exclusive (pandas.DataFrame): Dataframe with regions exclusive to the older dataframe.
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
    coincidence_plus, new_data_plus, old_data_exclusive_plus= bedops_coincidence(newest_data_plus, contrast_data_plus, "plus", genome_fasta)
    if not new_data_plus.empty:  # If the data frame has lines
        new_data_plus = df_columns_restore(new_data_plus, newest_data)  # Restore the columns
        new_data_plus = columns_to_numeric(new_data_plus)  # Convert columns to numeric
    if not coincidence_plus.empty:  # If the data frame has lines
        coincidence_plus = df_columns_restore(coincidence_plus, newest_data)  # Restore the columns
        coincidence_plus = columns_to_numeric(coincidence_plus)  # Convert columns to numeric
    if not old_data_exclusive_plus.empty:  # If the data frame has lines
        old_data_exclusive_plus = df_columns_restore(old_data_exclusive_plus, contrast_data)  # Restore the columns
        old_data_exclusive_plus = columns_to_numeric(old_data_exclusive_plus)  # Convert columns to numeric
    # -----------------------------------------------------------------------------
    # Call BEDOPS on minus.
    # And now call BEDOPS on minus
    print("")
    print("\t\t- '-' strand:")
    coincidence_minus, new_data_minus, old_data_exclusive_minus = bedops_coincidence(newest_data_minus, contrast_data_minus, "minus", genome_fasta)
    # Restore the coordinates
    if not new_data_minus.empty:  # If the data frame is not empty
        new_data_minus = df_columns_restore(new_data_minus, newest_data)  # Restore the columns
        new_data_minus = columns_to_numeric(new_data_minus)
    if not coincidence_minus.empty:
        coincidence_minus = df_columns_restore(coincidence_minus, newest_data)
        coincidence_minus = columns_to_numeric(coincidence_minus)
    if not old_data_exclusive_minus.empty:
        old_data_exclusive_minus = df_columns_restore(old_data_exclusive_minus, contrast_data)
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

    
    


