import pandas as pd
import os

from modules.bedops import bedops_coincidence
from modules.files_manager import df_columns_restore, columns_to_numeric

def compare_main(last_df, old_df, folder_path, genome_fasta):
    # -----------------------------------------------------------------------------
    # Select data
    last_df_plus = last_df[last_df["sstrand"] == "plus"].copy()
    last_df_minus = last_df[last_df["sstrand"] == "minus"].copy()
    old_df_plus = old_df[old_df["sstrand"] == "plus"].copy()
    old_df_minus = old_df[old_df["sstrand"] == "minus"].copy()
    # -----------------------------------------------------------------------------
    # Prepare paths
    plus_path = os.path.join(folder_path, "bedops_coincidence_plus")
    minus_path = os.path.join(folder_path, "bedops_coincidence_minus")
    # -----------------------------------------------------------------------------
    # Call BEDOPS on plus
    print("")
    print("\t\t- '+' strand:")
    coincidence_plus, new_data_plus, old_data_exclusive_plus= bedops_coincidence(last_df_plus, old_df_plus, plus_path, "plus", genome_fasta)
    if not new_data_plus.empty:  # If the data frame is not empty
        new_data_plus = df_columns_restore(new_data_plus, last_df)  # Restore the columns
        new_data_plus = columns_to_numeric(new_data_plus)  # Convert columns to numeric
    if not coincidence_plus.empty:  # If the data frame is not empty
        coincidence_plus = df_columns_restore(coincidence_plus, last_df)  # Restore the columns
        coincidence_plus = columns_to_numeric(coincidence_plus)  # Convert columns to numeric
    if not old_data_exclusive_plus.empty:  # If the data frame is not empty
        old_data_exclusive_plus = df_columns_restore(old_data_exclusive_plus, old_df)  # Restore the columns
        old_data_exclusive_plus = columns_to_numeric(old_data_exclusive_plus)  # Convert columns to numeric
    # -----------------------------------------------------------------------------
    # Call BEDOPS on minus. Special case, because BEDOPS reads the coordinates like the "+" strand.
    # First modify the coordinates.
    last_df_minus[["sstart", "send"]] = last_df_minus[["send", "sstart"]]
    old_df_minus[["sstart", "send"]] = old_df_minus[["send", "sstart"]]

    # And now call BEDOPS on minus
    print("")
    print("\t\t- '-' strand:")
    coincidence_minus, new_data_minus, old_data_exclusive_minus = bedops_coincidence(last_df_minus, old_df_minus, minus_path, "minus", genome_fasta)
    # Restore the coordinates
    if not new_data_minus.empty:  # If the data frame is not empty
        new_data_minus[["sstart", "send"]] = new_data_minus[["send", "sstart"]]  # restore "new_data_minus" coordinates
        new_data_minus = df_columns_restore(new_data_minus, last_df)  # Restore the columns
        new_data_minus = columns_to_numeric(new_data_minus)
    if not coincidence_minus.empty:
        coincidence_minus[["sstart", "send"]] = coincidence_minus[["send", "sstart"]]
        coincidence_minus = df_columns_restore(coincidence_minus, last_df)
        coincidence_minus = columns_to_numeric(coincidence_minus)
    if not old_data_exclusive_minus.empty:
        old_data_exclusive_minus[["sstart", "send"]] = old_data_exclusive_minus[["send", "sstart"]]
        old_data_exclusive_minus = df_columns_restore(old_data_exclusive_minus, old_df)
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

    
    


