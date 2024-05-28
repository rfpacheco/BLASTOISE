import pandas as pd
import os

from modules.bedops import bedops_coincidence

# Compare the coordinates
def coincidence_counter (df1, df2):
    coincidence = 0
    for _, row in df2.iterrows():
        if row["sstart"] in list(df1["sstart"].values) and row["send"] in list(df1["send"].values):
            coincidence += 1
    return coincidence

def stopping_main(data_df1, data_df2):
    data_df1_plus = data_df1[data_df1["sstrand"] == "plus"].copy()
    data_df1_minus = data_df1[data_df1["sstrand"] == "minus"].copy()
    data_df2_plus = data_df2[data_df2["sstrand"] == "plus"].copy()
    data_df2_minus = data_df2[data_df2["sstrand"] == "minus"].copy()

    coincidence_plus = coincidence_counter(data_df1_plus, data_df2_plus)
    coincidence_minus = coincidence_counter(data_df1_minus, data_df2_minus)

    total_coincidence = coincidence_plus + coincidence_minus
    perc_coincidence = total_coincidence / data_df2.shape[0] * 100
    print(f"\t\t\t- Coincidence with last corrected sequences:\n",
          f"\t\t\t\t- {total_coincidence}/{data_df2.shape[0]} - {perc_coincidence:.2f}%")

    if total_coincidence == data_df2.shape[0]:
        print(f"\t\t\t\t- TRUE")
        return True
    else:  # If the the coincidence is not the 100%
        print(f"\t\t\t\t- FALSE")
        return False
    
def stopping_bedops(data_df1, data_df2, folder_path, genome_fasta):
    data_df1_plus = data_df1[data_df1["sstrand"] == "plus"].copy()
    data_df1_minus = data_df1[data_df1["sstrand"] == "minus"].copy()
    data_df2_plus = data_df2[data_df2["sstrand"] == "plus"].copy()
    data_df2_minus = data_df2[data_df2["sstrand"] == "minus"].copy()

    # Prepare paths
    plus_path = os.path.join(folder_path, "bedops_coincidence_plus")
    minus_path = os.path.join(folder_path, "bedops_coincidence_minus")

    # -----------------------------------------------------------------------------
    ## Call BEDOPS on plus
    coincidence_plus, data_plus = bedops_coincidence(data_df1_plus, data_df2_plus, plus_path, "plus", genome_fasta)

    # Call BEDOPS on minus. Special case, because BEDOPS reads the coordinates like the "+" strand.
    ## First modify the coordinates.
    data_df1_minus[["sstart", "send"]] = data_df1_minus[["send", "sstart"]]
    data_df2_minus[["sstart", "send"]] = data_df2_minus[["send", "sstart"]]

    # -----------------------------------------------------------------------------
## And now call BEDOPS on minus
    coincidence_minus, data_minus = bedops_coincidence(data_df1_minus, data_df2_minus, minus_path. "minus", genome_fasta)
    if not data_minus.empty:  # If the data frame is not empty
        data_minus[["sstart", "send"]] = data_minus[["send", "sstart"]]  # restore "data_minus" coordinates

    recapture_data = pd.concat([data_plus, data_minus], ignore_index=True)  # joins both Data Frames
    if not recapture_data.empty:  # If the data frame is not empty
        recapture_data = recapture_data.sort_values(by=["sstrand", "sseqid", "sstart"])  # Sort the data frame by the start coordinate
 
    # -----------------------------------------------------------------------------
    ## math part
    total_coincidence = coincidence_plus + coincidence_minus
    perc_coincidence = total_coincidence / data_df2.shape[0] * 100
    print(f"\t\t\t- Coincidence (BEDOPS version):\n",
          f"\t\t\t\t- {total_coincidence}/{data_df2.shape[0]} - {perc_coincidence:.2f}%")
    
    # -----------------------------------------------------------------------------
    ## Now recapture data

    
    if total_coincidence == data_df2.shape[0]:
        print(f"\t\t\t\t- TRUE")
        return True, recapture_data
    else:  # If the the coincidence is not the 100%
        print(f"\t\t\t\t- FALSE")
        return False, recapture_data

