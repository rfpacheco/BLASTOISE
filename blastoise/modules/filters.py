import time
import pandas as pd

from modules.genomic_ranges import get_merge_stranded


def global_filters_main(data_input: pd.DataFrame, min_length: int) -> pd.DataFrame:
    """
    Processes a given DataFrame by applying filtering and modification operations. It adjusts lengths,
    filters data by a specified minimum length, removes dashes from coordinates, and optionally processes
    with the BEDOPS method if data remains after filtering.

    Parameters
    ----------
    data_input : pd.DataFrame
        The input DataFrame containing columns 'send' and 'sstart'.
    min_length : int 
        The minimum acceptable length for filtering rows.

    Returns
    -------
    pd.DataFrame
        A DataFrame that has been processed and filtered according to the specified rules.
    """
    print("\t\t\t- Setting length:")
    data_input["length"] = abs(data_input["send"] - data_input["sstart"]) + 1

    print("\t\t\t- Filtering by length:")
    data_filtered = data_input[data_input["length"].astype(int) >= min_length]  # Filter by length

    print("\t\t\t- Removing dashes from coordinates:")
    data_filtered = data_filtered.apply(lambda x: x.replace("-", ""))  # Filter dashes

    final_data = data_filtered.copy()
    if not data_filtered.empty:  # Checks if the data is empty. If it is, it will skip the next part of the code
        tic = time.perf_counter()
        print("\t\t\t- Processing with BEDOPS:")
        data_bedops = get_merge_stranded(data_filtered)
        final_data = data_bedops.copy()
        toc = time.perf_counter()
        print(f"\t\t\t\t- Execution time: {toc - tic:0.2f} seconds")
    else:
        pass

    return final_data