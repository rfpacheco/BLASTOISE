import pandas as pd
import pyranges as pr


# Column mapping constants
BLAST_TO_PYRANGES = {
    "sseqid": "Chromosome",
    "sstart": "Start",
    "send": "End"
}

PYRANGES_TO_BLAST = {
    "Chromosome": "sseqid",
    "Start": "sstart",
    "End": "send"
}


def to_pyranges_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Renames specific columns in a Pandas DataFrame to standardized names for PyRanges.

    Parameters:
    -----------
    df: pd.DataFrame
        A Pandas DataFrame containing the columns 'sseqid', 'sstart', and 'send'.

    Returns:
    --------
    pd.DataFrame
        A new DataFrame with renamed columns for PyRanges compatibility.
    """
    return df.rename(columns=BLAST_TO_PYRANGES)


def from_pyranges_format(df: pd.DataFrame) -> pd.DataFrame:
    """
    Renames PyRanges columns back to the original BLAST format.

    Parameters:
    -----------
    df: pd.DataFrame
        A Pandas DataFrame with PyRanges column names.

    Returns:
    --------
    pd.DataFrame
        A DataFrame with columns renamed back to BLAST format.
    """
    return df.rename(columns=PYRANGES_TO_BLAST)


def get_interval_overlap(df: pd.DataFrame, interval_df: pd.DataFrame, invert: bool = False) -> pd.DataFrame:
    """
    Finds intervals between two dataframes that either overlap or don't overlap using PyRanges.

    This function takes two dataframes as input, converts them into PyRanges
    objects, and calculates either the overlapping or non-overlapping intervals
    between them based on the 'invert' parameter.

    Parameters
    ----------
    df : pd.DataFrame
        The main dataframe containing genomic intervals for which overlap
        will be calculated.
    interval_df : pd.DataFrame
        The dataframe containing intervals to check for overlap with the main
        dataframe.
    invert : bool, default=False
        If True, returns intervals in df that do not overlap with interval_df.
        If False, returns intervals in df that do overlap with interval_df.

    Returns
    -------
    pd.DataFrame
        A dataframe containing the requested genomic intervals with standardized
        column names.
    """
    # Convert dataframes to PyRanges format
    pr_df = pr.PyRanges(to_pyranges_format(df))
    pr_intervals = pr.PyRanges(to_pyranges_format(interval_df))

    # Perform an overlap operation
    result_pr = pr_df.overlap(pr_intervals, strandedness=False, invert=invert)

    # Convert back to DataFrame with original column names
    if hasattr(result_pr, 'df'):
        return from_pyranges_format(result_pr.df)
    else:
        return pd.DataFrame(columns=['sseqid', 'sstart', 'send'])


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping or adjacent intervals within a DataFrame.

    This function takes a DataFrame containing genomic interval information, formats
    it to conform with PyRanges input requirements, merges overlapping or adjacent
    intervals, and returns the processed intervals in the original DataFrame structure.

    Parameters:
        df (pd.DataFrame): A DataFrame containing the input intervals.

    Returns:
        pd.DataFrame: A DataFrame with merged genomic intervals, formatted back to
        the original structure with renamed column labels.
    """
    # Convert to PyRanges, merge, and convert back
    pr_df = pr.PyRanges(to_pyranges_format(df))
    merged = pr_df.merge()

    # Handle an empty result
    if hasattr(merged, 'df'):
        return from_pyranges_format(merged.df)
    else:
        return pd.DataFrame(columns=['sseqid', 'sstart', 'send'])


def get_merge_stranded(data_input: pd.DataFrame) -> pd.DataFrame:
    """
    Processes genomic data to merge overlapping ranges based on strand information
    using PyRanges technology.

    This function replaces the previous bedops_main function, providing the same
    functionality but using PyRanges instead of BEDOPS tools.

    Arguments:
        data_input (pd.DataFrame): Input data containing genome coordinates and strand
            information. Expected to have specific columns like 'sstrand', 'sseqid',
            'sstart', 'send', etc.

    Returns:
        pd.DataFrame: A processed DataFrame containing genomic ranges with associated
            sequence data and additional calculated attributes such as sequence length.
    """
    result_dfs = []

    # Process each strand separately
    for strand in ['plus', 'minus']:
        # Filter and sort data for the current strand
        strand_df = data_input[data_input['sstrand'] == strand].copy()

        if not strand_df.empty:
            # Sort and merge intervals
            strand_df = strand_df.sort_values(by=['sseqid', 'sstart'])
            merged_df = merge_intervals(strand_df)
            merged_df['sstrand'] = strand
            result_dfs.append(merged_df)

    # Combine results from both strands
    if result_dfs:
        all_data = pd.concat(result_dfs, ignore_index=True)
        # Add a length column
        all_data['len'] = all_data['send'] - all_data['sstart'] + 1
        return all_data
    else:
        # Return empty DataFrame with expected columns
        return pd.DataFrame(columns=['sseqid', 'sstart', 'send', 'sstrand', 'len'])


def _print_dataset_stats(description: str, data: pd.DataFrame, total: int, label: str) -> None:
    """
    Helper function to print dataset statistics.

    Parameters:
    -----------
    description : str
        Description of the data being reported.
    data : pd.DataFrame
        The dataset to report statistics for.
    total : int
        The total number of items for percentage calculation.
    label : str
        Label for the data (e.g., "New data", "Previous data").
    """
    count = data.shape[0]
    if total > 0:
        percentage = count / total * 100
        print(f"\t\t\t\t- {label}: {count}/{total} - {percentage:.2f}%")
    else:
        print(f"\t\t\t\t- {label}: {count}/{total}")


def compare_genomic_datasets(
        main_data: pd.DataFrame,
        data_for_contrast: pd.DataFrame,
        strand: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compares two genomic datasets to identify overlapping and non-overlapping regions.

    This function processes and identifies overlapping and non-overlapping genomic data between
    two datasets using PyRanges tools, sorts and merges the results.

    Parameters:
    -----------
    main_data : pd.DataFrame
        The primary dataset containing genomic data to be analyzed. Expected to have
        columns 'sseqid', 'sstart', 'send', and 'sstrand'.
    data_for_contrast : pd.DataFrame
        The secondary dataset used for comparison against the main dataset. Expected
        to have the same structure as main_data.
    strand : str
        The genomic strand orientation ('plus' or 'minus').

    Returns:
    --------
    tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        A tuple containing:
        - coincidence_data (pd.DataFrame):
            Data overlapping in both `main_data` and `data_for_contrast` after merging.
        - new_data (pd.DataFrame):
            Data unique to `main_data` (not found in `data_for_contrast`).
        - only_in_contrast_data (pd.DataFrame):
            Data unique to `data_for_contrast` (not found in `main_data`).
    """
    # Sort the data frames by the start coordinate
    main_data = main_data.sort_values(by=['sseqid', 'sstart'])
    data_for_contrast = data_for_contrast.sort_values(by=['sseqid', 'sstart'])

    # Get length of the datasets for reporting
    main_data_len = main_data.shape[0]
    data_contrast_len = data_for_contrast.shape[0]

    # Find overlapping regions
    print("\n\t\t\t- Coincidence data:")
    main_exists_in_contrast = get_interval_overlap(main_data, data_for_contrast, invert=False)
    _print_dataset_stats("Coincidence", main_exists_in_contrast, main_data_len, 
                        "New data in Previous data")

    contrast_exists_in_main = get_interval_overlap(data_for_contrast, main_data, invert=False)
    _print_dataset_stats("Coincidence", contrast_exists_in_main, data_contrast_len, 
                        "Previous data in New data")

    # Merge overlapping regions
    overlapping_data = pd.concat([main_exists_in_contrast, contrast_exists_in_main], ignore_index=True)
    merged_data = merge_intervals(overlapping_data)
    print(f"\t\t\t\t- Merged data: {merged_data.shape[0]}")

    # Add strand information to coincidence data
    coincidence_data = merged_data.copy()
    coincidence_data['sstrand'] = strand

    # Find non-overlapping regions
    print("\n\t\t\t- NOT coincidence data:")
    main_not_in_contrast = get_interval_overlap(main_data, data_for_contrast, invert=True)
    _print_dataset_stats("Non-coincidence", main_not_in_contrast, main_data_len, 
                        "New data NOT in Previous data")

    contrast_not_in_main = get_interval_overlap(data_for_contrast, main_data, invert=True)
    _print_dataset_stats("Non-coincidence", contrast_not_in_main, data_contrast_len, 
                        "Previous data NOT in New data")

    # Prepare return values, handling empty dataframes
    new_data = main_not_in_contrast if not main_not_in_contrast.empty else pd.DataFrame()
    only_in_contrast_data = contrast_not_in_main if not contrast_not_in_main.empty else pd.DataFrame()

    return coincidence_data, new_data, only_in_contrast_data
