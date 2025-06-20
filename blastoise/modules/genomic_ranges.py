import pandas as pd
import pyranges as pr


def pyranges_column_name_change(df: pd.DataFrame) -> pd.DataFrame:
    """
    Renames specific columns in a Pandas DataFrame to standardized names related
    to genomic data representation. The function is expected to handle a DataFrame
    containing columns such as 'sseqid', 'sstart', and 'send', and convert these
    to 'Chromosome', 'Start', and 'End', respectively.

    Parameters:
    -----------
    df: pd.DataFrame
        A Pandas DataFrame containing the columns 'sseqid', 'sstart', and 'send'.
        These columns are expected to represent genomic data that will be renamed
        to standardized identifiers for further processing.

    Returns:
    --------
    pd.DataFrame
        A new DataFrame with renamed columns:
        - 'sseqid' renamed to 'Chromosome'
        - 'sstart' renamed to 'Start'
        - 'send' renamed to 'End'

    """
    df_renamed = df.rename(columns={
        "sseqid": "Chromosome",
        "sstart": "Start",
        "send": "End"
    })

    return df_renamed


def get_interval_coincidence(df: pd.DataFrame, interval_df: pd.DataFrame) -> pd.DataFrame:
    """
    Finds overlap intervals between two dataframes using PyRanges.

    This function takes two dataframes as input, converts them into PyRanges
    objects, calculates the overlapping intervals between them, and then returns
    the result as a dataframe. Input dataframes should contain genomic interval
    information. The column names of the dataframes will be preprocessed to match
    the requirements of PyRanges, and the result will have its columns renamed
    back before returning.

    Parameters
    ----------
    df : pd.DataFrame
        The main dataframe containing genomic intervals for which overlap
        will be calculated.
    interval_df : pd.DataFrame
        The dataframe containing intervals to check for overlap with the main
        dataframe.

    Returns
    -------
    pd.DataFrame
        A dataframe containing overlapping genomic intervals with standardized
        column names. The columns of the resulting dataframe will include
        'sseqid', 'sstart', and 'send'.
    """

    # Rename columns to match PyRanges expectations
    df_renamed = pyranges_column_name_change(df)
    interval_df_renamed = pyranges_column_name_change(interval_df)

    # Convert to PyRanges
    pr_df = pr.PyRanges(df_renamed)
    pr_intervals = pr.PyRanges(interval_df_renamed)

    # Perform overlap (intersection)
    overlaps = pr_df.overlap(pr_intervals, strandedness=False) # TODO: strandedness could be a parameter?

    # Return as DataFrame and rename back
    result = overlaps.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result


def get_interval_not_coincidence(df: pd.DataFrame, interval_df: pd.DataFrame) -> pd.DataFrame:
    """
    Gets intervals from the first DataFrame that do not coincide with intervals in
    the second DataFrame. This function converts input data into PyRanges objects to
    efficiently perform interval operations and identifies intervals from the first
    DataFrame that do not overlap with those in the second DataFrame.

    Parameters:
    ----------
        df: DataFrame
            The first DataFrame containing intervals to compare. Expected columns
            are compatible with PyRanges, such as Chromosome, Start, and End.
        interval_df: DataFrame
            The second DataFrame representing intervals for comparison. Similar
            structure to the first DataFrame.

    Returns:
    ----------
        DataFrame
            A DataFrame of intervals from the first DataFrame that do not overlap
            with any intervals in the second DataFrame. The resulting DataFrame has
            columns renamed back to match the original format, such as sseqid,
            sstart, and send.
    """
    # Rename columns to match PyRanges expectations
    df_renamed = pyranges_column_name_change(df)
    interval_df_renamed = pyranges_column_name_change(interval_df)

    # Convert to PyRanges
    pr_df = pr.PyRanges(df_renamed)
    pr_intervals = pr.PyRanges(interval_df_renamed)

    # Perform overlap (intersection)
    non_overlapping = pr_df.overlap(pr_intervals, strandedness=False, invert=True) # TODO: strandedness could be a parameter?

    # Return as DataFrame and rename back
    result = non_overlapping.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result

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

    # Rename to PyRanges expected format
    pr_df = pyranges_column_name_change(df)

    # Get `pr_df` as a PyRanges obejct
    pr_df = pr.PyRanges(pr_df)

    # Merge (collapse overlapping/adjacent intervals)
    merged = pr_df.merge()

    # Convert back to DataFrame and rename
    result = merged.df.rename(columns={
        "Chromosome": "sseqid",
        "Start": "sstart",
        "End": "send"
    })

    return result


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
    # -----------------------------------------------------------------------------
    # 1) Filter and sort data
    # -----------------------------------------------------------------------------
    df_plus = data_input[data_input['sstrand'] == 'plus'].copy()  # filters the "+" strand
    df_minus = data_input[data_input['sstrand'] == 'minus'].copy()  # filters the "-" strand

    # Sort the data by the start coordinate
    df_plus = df_plus.sort_values(by=['sseqid', 'sstart'])  # sorts the "+" strand by the start coordinate
    df_minus = df_minus.sort_values(by=['sseqid', 'sstart'])  # sorts the "-" strand by the start coordinate

    # -----------------------------------------------------------------------------
    # 2) Merge overlapping intervals using PyRanges for each strand
    # -----------------------------------------------------------------------------
    if not df_plus.empty:
        df_plus_merged = merge_intervals(df_plus)
        df_plus_merged['sstrand'] = 'plus'
    else:
        df_plus_merged = pd.DataFrame(columns=['sseqid', 'sstart', 'send', 'sstrand'])

    if not df_minus.empty:
        df_minus_merged = merge_intervals(df_minus)
        df_minus_merged['sstrand'] = 'minus'
    else:
        df_minus_merged = pd.DataFrame(columns=['sseqid', 'sstart', 'send', 'sstrand'])

    # -----------------------------------------------------------------------------
    # 3) Combine results and add length column
    # -----------------------------------------------------------------------------
    # Join both data frames
    all_data = pd.concat([df_plus_merged, df_minus_merged], ignore_index=True)

    # Add "len" column
    all_data['len'] = all_data['send'] - all_data['sstart'] + 1

    return all_data


def compare_genomic_datasets(
        main_data: pd.DataFrame,
        data_for_contrast: pd.DataFrame,
        strand: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compares two genomic datasets to identify overlapping and non-overlapping regions.

    This function processes and identifies overlapping and non-overlapping genomic data between
    two datasets using PyRanges tools, sorts and merges the results. It replaces the previous
    bedops_coincidence function, providing the same functionality but using PyRanges instead
    of BEDOPS tools.

    Parameters:
    -----------
    main_data : pd.DataFrame
        The primary dataset containing genomic data to be analyzed. Expected to have
        columns 'sseqid', 'sstart', 'send', and 'sstrand'.
    data_for_contrast : pd.DataFrame
        The secondary dataset used for comparison against the main dataset. Expected
        to have the same structure as main_data.
    strand : str
        The genomic strand orientation ('+' or '-').

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

    # Check elements in the "main data" that overlap in the "contrast" data
    main_exists_in_contrast_data = get_interval_coincidence(main_data, data_for_contrast)
    print("")
    print("\t\t\t- Coincidence data:")
    if main_data_len > 0:
        print(f"\t\t\t\t- New data in Previous data: {main_exists_in_contrast_data.shape[0]}/{main_data_len} - {main_exists_in_contrast_data.shape[0]/main_data_len*100:.2f}%")
    else:  # In this case `main_data_len == 0`
        print(f"\t\t\t\t- New data in Previous data: {main_exists_in_contrast_data.shape[0]}/{main_data_len}")

    contrast_exists_in_main_data = get_interval_coincidence(data_for_contrast, main_data)
    if data_contrast_len > 0:
        print(f"\t\t\t\t- Previous data in New data: {contrast_exists_in_main_data.shape[0]}/{data_contrast_len} - {contrast_exists_in_main_data.shape[0]/data_contrast_len*100:.2f}%")
    else:  # Then `data_contrast_len == 0`
        print(f"\t\t\t\t- Previous data in New data: {contrast_exists_in_main_data.shape[0]}/{data_contrast_len}")

    # There would be elements that are in both datasets. The next step is to merge them.
    # Combine the overlapping data from both directions and merge
    overlapping_data = pd.concat([main_exists_in_contrast_data, contrast_exists_in_main_data], ignore_index=True)
    merged_data = merge_intervals(overlapping_data)
    print(f"\t\t\t\t- Merged data: {merged_data.shape[0]}")

    coincidence_data = merged_data.copy()
    coincidence_data['sstrand'] = strand

    # Check elements from `main_data` that DO NOT overlap with `data_for_contrast`
    # Because these elements are not in the contrast data, they will be novice elements.
    print("")
    print("\t\t\t- NOT coincidence data:")
    main_not_exists_in_contrast_data = get_interval_not_coincidence(main_data, data_for_contrast)
    if main_data_len > 0:
        print(f"\t\t\t\t- New data NOT in Previous data: {main_not_exists_in_contrast_data.shape[0]}/{main_data_len} - {main_not_exists_in_contrast_data.shape[0]/main_data_len*100:.2f}%")
    else:
        print(f"\t\t\t\t- New data NOT in Previous data: {main_not_exists_in_contrast_data.shape[0]}/{main_data_len}")

    if not main_not_exists_in_contrast_data.empty:  # If the data frame has data
        new_data = main_not_exists_in_contrast_data.copy()
    else:  # If the data frame is empty
        new_data = pd.DataFrame()

    # Now check the elements in Old that are not in Last
    contrast_not_exists_in_main_data = get_interval_not_coincidence(data_for_contrast, main_data)
    if data_contrast_len > 0:
        print(f"\t\t\t\t- Previous data NOT in New data: {contrast_not_exists_in_main_data.shape[0]}/{data_contrast_len} - {contrast_not_exists_in_main_data.shape[0]/data_contrast_len*100:.2f}%")
    else:
        print(f"\t\t\t\t- Previous data NOT in New data: {contrast_not_exists_in_main_data.shape[0]}/{data_contrast_len}")

    if not contrast_not_exists_in_main_data.empty:  # If the data frame has lines
        only_in_contrast_data = contrast_not_exists_in_main_data.copy()
    else:  # If the data frame is empty
        only_in_contrast_data = pd.DataFrame()

    return coincidence_data, new_data, only_in_contrast_data
