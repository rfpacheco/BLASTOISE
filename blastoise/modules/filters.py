"""
BLASTOISE Module: Data Filtering and Matching Operations
======================================================

This module provides functions for performing operations on genomic data stored
in pandas DataFrames. The primary focus is on matching and filtering rows based
on shared genomic coordinate information across DataFrames.

The module contains two main functions:
1. `match_data`: Finds matching rows between DataFrames based on genomic coordinates
2. `match_data_and_remove`: Removes rows from a DataFrame that match coordinates in another

These functions enable efficient handling of genomic data by comparing coordinates 
(sseqid, sstart, send, sstrand) between DataFrames to either identify matching 
regions or filter out unwanted entries based on coordinate matches.

Author: R. Pacheco
"""

import pandas as pd


def match_data(data_input: pd.DataFrame, to_discard: pd.DataFrame) -> pd.DataFrame:
    """
    Match data from two DataFrames based on specific genomic coordinate columns.

    This function performs an inner join operation between the input DataFrame
    `data_input` and a subset of columns from the DataFrame `to_discard`. The
    matching is conducted over the genomic coordinate columns 'sseqid', 'sstart',
    'send', and 'sstrand'. The resulting DataFrame contains rows where the values
    in these columns align in both DataFrames.

    Parameters
    ----------
    data_input : pd.DataFrame
        The primary DataFrame to be matched against. Must contain the columns
        'sseqid', 'sstart', 'send', and 'sstrand'.
    to_discard : pd.DataFrame
        The DataFrame from which specific columns are used to find matches in
        `data_input`. Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing rows from `data_input` that share the same values
        in the specified columns with the corresponding rows in `to_discard`.

    Raises
    ------
    KeyError
        If any of the required columns are missing from either DataFrame.
    """

    matches: pd.DataFrame = data_input.merge(
        to_discard[['sseqid', 'sstart', 'send', 'sstrand']],
        on=['sseqid', 'sstart', 'send', 'sstrand'],
        how='inner'
    )

    return matches


def match_data_and_remove(data_input: pd.DataFrame, to_discard: pd.DataFrame) -> pd.DataFrame:
    """
    Remove rows from a DataFrame based on matching genomic coordinates from another DataFrame.

    This function identifies rows from the given `data_input` DataFrame that matches rows in the
    `to_discard` DataFrame based on genomic coordinate columns and removes these matched rows.
    The matching is performed using the columns 'sseqid', 'sstart', 'send', and 'sstrand'
    as the indices for comparison.

    Parameters
    ----------
    data_input : pd.DataFrame
        The input DataFrame from which rows should be evaluated and potentially removed.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.
    to_discard : pd.DataFrame
        A DataFrame containing rows to be matched and removed from the `data_input` DataFrame.
        Must contain the columns 'sseqid', 'sstart', 'send', and 'sstrand'.

    Returns
    -------
    pd.DataFrame
        A copy of the `data_input` DataFrame with rows matching the criteria removed.

    Raises
    ------
    KeyError
        If any of the required columns are missing from either DataFrame.
    """

    matches: pd.DataFrame = match_data(data_input, to_discard)

    # Remove only the matching data
    data_input: pd.DataFrame = data_input.loc[
        ~data_input.set_index(['sseqid', 'sstart', 'send', 'sstrand']).index.isin(
            matches.set_index(['sseqid', 'sstart', 'send', 'sstrand']).index
        )
    ].copy()

    return data_input


