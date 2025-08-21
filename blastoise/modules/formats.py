"""
BLASTOISE Module: Data Formatting and Export Operations
======================================================

This module provides core functionality for the BLASTOISE pipeline, handling data formatting
and export operations. It serves as the output processor for genomic alignment results,
converting BLAST output into standardized formats and exporting data for downstream analysis.

Main functions:
1. format_output_dataframe: Formats DataFrame by mapping columns to standardized names
2. write_gff_from_formatted: Exports formatted data to General Feature Format (GFF) files

Author: R. Pacheco
"""

import pandas as pd
import logging
from typing import Dict, List

logger = logging.getLogger(__name__)


def format_output_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Formats the input DataFrame by mapping specific column names to a standardized set of column names
    and returns the formatted DataFrame. If the input DataFrame is empty or None, a new DataFrame
    with predefined column names is returned.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame to be formatted. Should contain specific column names that are required
        for mapping to standardized column names.

    Returns
    -------
    pandas.DataFrame
        A formatted DataFrame with columns renamed according to the specified mapping. If the input
        DataFrame is empty or None, it returns a new DataFrame with predefined column names but no data.

    """
    if df is None or df.empty:
        logger.info("format_output_dataframe: empty input, returning predefined columns")
        return pd.DataFrame(columns=['chromosome', 'start', 'end', 'strand', 'len', 'seq'])

    cols_map: Dict[str, str] = {
        'sseqid': 'chromosome',
        'sstart': 'start',
        'send': 'end',
        'sstrand': 'strand',
        'len': 'len',
        'sseq': 'seq',
    }
    # Select only available required columns and rename
    available: List[str] = ['sseqid', 'sstart', 'send', 'sstrand', 'len', 'sseq']
    formatted: pd.DataFrame = df[available].rename(columns=cols_map).copy()
    logger.info("format_output_dataframe: formatted %d records", len(formatted))

    return formatted


def write_gff_from_formatted(df_formatted: pd.DataFrame, gff_path: str) -> None:
    """
    Writes a GFF (General Feature Format) file from a formatted DataFrame.

    This function converts a given formatted pandas DataFrame containing genomic features 
    into a GFF output file. The function allows translation of strand data and assigns 
    unique identifiers to each feature. If the input DataFrame is empty or None, an 
    empty GFF file is created.

    Parameters
    ----------
    df_formatted : pd.DataFrame
        A formatted DataFrame containing genomic feature data, including columns such 
        as 'chromosome', 'start', 'end', and optionally 'strand'.
    gff_path : str
        The path where the output GFF file will be written.

    Returns
    -------
    None
    """
    # Handle empty
    if df_formatted is None or df_formatted.empty:
        open(gff_path, "w").close()
        return

    strand_series: pd.Series = df_formatted.get('strand')
    strand: pd.Series
    if strand_series is not None:
        strand = strand_series.map({'plus': '+', 'minus': '-'}).fillna('.')
    else:
        strand = pd.Series(['.'] * len(df_formatted))

    ids: List[str] = [f'BLASTOISE_{i + 1}' for i in range(len(df_formatted))]

    gff_df: pd.DataFrame = pd.DataFrame({
        'seqname': df_formatted['chromosome'],  # type: str
        'source': 'BLASTOISE',  # type: str
        'feature': '.',  # type: str
        'start': df_formatted['start'],  # type: int
        'end': df_formatted['end'],  # type: int
        'score': '.',  # type: str
        'strand': strand,  # type: str
        'frame': '.',  # type: str
        'attribute': [f'ID={x}' for x in ids],  # type: str
    })

    gff_df.to_csv(gff_path, sep='\t', index=False, header=False)
    logger.info("write_gff_from_formatted: wrote %d features to %s", len(gff_df), gff_path)