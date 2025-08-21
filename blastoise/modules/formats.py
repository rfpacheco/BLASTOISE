


import pandas as pd

def format_output_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Formats the given DataFrame to ensure it has the required columns with standardized names.

    The method renames the columns of the given DataFrame to adhere to a predefined mapping, and
    it retains only the required columns specified in the mapping. The resulting DataFrame will
    contain the columns `chromosome`, `start`, `end`, `strand`, `len`, and `seq`. For empty or
    `None` input DataFrame, a new empty DataFrame with the appropriate column names is returned.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing sequence data. It is expected to have columns such as
        'sseqid', 'sstart', 'send', 'sstrand', 'len', and 'sseq'.

    Returns
    -------
    pd.DataFrame
        A DataFrame with standardized column names and only the required columns. If the input
        DataFrame is empty or None, an empty DataFrame with the expected column names is returned.
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=['chromosome', 'start', 'end', 'strand', 'len', 'seq'])

    cols_map = {
        'sseqid': 'chromosome',
        'sstart': 'start',
        'send': 'end',
        'sstrand': 'strand',
        'len': 'len',
        'sseq': 'seq',
    }
    # Select only available required columns and rename
    available = ['sseqid', 'sstart', 'send', 'sstrand', 'len', 'sseq']
    formatted = df[available].rename(columns=cols_map).copy()

    return formatted


def write_gff_from_formatted(df_formatted: pd.DataFrame, gff_path: str) -> None:
    """
    Write a GFF file from a formatted DataFrame.

    This function takes a formatted DataFrame and writes its contents to a file
    in GFF (General Feature Format) format. The function ensures that required
    GFF fields are populated, and missing fields such as strand are handled
    appropriately. The DataFrame is expected to contain specific columns, and
    errors in numeric parsing of the 'start' and 'end' columns will coerce invalid
    entries to NaN, which may raise conversion errors.

    Parameters
    ----------
    df_formatted : pd.DataFrame
        A pandas DataFrame containing data to be converted to GFF format.
        Assumes mandatory columns such as 'chromosome', 'start', and 'end'.
    gff_path : str
        The file path where the GFF output will be written.
    """
    # Handle empty
    if df_formatted is None or df_formatted.empty:
        open(gff_path, "w").close()
        return

    strand_series = df_formatted.get('strand')
    if strand_series is not None:
        strand = strand_series.map({'plus': '+', 'minus': '-'}).fillna('.')
    else:
        strand = pd.Series(['.'] * len(df_formatted))

    ids = [f'BLASTOISE_{i + 1}' for i in range(len(df_formatted))]

    gff_df = pd.DataFrame({
        'seqname': df_formatted['chromosome'].astype(str),
        'source': 'BLASTOISE',
        'feature': '.',
        'start': pd.to_numeric(df_formatted['start'], errors='coerce').astype(int),
        'end': pd.to_numeric(df_formatted['end'], errors='coerce').astype(int),
        'score': '.',
        'strand': strand,
        'frame': '.',
        'attribute': [f'ID={x}' for x in ids],
    })

    gff_df.to_csv(gff_path, sep='\t', index=False, header=False)