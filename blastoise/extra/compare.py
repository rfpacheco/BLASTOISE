from tabnanny import check  # type: ignore
from typing import Dict, Union
from pandas.core.groupby.generic import DataFrameGroupBy
from pandas import DataFrame

import pandas as pd
import argparse
import os
import subprocess


def bedods_contrast(base_df_path: str, contrast_df_path: str, bedops_mode: str) -> DataFrame:
    """
    :param base_df_path: Path to the base dataframe file. What we want to find.
    :param contrast_df_path: Path to the contrast dataframe file. Where we want to find the elements in 'base_df'.
    :param bedops_mode: Mode for the BEDOPS operation; can be 'coincidence', 'opposite', or 'merge'.
    :return: DataFrame with elements that match the specified BEDOPS operation.
    """
    bedops_mode_map = {'coincidence': '--element-of 1',
                       'opposite': '--not-element-of 1',
                       'merge': '--merge'}
    cmd_mode = bedops_mode_map.get(bedops_mode)

    # Check which elements in 'base_df' are inside 'contrast_df'
    cmd_coincidence = f"bedops {cmd_mode} {base_df_path} {contrast_df_path}"
    check_coincidence = subprocess.check_output(cmd_coincidence, shell=True, universal_newlines=True)
    check_coincidence = pd.DataFrame([x.split("\t") for x in check_coincidence.split("\n") if x],
                             columns=["sseqid", "sstart", "send"])
    check_coincidence[['sstart', 'send']] = check_coincidence[['sstart', 'send']].apply(pd.to_numeric)
    return check_coincidence


# Let's define the function since it seems that we will use it more than once
def compare_sequences(df_group1: DataFrameGroupBy, df_group2: DataFrameGroupBy,
                      df_group1_og_data: pd.DataFrame,
                      path: str) -> [Dict[str, Union[str, float]], DataFrame]:
    """
    :param df_group1: First dataframe group comprising sequence data to be compared as true positives.
    :param df_group2: Second dataframe group comprising reference sequence data to check against df_group1.
    :param df_group1_og_data: Original dataset corresponding to df_group1 used for final comparison output.
    :param path: Directory path where intermediate and result files are saved.
    :return: A dictionary summarizing the comparison results and a dataframe with sequences not captured in the second dataframe.
    """
    compare_dict = {}
    total = 0
    not_captured_df = pd.DataFrame()

    for (name1, group1), (name2, group2) in zip(df_group1, df_group2):
        # group1 should be True Positive 0.1.data
        # ------------------------------------------------------------------------------
        path_chr = os.path.join(path, str(name1))
        os.makedirs(path_chr, exist_ok=True)
        # ------------------------------------------------------------------------------
        group1_total = group1[["sseqid", "sstart", "send"]].copy()
        group2_total = group2[["sseqid", "sstart", "send"]].copy()
        # ------------------------------------------------------------------------------
        group1_total.sort_values(by=["sstart", "send"], inplace=True)
        group2_total.sort_values(by=["sstart", "send"], inplace=True)
        # ------------------------------------------------------------------------------
        path_group1_total = os.path.join(path_chr, "group1_total.bed")
        path_group2_total = os.path.join(path_chr, "group2_total.bed")
        # ------------------------------------------------------------------------------
        group1_total.to_csv(path_group1_total, sep="\t", header=False, index=False)  # tabular sep because of bed format
        group2_total.to_csv(path_group2_total, sep="\t", header=False, index=False)  # tabular sep because of bed format

        # ------------------------------------------------------------------------------
        coincidence_result = bedods_contrast(path_group1_total, path_group2_total, "coincidence")
        compare_dict[name1] = [f"{coincidence_result.shape[0]}/{group1_total.shape[0]}",
                               f"{coincidence_result.shape[0] / group1_total.shape[0] * 100:.2f} %"]
        total += coincidence_result.shape[0]

        # ------------------------------------------------------------------------------
        reject_result = bedods_contrast(path_group1_total, path_group2_total, "opposite")
        not_captured_df = pd.concat([not_captured_df, reject_result])

    print(f"""
    There are {total} first DF sequences of {df_group1_og_data.shape[0]} in the second DF:
        - That's {total / df_group1_og_data.shape[0] * 100:.2f}% of the TP 0.data.
        - {df_group1_og_data.shape[0] - total} TP sequences are not in consensus 0.data, 
          which is {100 - total / df_group1_og_data.shape[0] * 100:.2f}%.
    """)
    return compare_dict, not_captured_df



def main():
    parser = argparse.ArgumentParser(description='Process some CSV files.')
    parser.add_argument('--file1', type=str, required=True, help='Path to the software result CSV file.')
    parser.add_argument('--file2', type=str, required=True, help='Path to the bringaud data CSV file.')
    parser.add_argument('--not_name', type=str, required=True, help='Name of the not captured sequences file.')

    args = parser.parse_args()

    df_1 = pd.read_csv(args.file1, sep=',', header=0)  # output from software
    df_2 = pd.read_csv(args.file2, sep=',', header=0)  # bringaud data to compare

    df_1 = df_1[['sseqid', 'sstart', 'send']].copy()  # Take only column 0 ('chrom name'), 1 ('start') and 2 ('end')

    # Create a tmp folder
    df_1_path = os.path.expanduser(args.file1)
    df_1_parent_path = os.path.dirname(df_1_path)
    tmp_folder_path = os.path.join(df_1_parent_path, "tmp_bedops")
    os.makedirs(tmp_folder_path, exist_ok=True)

    # Make the compare:
    ## Frits make a groupby
    df_1_grouped = df_1.groupby("sseqid")
    df_2_grouped = df_2.groupby("sseqid")

    # remember df_2 is the TP data we want to get in the result data df_1
    compare_dict, not_captured_df = compare_sequences(df_2_grouped, df_1_grouped, df_2, tmp_folder_path)
    for key, value in compare_dict.items():
        print(f"{key}: {value[0]} ({value[1]}%)")

    # Remove all the folder `tmp_folder_path`:
    if os.path.exists(tmp_folder_path):
        os.system(f"rm -rf {tmp_folder_path}")

    # Save the 'not_captured_df'
    not_captured_name = f'not_captured-{args.not_name}.csv'
    not_captured_df_path = os.path.join(df_1_parent_path, not_captured_name)
    not_captured_df.to_csv(not_captured_df_path, sep=',', header=True, index=False)

if __name__ == "__main__":
    main()
