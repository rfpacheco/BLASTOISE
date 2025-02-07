import argparse
import os
import pandas as pd

from blastoise.modules.blaster import blastn_dic
from blastoise.extra.functions import sider_filter


# ======================================================================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    # noinspection DuplicatedCode
    parser.add_argument("-f", "--file", type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument("-d", "--dict_path", type=str, required=True, help='Path to the genome fasta file.')
    parser.add_argument("-ws", "--word_size", type=str, required=True, help='word_size value of BLASTN')
    parser.add_argument("-rf", "--recaught_file", type=str, required=True, help='Path to the recaught file.')
    return parser.parse_args()
# ======================================================================================================================









# ======================================================================================================================
if __name__ == "__main__":
    # Load all arguments
    args = parse_arguments()
    csv_path = os.path.expanduser(args.file)
    dict_path = os.path.expanduser(args.dict_path)
    word_size = args.word_size
    recaught_file_path = os.path.expanduser(args.recaught_file)

    # Prepare subfolder
    csv_parent_path = os.path.dirname(csv_path)
    folder_path = os.path.join(csv_parent_path, "tmpFilterCoord")  # tmp folder to place all data before deleting

    # Prepare BLASTn dict
    dict_folder_path = os.path.join(folder_path, "blastn_dict")
    os.makedirs(dict_folder_path, exist_ok=True)
    dict_file_path_out = os.path.join(dict_folder_path, os.path.basename(dict_path))
    blastn_dic(path_input=dict_path, path_output=dict_file_path_out)

    # First SIDER filter
    data = pd.read_csv(csv_path, sep=",", header=0)
    yes_data, no_data = sider_filter(
        df=data,
        dict_path=dict_file_path_out,
        folder_path=folder_path,
        word_size=word_size,
        recaught_file=recaught_file_path
    )

    print("Finished Test")





