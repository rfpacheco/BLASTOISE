import argparse
import pandas as pd
import os


# =============================================================================
# Defining argparser arguments
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Create a GFF file from a CSV file.')
    parser.add_argument("-csv", "--csv_file", type=str, required=True, help="Path to the CSV file with the data.")
    parser.add_argument("--source", type=str, required=False, default=".", help="Source of the data.")
    parser.add_argument("--feature", type=str, required=False, default=".", help="Feature of the data.")
    parser.add_argument("--attribute", type=str, required=False, default=".", help="Column name of the attribute.")
    return parser.parse_args()


# =============================================================================
# Main function
# =============================================================================
def csv_to_gff(csv_file, source=".", feature=".", attribute="."):
    csv_file = os.path.expanduser(csv_file)
    csv_file_name = os.path.basename(csv_file).replace(".csv", "")
    csv_folder = os.path.dirname(csv_file)

    # Read the CSV file
    data = pd.read_csv(csv_file, sep=",", header=0)

    # Strand input
    if "strand" in data.columns:
        data.loc[data["strand"] == "minus", "strand"] = "-"
        data.loc[data["strand"] == "plus", "strand"] = "+"

    # Prepare GFF data
    pre_gff = []
    column_names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    for index, row in data.iterrows():
        pre_gff.append(
            {
                "seqname": row["chromosome"],
                "source": source,
                "feature": feature,
                "start": row["start"],
                "end": row["end"],
                "score": ".",
                "strand": row["strand"] if "strand" in row and pd.notna(row["strand"]) else ".",
                "frame": '.',
                "attribute": f"ID={row['qseqid'] if 'qseqid' in row else index}_{hash(index) + 1}_{data.shape[0]}"
            }
        )

    # Create the GFF file
    final_gff_df = pd.DataFrame(pre_gff, columns=column_names)
    # File path  
    gff_save_path = os.path.join(csv_folder, f"{csv_file_name}.gff")

    # Save the GFF file
    final_gff_df.to_csv(gff_save_path, sep='\t', index=False, header=False)


if __name__ == '__main__':
    args = parse_arguments()
    csv_to_gff(args.csv_file, args.source, args.feature, args.attribute)