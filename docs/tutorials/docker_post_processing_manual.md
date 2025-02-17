# BLASTOISE Docker Image Guide - Post Processing

After running the script `blastoise.py` inside the container, the structure would be the next one inside the docker container:

```text
blastoise docker image
│
└── app/  # Path inside the docker container with all the data
    └── data/  # Mounted container docker path with the host path
        └── <folder_name>/
           └── execution_data/
                └── Last_Data.csv  # Final data after running blastoise.py
```

## 1. Joining strands data

The `Last_Data.csv` represents the distributions of elements in both DNA strands, as BLAST search outputs. However, since the _Leishmania spp._ genome is made out of DGC (Directional Gene Clusters), this matter should not be stated yet. This is both strands shall be joined in a common strand:

```bash
python3 extra/joinStrands.py \
<path_to_Last_Data.csv> \
<genome_fasta_file_path> \
plus  # for the strand, can use "minus" as well
```
After finishing, a file with the name `merged_sequences.csv` will be generated.

```text
blastoise docker image
│
└── app/  
    └── data/  
        └── <folder_name>/
           └── execution_data/
                ├── Last_Data.csv  
                └── merged_sequences.csv # new file
```
## 2.SIDER-specific filter

Since the way the software works, not only did it caught SIDER elements, but other repetitive elements. Due that, we consider that every SIDER element $\forall x$, when searched against the whole genome, needs to appear in at least 5 different chromosomes $C(x)$ and to have an expected value < 1.0E-09 $E(x)$:

$$
\forall x \ (C(x) \ \land \ E(x) )
$$

This is why we apply the following to `merged_sequences.py`:
```bash
python3 extra/true_sider_filter.py \
-f <path_to_merged_secuences.csv> \
-d <genome_fasta_file_path> \
--word-size 15 \  # The number depends on the value used in blastoise.py
--recaught_file <query_fasta_file_path>
```
After executing the script, it leaves with two new files:
- `final_no_data.csv` --> those elements which did not pass the SIDER test.
- `final_yes_data.csv` --> those elements which did pass the SIDER test.

```text
blastoise docker image
│
└── app/  
    └── data/  
        └── <folder_name>/
           └── execution_data/
                ├── Last_Data.csv  
                ├── merged_sequences.csv
                ├── final_no_data.csv  # new file
                └── final_yes_data.csv  # new file
```

## 3. Coordinates correction.

A coordinate correction step is implemented next. This step addresses the potential over-extension of sequences that can occur during the iterative search process while using BLASTn. This coordinate method refines the boundaries of each element, potentially fragmenting sequence si if necessary:

```bash
python3 extra/correct_coor_to_json.py \
-f <path final_yes_data.csv> \
-d <genome_fasta_file_path> \
-o <coordinate_folder_path>  # Place it in the same path as final_yes_data.csv, write the folder name
```
Inside that folder, a file named `main_dict.json` is created with all the needed information for this process.

```text
blastoise docker image
│
└── app/  
    └── data/  
        └── <folder_name>/
           └── execution_data/
                ├── Last_Data.csv  
                ├── merged_sequences.csv
                ├── final_no_data.csv  
                ├── final_yes_data.csv  
                └── <coordinate_folder_path>/ # new directory
                    └── main_dict.json  # new file
```

The next step is using that JSON file to filter the SIDER elements:

```bash
python3 extra/correct_coor_json_true_sider_test.py \
-f <coordinate_folder_path/main_dict.json>  \
-o <coordinate_folder_path> \
-db <genome_fasta_file_path> 
```

At the end of the process, it will output a new JSON file `filtered_data.json` with the needed information.
```text
blastoise docker image
│
└── app/  
    └── data/  
        └── <folder_name>/
           └── execution_data/
                ├── Last_Data.csv  
                ├── merged_sequences.csv
                ├── final_no_data.csv  
                ├── final_yes_data.csv  
                └── <coordinate_folder_path>/ 
                    ├── main_dict.json  
                    └── filtered_data.json  # new file
```

The last step, is getting a CSV files with the sequences from the JSON file:

```bash
python3 extra/correct_coor_json_to_csv.py \
-f  <coordinate_folder_path/filtered_data.json>\
-o <execution_data_path>\
-db <genome_fasta_file_path>\
-ndb <path_final_no_data.csv>
```
In this last step, it will output our final data:

- **positive_database.csv** --> the our SIDER elements
- **negative_database.csv** --> with the no-SIDER elements