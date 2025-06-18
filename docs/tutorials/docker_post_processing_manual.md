# BLASTOISE Docker Image Guide - Post Processing

After running the script `blastoise.py` inside the container, the structure would be the next one inside the docker container:

```text
blastoise docker image
│
└── app/  # Path inside the docker container with all the data
    └── data/  # Mounted container docker path with the host path
        └── <folder_name>/
           ├── original_data/  # Folder where the inputs files will be stored
           ├── tmpBlastoise/   # Folder where all the data needed for blastoise to work, will be placed
           └── blastoise_df.csv # Final data after running blastoise.py

```

## 1. Joining strands data

The `blastoise_df.csv` represents the distributions of elements in both DNA strands, as BLAST search outputs. However, since the _Leishmania spp._ genome is made out of DGC (Directional Gene Clusters), this matter should not be stated yet. This is both strands shall be joined in a common strand:

TODO: REMOVE
```bash
python3 extra/.py \
-f <path_to_blastoise_df.csv> \
-d <genome_fasta_file_path> \
--strand plus  # for the strand, can use "minus" as well
```
After finishing, a file with the name `blastoise_df_merged.csv` will be generated.

```text
blastoise docker image
│
└── app/  
    └── data/  
        └── <folder_name>/
           ├── original_data/  
           ├── tmpBlastoise/   
           ├── blastoise_df.csv 
           └── blastoise_df_merged.csv  # Data with merged strands
```
## 2.SIDER-specific filter & Coordinate correction

### 2.1. SIDER-specific filter

Since the way the software works, not only did it caught SIDER elements, but other repetitive elements. Due that, we consider that every SIDER element $\forall x$, when searched against the whole genome, needs to appear in at least 5 different chromosomes $C(x)$ and to have an expected value < 1.0E-09 $E(x)$:

$$
\forall x \ (C(x) \ \land \ E(x) )
$$



### 2.2. Coordinates correction.

A coordinate correction step is implemented next. This step addresses the potential over-extension of sequences that can occur during the iterative search process while using BLASTn. This coordinate method refines the boundaries of each element, potentially fragmenting sequence si if necessary.


### 2.3. Implementing script

```bash
python3 extra/ \
-f <path_to_blastoise_df_merged.csv> or <path_to_blastoise_df.csv> \
-d <genome_fasta_file_path> \
-ws 15 \  # If you used 15 in blastoise.py, if not, use the 'word_size' number as preferred.
-rf <query_fasta_file_path> 
```

At the end of the process, it will leave with the next datasets:

- **siders_df.csv** --> SIDER elements
- **non_siders_df.csv** --> with the non-SIDER elements

```text
blastoise docker image
│
└── app/  
    └── data/  
        └──  <folder_name>/
           ├── original_data/  
           ├── tmpBlastoise/   
           ├── tmpFilterCoord/  # Needed for the module to work
           ├── blastoise_df.csv 
           ├── blastoise_df_merged.csv
           ├── siders_df.csv  # new file
           └── non_siders_df.csv  # new file
```