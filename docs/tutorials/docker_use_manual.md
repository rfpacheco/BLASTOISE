# BLASTOISE Docker image use

After downloading the docker image, you can check its existance with:

```bash
docker images
```
where you should get:

```text
REPOSITORY            TAG       IMAGE ID       CREATED      SIZE
rfpacheco/blastoise   latest    8d62dad40d07   9 days ago   2.88GB
```

## How to execute the image

For executing the **BLASTOISE** software, you can do it in two different ways:

- With mounted volumes (recommended)
- Inside the docker container.

## With mounted volumes

With _mounted volumes_ you will sync the folder of interest in your host to the execution folder inside the docker container. For the **blastoise** software you will need:

```bash
docker run -it --rm \
-v <path to your data files>:/app/data \
rfpacheco/blastoise
```

This way you will link your \<path to your data files\> where your genome and query sequence is located to the path /app/data inside the **blastoise** docker container.

The next step is to call the software inside the container:

```bash
python3 blastoise.py \
-d data/<query fasta file name> \
-g data/<genome fasta file name>
```

After doing that the software will start running and it will ask the next questions:

- All the data will be saved inside a folder, so you will need to input the name of the folder:
```text
Enter folder name:
```

- The place where you can to save that folder. Since we mounted your host <path to data> folder to the container data/ folder, you will need to save it there, in data/
```text
Enter path where you want to place all data:
```

- The identity threshold for the minimum percentage of sequence similarity required for a match. For SIDER elements, use a 60 (i.e., 60%) identity.
```text
Enter the identity for the first BLASTn step: 
```

- BLASTn `word_size` parameter, i.e., length of short exact matched used to seed BLAST alignment. For SIDER elements use 15.
```text
Enter the `word_size` value:
```

- The minimum length is the shortest acceptable length for a candidate sequence to be considered. Use 100 for SIDER elements.
```text 
Enter the `min_length` value:
```

- The extension length specifies the minimum length required for a sequence to be considered in subsequent iterations of the software.
```text 
Enter the `extend_number` value:
```

- The number will be used to number all the software runs, the first numbered run will be the one placed here, so if the argument is "1", it will start in 1.
```text 
Enter the number of the first run: 
```

## End Results

The end results will be saved in the defined destiny folder, inside the subfolder "execution_data", with the name "Final_Data.csv".