# BLASTOISE Docker Image Guide - blastoise.py

This guide will help you use the **BLASTOISE** software via its Docker image. Follow the steps below to successfully set up and execute **BLASTOISE**.

---

## Checking the Docker Image

After downloading the Docker image, you can verify its presence by running:

```bash
docker images
```

You should see an entry similar to this:

```text
REPOSITORY            TAG       IMAGE ID       CREATED      SIZE
rfpacheco/blastoise   latest    8d62dad40d07   9 days ago   2.88GB
```

---

## Running the BLASTOISE Image

There are two ways to execute the **BLASTOISE** software:

1. **With mounted volumes** (recommended)  
2. **Directly inside the Docker container**

### 1. Running with Mounted Volumes (Recommended)

Mounted volumes allow you to sync a folder from your host machine to the Docker container's execution folder. This method ensures that the data remains easily accessible on the host system.

Run the following command to start the container with your data directory linked:

```bash
docker run -it --rm \
-v <path_to_your_data_files>:/app/data \
rfpacheco/blastoise
```

Replace `<path_to_your_data_files>` with the absolute path to your directory containing the genome and query sequence. This command will map your host directory to `/app/data` inside the **BLASTOISE** container.

Once the container is running, execute the software using:

```bash
python3 blastoise.py \
-d data/<query_fasta_file_path> \
-g data/<genome_fasta_file_path>
```

Replace `<query_fasta_file_path>` and `<genome_fasta_file_path>` as appropriate.

---

### 2. Running Inside the Docker Container

If you prefer to work directly inside the Docker container, you can start the container interactively:

```bash
docker run -it --rm rfpacheco/blastoise
```

From here, inside the container `bash` you will have to copy and move all the data manually, with these types of commands, in another terminal different from the docker running image:

```bash
docker cp <source_path> <container_name OR container_id>:<destination_path>
```

To check the `container_name` OR `container_id` run:
```bash
docker ps --all
```

This method can be tedious so we don't recommend it. The mounted volumes are preferred.


---

## User Input During BLASTOISE Execution

Once the software starts, you will need to provide the following input as prompted:

#### 1. **Folder Name**
The name of the folder where all results will be saved:
```text
Enter folder name:
```

#### 2. **Save Path**
The target directory for saving the folder. If you used mounted volumes, save the data in `/app/data` (mapped to your host directory). For example:
```text
Enter path where you want to place all data:
```

#### 3. **Identity Threshold**
The minimum percentage of sequence similarity required for a match. For SIDER elements, use `60`:
```text
Enter the identity for the first BLASTn step: 
```

#### 4. **BLASTn `word_size`**
The `word_size` parameter defines the length of the short exact matches used to seed BLAST alignments. For SIDER elements, use `15`:
```text
Enter the `word_size` value:
```

#### 5. **Minimum Length**
The shortest acceptable length for a candidate sequence. For SIDER elements, use `100`:
```text 
Enter the `min_length` value:
```

#### 6. **Extension Length**
The minimum length required for a sequence to progress in subsequent iterations of the software:
```text 
Enter the `extend_number` value:
```

#### 7. **Starting Run Number**
The starting number for the software runs. If "1" is entered, the first run will be labeled as `1`:
```text 
Enter the number of the first run:
```

---

## Results

The results will be saved in the destination folder you specified, inside the subfolder `execution_data`, with the output file named:

```text
Last_Data.csv
```

You can find all data generated from the software in this folder.

---

## Example

Below is an example of running **BLASTOISE** with mounted volumes:

1. Start the Docker container:
    ```bash
    docker run -it --rm -v /path/to/data:/app/data rfpacheco/blastoise
    ```

2. Run the software:
    ```bash
    python3 blastoise.py \
    -d data/query.fasta \
    -g data/genome.fasta
    ```

3. Provide the required inputs as prompted (refer to the **User Input** section).

4. Retrieve your results from `/<path>/data/execution_data/Last_Data.csv`.

---

## Additional Notes

- Ensure that your host data folder contains the required `query` and `genome` FASTA files before starting the container.
- If you exit the container, all the data inside the Docker container will be removed (except if you used the mounted volumes).
- The identity threshold, word_size, and minimum/extension lengths depend on your use case and input data requirements. Follow best practices for the data you're working with.

---

Happy Hello World! 🎉