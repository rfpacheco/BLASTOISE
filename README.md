# BLASTOISE: BLAST Oriented to the Identification of SIDER Elements

**BLASTOISE** is a Python-based software tool engineered to automate the identification of SIDER (Short Interspersed 
Degenerated Retroposon) elements within the genomes of Leishmania spp. species. The inherent heterogeneity and 
degeneration of these retroposons present considerable identification challenges, which BLASTOISE is designed to 
overcome.

The software is distributed through several channels: 
- Docker image on [DockerHub](https://hub.docker.com/r/rfpacheco/blastoise).
- Anaconda package on [Anaconda Cloud](https://anaconda.org/rfpacheco/blastoise).
- Open-source code on [GitHub](https://github.com/rfpacheco/BLASTOISE).

## Table of Contents  

* [Use Methods](#use-methods)
* [Installation](#installation)
* [Usage](#usage)
* [Project Structure](#project-structure)
* [Dependencies](#dependencies)
* [License](#license)
* [Citation](#citation)
* [Contact](#contact)

## Use Methods
BLASTOISE offers two distinct operational interfaces to accommodate different user requirements:
* **Web Application:** A comprehensive web application, powered by the FastAPI framework, that provides an intuitive 
    graphical user interface for interactive analysis. This modality is recommended for users deploying the software via 
    Docker.
* **Command-Line Interface (CLI):** A suite of robust command-line tools designed for integration into automated 
    bioinformatics workflows and for users preferring a terminal-based environment.

## Installation
BLASTOISE can be deployed using one of several methods, depending on the desired use case.

### Method 1: Docker (Recommended for Web Application)
This approach is recommended for most users as it simplifies dependency management and ensures a consistent runtime
environment. The Docker container runs on **Windows, macOS, and Linux**.

1.  **Install Docker:** Ensure that **Docker Desktop** (for Windows and macOS) or **Docker Engine** (for Linux) is 
    installed on your system by following the official instructions.
2.  **Get the Docker Image:** You can get the image in two ways:

    *   **Option A: Docker Desktop GUI (Windows/macOS):** Open Docker Desktop, search for `rfpacheco/blastoise` in the 
        images search bar, and click "Pull".
    *   **Option B: Command Line (All Platforms):** Download the official BLASTOISE image from Docker Hub. Users may 
        pull the `:latest` tag for the most recent version or a specific version tag to ensure long-term reproducibility.

    ```bash
    # Pull the latest version
    docker pull rfpacheco/blastoise:latest

    # Or, pull a specific version (e.g., 0.4.3)
    docker pull rfpacheco/blastoise:0.4.3
    ```

### Method 2: Conda Package (Recommended for CLI)
This is the most direct method for installing the command-line tools on a local system.
   1. **Prerequisites:** A functional installation of the Conda or Mamba package manager is required.
   2. **Install the Package:** Execute the following command to install BLASTOISE and all its dependencies directly 
        from the Anaconda Cloud repository. This makes the `blastoise` and `blastoise-sider-filter` commands immediately 
        available in your terminal.

```bash
conda install rfpacheco::blastoise
```

### Method 3: From Source (for Development)
This method facilitates a direct local installation from the source code, which is particularly suitable for development 
purposes.
1. **Prerequisites:** A functional installation of the Conda or Mamba package manager is required.
2. **Clone the Source Repository:**

```bash
git clone https://github.com/rfpacheco/BLASTOISE.git
cd BLASTOISE
```

3. **Create and Activate the Conda Environment:** This command interprets the `environment.yaml` file to provision an 
  environment with all the necessary dependencies.

```bash
# Using Conda
conda env create --file environment.yaml
conda activate blastoise

# Or, using the Mamba package manager for faster resolution
mamba env create --file environment.yaml
mamba activate blastoise
```

4. **Install the BLASTOISE Package:** Install the package in editable mode via pip to complete the setup.

```bash
pip install -e .
```

## Usage
BLASTOISE may be operated through either its web application or its command-line interface.

### 1. Docker: Web Application (Primary Method)

The primary and recommended method for using the Docker image is through the web-based graphical user interface.

*   **Option A: Docker Desktop GUI (Windows/macOS):**
     1.  **Run the image:** In Docker Desktop, go to the "Images" section, find the `rfpacheco/blastoise` image, and 
         click the "Run" button.
     2.  **Configure port:** In the "Optional settings", set the "Host port" to `8000`. This maps the container's port 
         8000 to your local machine's port 8000.
     3.  **Run the container:** Click "Run".
     4.  **Access the Web Interface:** Navigate to `http://localhost:8000` in a web browser.

*   **Option B: Command Line (All Platforms):**
    1.  Execute the following command to start the container. The `-p 8000:8000` flag maps the container's internal port
        to the host machine.
    ```bash
    docker run --rm -p 8000:8000 rfpacheco/blastoise:latest
    ```

> **Note**  
> If the image was built locally under a different tag (e.g., `my-blastoise-app`), that tag should be used
> instead of `rfpacheco/blastoise:latest`.

Once the container is running, you can access the web interface.

2. **Access the Web Interface:** Navigate to the following address in a web browser: http://localhost:8000 (or replace
   localhost with the host IP and 8000 with the port number specified when running the container).
3. **Execute an Analysis:**
   * Upload the requisite input files (query data and genome).
   * Configure the analysis parameters as necessary.
   * Initiate the job by selecting "Run BLASTOISE" or "Run SIDER Filter".
   * Upon completion, the results will be available for download.

### 2. Command-Line Interface (CLI)
The command-line interface is designed for users who require automation, scripting capabilities, or prefer a 
terminal-based workflow. It is the principal mode of interaction for local Conda-based installations.

#### `blastoise` (Main Discovery Tool)
This executable is the primary entry point for executing the SIDER discovery pipeline.

**Basic Example:**
```bash
blastoise -d path/to/query.fasta -g path/to/genome.fasta -o path/to/output_dir
```

**Arguments:**

| Flag   | Full Argument  | Description                                                                             | Default | 
|--------|----------------|-----------------------------------------------------------------------------------------|---------|
| `-h`   | `--help`       | Display the help message and exit.                                                      | -       |
| `-d`   | `--data`       | **Required**. Specifies the path to the input data file (FASTA format query).           | None    |
| `-g`   | `--genome`     | **Required**. Specifies the path to the reference genome file (FASTA format).           | None    |
| `-o`   | `--output`     | **Required**. Specifies the path for the output working directory.                      | None    |
| `-i`   | `--identity`   | Sets the identity percentage for the initial BLASTn step.                               | 60      |
| `-ws`  | `--word_size`  | Defines the word size for the BLASTn algorithm.                                         | 11      |
| `-min` | `--min_length` | Sets the minimum sequence length for filtering operations.                              | 100     |
| `-ext` | `--extend`     | Defines the number of nucleotides for sequence extension.                               | 100     |
| `-lim` | `--limit`      | Sets the length limit that triggers sequence extension.                                 | 1000    |
| `-j`   | `--jobs`       | Specifies the number of jobs for parallel processing (-1 utilizes all available cores). | -1      |

#### `blastoise-sider-filter` (Validation Tool)
This utility is used to validate candidate sequences that have been identified by the main discovery tool.

**Basic Example:**
```bash
blastoise-sider-filter -d candidates.csv -g genome.fasta -rf query.fasta
```

**Arguments:**

| Flag  | Full Argument           | Description                                                                             | Default | 
|-------|-------------------------|-----------------------------------------------------------------------------------------|---------|
| `-h`  | `--help `               | Display the help message and exit.                                                      | -       |
| `-d ` | `--data `               | **Required**. Specifies the path to the input CSV file containing sequence data.        | None    |
| `-g`  | `--genome  `            | **Required**. Specifies the path to the genome FASTA file for BLASTn database creation. | None    |
| `-rf` | `--recaught_file`       | **Required**. Specifies the path to the FASTA file used for recapturing sequences.      | None    |
| `-rt` | `--recaught_threshold ` | Sets the e-value threshold for recapturing sequences.                                   | 0.001   |
| `-ws` | `--word_size`           | Defines the word size for the BLASTn algorithm.                                         | 11      |
| `-e`  | `--evalue`              | Sets the e-value threshold for the initial BLASTn filtering step.                       | 1e-09   |
| `-i`  | `--identity`            | Defines the minimum percentage identity for sequence recapturing.                       | 60      |
| `-ms` | `--min_subjects`        | Sets the minimum number of unique subjects required for a sequence to be accepted.      | 5       |
| `-j`  | `--jobs `               | Specifies the number of jobs for parallel processing (-1 utilizes all available cores). | -1      |

### 3. Advanced: Docker CLI Usage
It is also possible to execute CLI commands directly via Docker by overriding the container's default command. This 
procedure requires mounting local data directories into the container's filesystem using the `-v` flag.

```bash
# Example: Execute the main 'blastoise' tool via Docker
docker run --rm \
 -v /path/to/your_data:/app/data \
 -v /path/to/your_output:/app/output \
 rfpacheco/blastoastoise:latest \
 blastoise -d /app/data/query.fasta -g /app/data/genome.fasta
```

## Project Structure

```
BLASTOISE/
├── blastoise/                  # Main package directory
│   ├── __init__.py            # Package initialization
│   ├── main.py                # Main CLI entry point
│   ├── modules/               # Core functionality modules
│   └── extra/                 # Additional tools (SIDER filter)
├── webapp/                     # FastAPI web application source
├── environment.yaml            # Conda environment specification
├── pyproject.toml              # Project metadata and build configuration
├── meta.yaml                   # Conda recipe for Anaconda distribution
├── Dockerfile                  # Docker image definition
└── README.md                   # This file
```

## Dependencies
The operation of BLASTOISE is contingent upon Python version 3.12 or higher, in addition to several external packages. 
All dependencies are managed automatically through the `environment.yaml` file for Conda-based installations and are 
pre-installed in the official Docker image.
* Core Tools:
  * blast
* Python Libraries:
  * biopython
  * pandas
  * pyranges
  * joblib
  * fastapi (for web application)
  * uvicorn (for web application)

## License
BLASTOISE is open-source software distributed under the terms of the **MIT License**.

## Citation
Should you use BLASTOISE in your research, we kindly request that you cite the following publication:

```
@software{Pacheco_BLASTOISE_2025,
author = {Pacheco, R.},
title = {{BLASTOISE: BLAST Oriented to the Identification of SIDER Elements}},
year = {2025},
publisher = {GitHub},
journal = {GitHub repository},
url = {https://github.com/rfpacheco/BLASTOISE}
}
```

## Contact
For inquiries, bug reports, or contributions, please use the following channels:
* **Author:** R. Pacheco (ronnyfph@gmail.com)
* **Issue Tracker:** GitHub Issues Page
