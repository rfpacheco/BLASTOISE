# BLASTOISE: BLAST Oriented to the Identification of SIDER Elements

**BLASTOISE** (BLAST Oriented to the Identification of SIDER Elements) is a Python-based Linux software designed to automate the identification of **SIDER elements** within *Leishmania spp.* species. These retroposons, known for their heterogeneity and degeneration, pose significant challenges in their identification, which BLASTOISE effectively addresses.

**Version:** 0.4.2

BLASTOISE is available as:
- A Docker image on [DockerHub](https://hub.docker.com/repository/docker/rfpacheco/blastoise/general)
- Open-source code on [GitHub](https://github.com/rfpacheco/BLASTOISE)

---

## Key Features

- **Automated SIDER Element Identification**: Tailored for the analysis of *Leishmania spp.* species
- **Flexible Installation Options**: Supports Docker-based deployment for easy containerized usage, as well as installation via `git clone` and `conda` environments
- **Robust Libraries**: Leverages a range of powerful libraries and tools, including:
    - [Pandas](https://pandas.pydata.org/) for data manipulation (McKinney et al., 2010)
    - [PyRanges](https://pyranges.readthedocs.io/) for genomic feature processing
    - [BioPython](https://biopython.org/) for biological computations (Cock et al., 2009)
    - [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) for sequence alignment (Camacho et al., 2009)
- **Additional Tools**:
    - Coordinate Corrector: Adjust genomic coordinates
    - SIDER Filter: Filter and process SIDER elements

---

## Installation

BLASTOISE can be installed and used in multiple ways:

### Using Docker (Recommended)

To run the Docker image, first you will need to install [Docker engine](https://docs.docker.com/engine/install/ubuntu/) for Linux, including the [Linux postinstall instructions](https://docs.docker.com/engine/install/linux-postinstall/) if you want to omit the `sudo` command. Then, to download **BLASTOISE** software, follow these instructions:

```bash
docker pull rfpacheco/blastoise:0.4.2
```

### Using Conda (Alternative)

1. Clone the repository:
   ```bash
   git clone https://github.com/rfpacheco/BLASTOISE.git
   cd BLASTOISE
   ```

2. Create and activate the conda environment:
   ```bash
   conda env create --name blastoise --file environment.yaml
   conda activate blastoise
   ```

3. Install the package:
   ```bash
   pip install -e .
   ```

---

## Usage

### Docker Usage

To run BLASTOISE using Docker:

```bash
# Run interactively with data mounting
docker run -it -v /path/to/your/data:/app/data rfpacheco/blastoise:0.4.2

# Run a specific BLASTOISE command
docker run --rm -v /path/to/data:/app/data -v /path/to/output:/app/output rfpacheco/blastoise:0.4.2 \
conda run -n blastoise blastoise -d /app/data/input.fasta -g /app/data/genome.fasta -o /app/output

# Check available commands
docker run --rm rfpacheco/blastoise:0.4.2 conda run -n blastoise blastoise --help
```

### Command Line Usage (when installed via conda)

```bash
# Main BLASTOISE command
blastoise -d input.fasta -g genome.fasta -o output_directory

# Additional tools
blastoise-coordinate-corrector --help
blastoise-sider-filter --help
```

---

## Project Structure

```
BLASTOISE/
├── blastoise/                  # Main package directory
│   ├── __init__.py            # Package initialization
│   ├── main.py                # Main entry point
│   ├── modules/               # Core functionality modules
│   │   ├── aesthetics.py      # UI and output formatting
│   │   ├── blaster.py         # BLAST integration
│   │   ├── compare.py         # Sequence comparison
│   │   ├── files_manager.py   # File handling
│   │   ├── filters.py         # Filtering mechanisms
│   │   ├── genomic_ranges.py  # Genomic range operations
│   │   ├── seq_identifier.py  # Sequence identification
│   │   ├── seq_modifier.py    # Sequence modification
│   │   └── strand_location.py # Strand location utilities
│   └── extra/                 # Additional tools
│       ├── coordinate_corrector.py  # Coordinate correction tool
│       ├── sider_filter.py          # SIDER filtering tool
│       └── utils/                   # Utility functions
├── docs/                      # Documentation
├── environment.yaml            # Conda environment specification
├── pyproject.toml             # Project metadata and build configuration
├── Dockerfile                 # Docker image definition
└── README.md                  # This file
```

---

## Dependencies

BLASTOISE requires Python 3.12 or higher and depends on several libraries. Most Python libraries are installed automatically via the conda environment (and pip),

- pandas
- pyranges
- biopython
- blast
- Other dependencies specified in environment.yaml

---

## License

BLASTOISE is open-source software, licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

## Citation

If you use BLASTOISE in your research, please cite:

```
Pacheco, R. (2023). BLASTOISE: BLAST Oriented to the Identification of SIDER Elements. 
GitHub repository: https://github.com/rfpacheco/BLASTOISE
```

---

## Contact

For questions, issues, or contributions, please contact:
- R. Pacheco (ronnyfph@gmail.com)
- [GitHub Issues](https://github.com/rfpacheco/BLASTOISE/issues)


## BLASTOISE Graphical User Interface (Tkinter)

BLASTOISE now includes a simple Tkinter-based GUI so you can run the pipeline without using the command line.
It reuses the same core pipeline as the CLI and writes the same outputs (CSV and GFF).

### How to launch the GUI

You can install the project in two ways. Pick one and then launch the GUI.

1) Conda build/install

    # Build the package
    conda build -c conda-forge -c bioconda .

    # Install from local cache
    conda install --use-local -c bioconda -c conda-forge blastoise

    # Launch the GUI
    blastoise-gui

2) Conda environment + pip editable

    # Create environment with conda dependencies
    conda create -n blastoise python=3.12 blast -c bioconda -c conda-forge
    conda activate blastoise

    # Install your package with pip (editable mode)
    pip install -e .

    # Launch the GUI
    blastoise-gui
    # or
    python -m blastoise.gui

Notes:
- BLAST+ tools must be available in PATH (provided by the conda package `blast`).
- Tkinter is part of Python, but on conda it’s provided by the `tk` package (included in the conda recipe). If you are
  outside conda and get errors about tkinter, install your OS’s python3-tk package.

### Using the GUI

1. Paths section
   - Input Data (FASTA): The input query sequences file (e.g., FASTA). Equivalent to CLI `--data`.
   - Reference Genome (FASTA): The genome FASTA to search against. Equivalent to CLI `--genome`.
   - Output Directory: Folder where BLASTOISE will create a workspace and write results. Equivalent to CLI `--output`.

2. Parameters
   - Identity (%): Identity threshold for BLASTn in the first step. Default 60. CLI: `--identity`.
   - Word Size: BLASTn word size. Default 11. CLI: `--word_size`.
   - Min Length: Minimum length to keep alignments/elements. Default 100. CLI: `--min_length`.
   - Extend (nt): Number of nucleotides used in extension steps. Default 100. CLI: `--extend`.
   - Length Limit (nt): Maximum length that triggers/limits extension. Default 1000. CLI: `--limit`.
   - Jobs (-1=all): Parallel jobs. -1 uses all processors. CLI: `--jobs`.

3. Run/Cancel
   - Click Run to start. The process can take minutes to hours depending on genome size and parameters.
   - A live Log panel shows status messages from the pipeline and will also be saved to `blastoise.log` in the
     selected output directory.
   - When finished, you will get a popup dialog with the paths to the generated CSV and GFF files.
   - Cancel will request cancellation, but long-running external tasks (e.g., BLAST) may need to finish their current
     step before the app can close.

### Outputs

- CSV: BLASTOISE--<input_basename>--<genome_basename>.csv
- GFF: BLASTOISE--<input_basename>--<genome_basename>.gff
- Log: blastoise.log in the selected output directory.

### Troubleshooting

- If the GUI does not start and you see errors about display or Tk, ensure you are running in a desktop session and
  that the `tk` package is installed in your environment.
- If BLAST commands are not found, verify that the `blast` package is installed in the current conda env and
  `blastn` is on your PATH.
- For very large genomes, consider increasing available memory and using fewer Jobs if you notice instability.
