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