# BLASTOISE: BLAST Oriented to the Identification of SIDER Elements

**BLASTOISE** (BLAST Oriented to the Identification of SIDER Elements) is a Python-based Linux software designed to automate the identification of **SIDER elements** within *Leishmania spp.* species. These retroposons, known for their heterogeneity and degeneration, pose significant challenges in their identification, which BLASTOISE effectively addresses.

BLASTOISE is available as:
- Open-source code on [GitHub](https://github.com/rfpacheco/BLASTOISE).
- A Docker image on [DockerHub](https://hub.docker.com/repository/docker/rfpacheco/blastoise/general).

---

## Key Features
- **Automated SIDER Element Identification**: Tailored for the analysis of *Leishmania spp.* species.
- **Flexible Installation Options**: Supports Docker-based deployment for easy containerized usage, as well as installation via `conda`.
- **Robust Libraries**: Leverages a range of powerful libraries and tools, including:
    - [Pandas](https://pandas.pydata.org/) for data manipulation (McKinney et al., 2010)
    - [BEDOPS](https://bedops.readthedocs.io/en/latest/) for genomic feature processing (Neph et al., 2012)
    - [BioPython](https://biopython.org/) for biological computations (Cock et al., 2009)
    - [BLAST](https://blast.ncbi.nlm.nih.gov/) for sequence alignment (Camacho et al., 2009)

---

## Installation

You can install and use BLASTOISE via different methods:

<!-- # TODO: Improve the Docker command, e.g., mount volumes -->
### 1. **Using Docker**
To run the Docker image, ensure you have [Docker](https://www.docker.com/) installed, and execute:
```bash
docker pull rfpacheco/blastoise
docker run -it rfpacheco/blastoise
```

<!-- TODO: improve the conda installation -->
### 2. **Conda Installation**
Use the Conda-based installation approach:
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/).
2. Create the environment and install dependencies:
   ```bash
   conda create -n blastoise_env
   conda activate blastoise_env
   # Additional library installations here (e.g., Pandas, BioPython, etc.)
   ```
   
---

## Usage

1. **Input Preparation**: Provide input genomic data in the required format.
2. **Run the Software**:
   Use command-line arguments to configure BLASTOISE for specific datasets. Example:
   ```bash
   python blastoise.py --input <input-file> --output <output-dir>
   ```
3. **Output Analysis**: The software generates reports and data visualizations for identified SIDER elements.

For advanced usage and configuration, please check the [user guide](https://github.com).

---

## Citation

When using BLASTOISE for your research, please cite the following:


---

## Support and Contributions

- **Issues**: For bugs or feature requests, please open an issue on the [GitHub Issues page](https://github.com).
- **Contributions**: Contributions are welcome! Please follow the contribution guidelines provided in the repository.

---

## License

BLASTOISE is open-source software, licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

## Additional Resources



---