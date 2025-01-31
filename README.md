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

### 1. **Using Docker**

To run the Docker image, first you will need to install [Docker engine](https://docs.docker.com/engine/install/ubuntu/) for Linux, including the [Linux postinstall instructions](https://docs.docker.com/engine/install/linux-postinstall/) if you want to omit the `sudo` command. Then, to download **BLASTOISE** software, follow the next instructions:

```bash
docker pull rfpacheco/blastoise
```

### 2. **Conda Installation**

Follow these steps to set up BLASTOISE using Conda:

1. **Install Conda**:  
   - Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/).


2. **Clone the Repository using Git**: 

 ```bash
 git clone https://github.com/rfpacheco/BLASTOISE.git
 ```

3. **Set Up the Conda Environment**:

```bash
conda env create -f blastoise/environment.yml 
conda activate blastoise
```
   
---

## Usage

1. For the Docker use, follow the [**BLASTOISE** docker use manual](./docs/tutorials/docker_use_manual).
2. For the conda-git based installation, follow [**BLASTOISE** git use manual](./docs/tutorials/conda_git_use_manual.md)


---

## Citation


---

## Support and Contributions


---

## License

BLASTOISE is open-source software, licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

