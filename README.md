# BLASTOISE: BLAST Oriented to the Identification of SIDER Elements

**BLASTOISE** (BLAST Oriented to the Identification of SIDER Elements) is a Python-based Linux software designed to automate the identification of **SIDER elements** within *Leishmania spp.* species. These retroposons, known for their heterogeneity and degeneration, pose significant challenges in their identification, which BLASTOISE effectively addresses.

BLASTOISE is available as:
- A Docker image on [DockerHub](https://hub.docker.com/repository/docker/rfpacheco/blastoise/general).
- Open-source code on [GitHub](https://github.com/rfpacheco/BLASTOISE).

---

## Key Features
- **Automated SIDER Element Identification**: Tailored for the analysis of *Leishmania spp.* species.
- **Flexible Installation Options**: Supports Docker-based deployment for easy containerized usage, as well as installation via `git clone` and `conda` environments.
- **Robust Libraries**: Leverages a range of powerful libraries and tools, including:
    - [Pandas](https://pandas.pydata.org/) for data manipulation (McKinney et al., 2010)
    - [BEDOPS](https://bedops.readthedocs.io/en/latest/) for genomic feature processing (Neph et al., 2012)
    - [BioPython](https://biopython.org/) for biological computations (Cock et al., 2009)
    - [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/) for sequence alignment (Camacho et al., 2009)

---

## Installation

You can install and use BLASTOISE via docker:

**Using Docker**

To run the Docker image, first you will need to install [Docker engine](https://docs.docker.com/engine/install/ubuntu/) for Linux, including the [Linux postinstall instructions](https://docs.docker.com/engine/install/linux-postinstall/) if you want to omit the `sudo` command. Then, to download **BLASTOISE** software, follow the next instructions:

```bash
docker pull rfpacheco/blastoise
```

## Usage

For the Docker use, follow the [**BLASTOISE** docker use manual](./docs/tutorials/docker_use_manual.md).

And for the Post-Processing step, use [**BLASTOISE** docker post processing use manual](./docs/tutorials/docker_post_processing_manual.md)


## License

BLASTOISE is open-source software, licensed under the [MIT License](https://opensource.org/licenses/MIT).

---

