
# =============================================================================
# BASE IMAGE
# =============================================================================
# miniconda3 as the base image - provides Python and conda package manager
# This image is based on Debian and includes conda pre-installed
FROM continuumio/miniconda3

# =============================================================================
# METADATA
# =============================================================================
LABEL maintainer="R. Pacheco <ronnyfph@gmail.com>"
LABEL description="BLASTOISE: A Tool for SIDER Repetitive Sequence Discovery"
LABEL version="0.4.3"

# =============================================================================
# SYSTEM SETUP
# =============================================================================
# Set the working directory inside the container
# All subsequent commands will run from this directory
WORKDIR /app

# Update system packages and install essential tools
RUN apt-get update && \
    apt-get install -y \
        wget \
        build-essential \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# =============================================================================
# CONDA ENVIRONMENT SETUP
# =============================================================================
# Copy the conda environment file first (for better Docker layer caching)
# Docker caches layers, so if environment.yml doesn't change, this layer is reused
COPY environment.yaml .

# Create the conda environment from the environment.yml file
RUN conda env create --name blastoise --file environment.yaml

# Clean conda cache to reduce image size
RUN conda clean -all

# =============================================================================
# BLASTOISE PACKAGE INSTALLATION
# =============================================================================
# Copy entire project
# This includes pyproject.toml, README.md, and all source code
COPY . .

# Install BLASTOISE package in the conda environment
# We use conda run to execute commands in the specific environment
# -n blastoise specifies which environment to use
# pip install -e . installs the package in editable mode
RUN conda run -n blastoise pip install -e .

# =============================================================================
# DIRECTORY SETUP
# =============================================================================
# Create a data directory for mounting external volumes
# This allows users to map their data from the host to the container
RUN mkdir -p /app/data

# Create an output directory for results
RUN mkdir -p /app/output

# Set appropriate permissions (optional, but good for security)
RUN chmod 755 /app/data /app/output

# =============================================================================
# ENVIRONMENT CONFIGURATION
# =============================================================================
# Create a script that activates conda environment and starts bash
RUN echo '#!/bin/bash' > /start.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /start.sh && \
    echo 'conda activate blastoise' >> /start.sh && \
    echo 'exec "$@"' >> /start.sh && \
    chmod +x /start.sh

# Set environment variables
# CONDA_DEFAULT_ENV makes the environment active by default
ENV CONDA_DEFAULT_ENV=blastoise
# Add conda environment to PATH so commands are available
ENV PATH=/opt/conda/envs/blastoise/bin:$PATH

# =============================================================================
# CONTAINER CONFIGURATION
# =============================================================================
# Set the entrypoint to our activation script
ENTRYPOINT ["/start.sh"]

# Set the default command when the container starts
# This starts an interactive bash shell with the conda environment activated
CMD ["/bin/bash"]


# =============================================================================
# USAGE INFORMATION (as comments)
# =============================================================================
# To build this image:
#   docker build -t blastoise:0.4.3 .
#
# To run interactively with data mounting:
#   docker run -it -v /path/to/your/data:/app/data
#
# To run a specific BLASTOISE command:
#   docker run --rm -v /path/to/data:/app/data -v /path/to/output:/app/output blastoise:0.4.3 \
#   conda run -n blastoise blastoise -d /app/data/input.fasta -g /app/data/genome.fasta -o /app/output
#
# To check available commands:
#   docker run --rm blastoise:0.4.3 conda run -n blastoise blastoise --help