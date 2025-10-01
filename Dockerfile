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

# Create an uploads directory for web app file uploads
RUN mkdir -p /app/uploads

# Set appropriate permissions (optional, but good for security)
RUN chmod 755 /app/data /app/output /app/uploads

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
# Expose the FastAPI port
EXPOSE 8000

# Set the entrypoint to our activation script
ENTRYPOINT ["/start.sh"]

# Start the FastAPI app with uvicorn inside the blastoise conda environment by default
CMD ["uvicorn", "webapp.app:app", "--host", "0.0.0.0", "--port", "8000"]

# =============================================================================
# USAGE INFORMATION (as comments)
# =============================================================================
#
# ---------------------------
# --- 1. Get the Image ---
# ---------------------------
#
# === Pull from Docker Hub (Recommended) ===
# Pull the latest version of the image.
#   docker pull rfpacheco/blastoise:latest
#
# You can also pull a specific version for reproducibility. This project's
# current version is 0.4.3.
#   docker pull rfpacheco/blastoise:0.4.3
#
# === Build Locally (Alternative) ===
# Build the image from the project's root directory. You can name it anything
# you like using the `-t` (tag) flag.
#
#   docker build -t my-blastoise-app .
#
# ----------------------------------------
# --- 2. How to Run BLASTOISE ---
# ----------------------------------------
# BLASTOISE offers two main interfaces: a web application and command-line tools.
# The recommended approach depends on your needs.
#
# === Run the Web Application (Recommended for Docker) ===
# The default and primary way to use this Docker image is to run the FastAPI
# web interface for interactive analysis.
# The `--rm` flag removes the container on exit, and `-p` maps the ports.
#
#   docker run --rm -p 8000:8000 rfpacheco/blastoise:latest
#
# Then, open a web browser and navigate to http://localhost:8000
#
# === Run Command-Line (CLI) Tools ===
# For users who primarily want to use the command-line tools, the local
# installation method using Conda (described in README.md) is often more direct.
# However, you can still run CLI commands via Docker, thought is not recommended.
# To do this, you must mount your local data and output
# directories into the container using the `-v` flag.
#
#   -v /path/on/your/host:/path/in/container
#
# Example (main 'blastoise' tool):
#   docker run --rm \
#     -v /path/to/your_data:/app/data \
#     -v /path/to/your_output:/app/output \
#     rfpacheco/blastoise:latest \
#     blastoise -d /app/data/query.fasta -g /app/data/genome.fasta -o /app/output/my_results
#
# Example (getting help):
#   docker run --rm rfpacheco/blastoise:latest blastoise --help
#
# ------------------------------------
# --- 3. Start an Interactive Shell ---
# ------------------------------------
# To get a bash shell inside the container for debugging or manual operations.
# The `-it` flags make the session interactive.
#
#   docker run --rm -it -v /path/to/your_data:/app/data rfpacheco/blastoise:latest bash