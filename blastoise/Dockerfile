# Loading image miniconda for conda environments
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Copy the environment.yml file into the container
COPY environment.yml .

# Create the Conda environment from the environment.yml file
RUN conda env create -f environment.yml

# Run data save path for mounting volumes.
RUN mkdir /app/data

# Ensure Conda is properly initialized and available in the container
RUN echo "conda activate blastoise" > ~/.bashrc

# Copy the application code into the container
COPY . .

# Set the default command to start with the environment ready.
## --login will execute the .bashrc scripts.
CMD ["/bin/bash", "--login"]