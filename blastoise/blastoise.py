import argparse
import os
import shutil
import time  # to measure the time of the program
from datetime import datetime
# import subprocess 

from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.aesthetics import boxymcboxface
from modules.files_manager import fasta_creator, columns_to_numeric
from modules.bedops import bedops_main

# Initiate parser
parser = argparse.ArgumentParser(
    prog="BLASTOISE",
    description="This is a program to search for repetitive sequences in SIDERs elements in Leishmania spp.",
)

# Let's get the user input data
parser.add_argument("-d", "--data", type=str, required=True, help="Path to the input data file")
parser.add_argument("-g", "--genome", type=str, required=True, help="Path to the genome file")

# Parsing the arguments
args = parser.parse_args()

# =============================================================================
# Creating working folder place 
# =============================================================================
# Ask the user for the folder name and data location
folder_name = input("\n\nEnter folder name: ")
data_location = input("Enter path where you want to place all data: ")

data_location = os.path.normpath(data_location)  # Normalize the path to avoid problems with the OS
data_location = os.path.expanduser(data_location)
folder_location = os.path.join(data_location, folder_name)  # Create the folder location
folder_location = os.path.expanduser(folder_location)

# Create the folder with the given name
os.makedirs(folder_location, exist_ok=True)
print(f"{'.'*20} Folder {folder_name} created in {data_location}")

identity_1 = input("Enter the identity for the first BLASTn step: "); identity_1 = int(identity_1)  # user input for the identity #
word_size_param = input("Enter the `word_size` value: "); word_size_param = int(word_size_param)
min_length_param = input("Enter the `min_length` value: "); min_length_param = int(min_length_param)
extend_number_param = input("Enter the `extend_number` value: "); extend_number_param = int(extend_number_param)
first_run = input("Enter the number of the first run: "); first_run = int(first_run)  # user input for the number of runs


# =============================================================================
# Start time
# =============================================================================
start_time = datetime.now()
tic_main = time.perf_counter()  # Start the timer
formatted_start_time = start_time.strftime("%Y %B %d at %H:%M")
print(f"{'.'*20} Program started: {formatted_start_time}")
# =============================================================================
# Take original data used and copy it inside the folder
# =============================================================================
# Expand user's home directory symbol if present
data_path = os.path.expanduser(args.data)
genome_path = os.path.expanduser(args.genome)

# Create subdirectory for original data
original_data_folder = os.path.join(folder_location, "original_data")
os.makedirs(original_data_folder, exist_ok=True)

# Check if the files exist and copy them
if os.path.exists(data_path) and os.path.isfile(data_path):
    shutil.copy(data_path, original_data_folder)
else:
    print(f"Error: The data file '{data_path}' does not exist.")
    exit(1)

if os.path.exists(genome_path) and os.path.isfile(genome_path):
    shutil.copy(genome_path, original_data_folder)
else:
    print(f"Error: The genome file '{genome_path}' does not exist.")
    exit(1)


# print(f"Files copied to {original_data_folder}")
args_data_path = os.path.join(original_data_folder, os.path.basename(args.data))  # save the path so we can use this one instead of the original one
args_genome_path = os.path.join(original_data_folder, os.path.basename(args.genome))  # save the path so we can use this one instead of the original one

# =============================================================================
# First blaster automatization
# =============================================================================
# Create folder for main BLASTN dictionary
blastn_dict_path = os.path.join(folder_location, "dict_data")
os.makedirs(blastn_dict_path, exist_ok=True)
blastn_dict_path_out = os.path.join(blastn_dict_path, os.path.basename(args.genome))

# Create the BLASTn dictionary
blastn_dic(path_input=args.genome, 
           path_output=blastn_dict_path_out)

# Call the first BLASTn
boxymcboxface(message="First BLASTn step initiated")

tic = time.perf_counter()  # Start the timer
first_blaster = blastn_blaster(query_path=args_data_path,
                               dict_path=blastn_dict_path_out, 
                               perc_identity=identity_1,
                               word_size=word_size_param)  # It has the data frame for the first blaster
first_blaster = columns_to_numeric(first_blaster, ["pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"])
toc = time.perf_counter()  # Stop the timer
print(f"1. Initial data:\n",
      f"\t- Data row length: {first_blaster.shape[0]}\n",
      f"\t- Execution time: {toc - tic:0.2f} seconds")
# first_blaster.to_csv(os.path.join(folder_location, "First_Blaster.csv"), index=False, header=True, sep=",")  # Save the data frame to a CSV file

# =============================================================================
# Use BEDOPS to merge the first BLASTn data
# =============================================================================
print("\t- Filtering data:")
tic = time.perf_counter()  # Start the timer
first_blaster_bedops = bedops_main(data_input=first_blaster,
                                   genome_fasta=blastn_dict_path_out,
                                   writing_path_input=folder_location)
toc = time.perf_counter()  # Stop the timers
print(f"\t\t- Data row length: {first_blaster_bedops.shape[0]}\n",
      f"\t\t- Execution time: {toc - tic:0.2f} seconds")

# first_blaster_bedops.to_csv(os.path.join(folder_location, "First_Blaster_BEDOPS.csv"), index=False, header=True, sep=",")  # Save the data frame to a CSV file

# =============================================================================
# Call the second and last BLASTn
# =============================================================================
# 1) Create fasta file for the second BLASTn from first_blaster data frame
fasta_file_path = os.path.join(folder_location, "First_Blaster.fasta")  # Path to the fasta file to create

# Now let's create the fasta file
tic = time.perf_counter()  # Start the timer
fasta_creator(data_input=first_blaster_bedops,
              fasta_output_path=fasta_file_path)
toc = time.perf_counter()  # Stop the timer
print("")
print(f"2. Fasta file creation:\n",
      f"\t- Execution time: {toc - tic:0.2f} seconds")

# Create new folder for all the data
repetitive_blaster_folder = os.path.join(folder_location, "execution_data")
os.makedirs(repetitive_blaster_folder, exist_ok=True)

tic = time.perf_counter()  # Start the timer
repetitive_blaster(data_input=first_blaster_bedops,
                   genome_fasta=blastn_dict_path_out,  # path to the genome dict
                   folder_path=repetitive_blaster_folder,
                   numbering=first_run,
                   start_time=formatted_start_time,
                   identity_1 = identity_1,
                   tic_start = tic_main,
                   word_size=word_size_param,
                   min_length=min_length_param,
                   extend_number=extend_number_param)
toc = time.perf_counter()  # Stop the timer

# =============================================================================
# End time
# =============================================================================
toc_main = time.perf_counter()  # Stop the timer
end_time = datetime.now()
formatted_end_time = end_time.strftime("%Y %B %d at %H:%M")
boxymcboxface(message="END OF THE PROGRAM")
print(f"\t- Execution time: {toc_main - tic_main:0.2f} seconds\n",
      f"\t- Program started: {formatted_start_time}\n",
      f"\t- Program ended: {formatted_end_time}")
