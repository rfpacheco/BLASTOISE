import argparse
import os
import shutil
import time  # to measure the time of the program
from datetime import datetime
import subprocess

from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.aesthetics import boxymcboxface
from modules.files_manager import fasta_creator, columns_to_numeric, end_always_greater_than_start
from modules.bedops import bedops_main

# Initiate parser
parser = argparse.ArgumentParser(
    prog='BLASTOISE',
    description='This is a program to search for repetitive sequences in SIDERs elements in Leishmania spp.',
)

# Let's get the user input data
parser.add_argument('-d', '--data', type=str, required=True, help='Path to the input data file')
parser.add_argument('-g', "--genome", type=str, required=True, help='Path to the genome file')
parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output file')
parser.add_argument('-i', '--identity', type=int, required=True, help='Identity percentage for first BLASTn step')
parser.add_argument('-ws', '--word_size', type=int, required=True, help='Word size value for BLASTn')
parser.add_argument('-min', '--min_length', type=int, required=True, help='Minimum length value for filtering')
parser.add_argument('-ext', '--extend', type=int, required=True, help='Extension number value')
parser.add_argument('-lim', '--limit', type=int, required=True, help='Length limit value')

# Parsing the arguments
args = parser.parse_args()

# =============================================================================
# Creating working folder place 
# =============================================================================
# Ask the user for the folder name and data location
# folder_name = input('\n\nEnter folder name: '); folder_name = folder_name.strip()
# data_location = input('Enter path where you want to place all data: ')

# data_location = os.path.normpath(data_location)  # Normalize the path to avoid problems with the OS
# data_location = os.path.expanduser(data_location)
# folder_location = os.path.join(data_location, folder_name)  # Create the folder location
# folder_location = os.path.expanduser(folder_location)

output_path = args.output
output_path = os.path.normpath(output_path)
output_path = os.path.expanduser(output_path)

# Folder name
folder_name = os.path.basename(output_path)
data_location = os.path.dirname(output_path)
folder_location = os.path.join(data_location, folder_name)

# Create the folder with the given name
os.makedirs(folder_location, exist_ok=True)
print(f"{'.'*20} Folder {folder_name} created in {data_location}")

# identity_1 = input('Enter the identity for the first BLASTn step: ')
identity_1 = args.identity

# word_size_param = input('Enter the `word_size` value: ')
word_size_param = args.word_size

# min_length_param = input('Enter the `min_length` value: ')
min_length_param = args.min_length

# extend_number_param = input('Enter the `extend_number` value: ')
extend_number_param = args.extend

# limit_length_param = input('Enter the `limit_length` value: ')
limit_length_param = args.limit

# =============================================================================
# Start time
# =============================================================================
start_time = datetime.now()
tic_main = time.perf_counter()  # Start the timer
formatted_start_time = start_time.strftime('%Y %B %d at %H:%M')
print(f"{'.'*20} Program started: {formatted_start_time}")
# =============================================================================
# Take original data used and copy it inside the folder
# =============================================================================
# Expand user's home directory symbol if present
data_path = os.path.expanduser(args.data)
genome_path = os.path.expanduser(args.genome)

# Create a subdirectory for original data
original_data_folder = os.path.join(folder_location, 'original_data')
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
blastn_dict_path = os.path.join(folder_location, 'dict_data')
os.makedirs(blastn_dict_path, exist_ok=True)
blastn_dict_path_out = os.path.join(blastn_dict_path, os.path.basename(args.genome))

# Create the BLASTn dictionary
blastn_dic(path_input=args.genome, 
           path_output=blastn_dict_path_out)

# Call the first BLASTn
boxymcboxface(message='First BLASTn step initiated')

tic = time.perf_counter()  # Start the timer
first_blaster = blastn_blaster(query_path=args_data_path,
                               dict_path=blastn_dict_path_out, 
                               perc_identity=identity_1,
                               word_size=word_size_param)  # It has the data frame for the first blaster
# Take only the necessary columns
first_blaster = first_blaster[['qseqid', 'sseqid', 'sstart', 'send', 'sstrand']]

toc = time.perf_counter()  # Stop the timer
print(f"1. Initial data:\n",
      f"\t- Data row length: {first_blaster.shape[0]}\n",
      f"\t- Execution time: {toc - tic:0.2f} seconds")

# =============================================================================
# Use BEDOPS to merge the first BLASTn data
# =============================================================================
print('\t- Filtering data:')
tic = time.perf_counter()  # Start the timer
first_blaster_bedops = bedops_main(data_input=first_blaster,
                                   genome_fasta=blastn_dict_path_out)
toc = time.perf_counter()  # Stop the timers
print(f"\t\t- Data row length: {first_blaster_bedops.shape[0]}\n",
      f"\t\t- Execution time: {toc - tic:0.2f} seconds")

# =============================================================================
# Call the Repetitive BLASTn step
# =============================================================================
# Create a new folder for all the data
repetitive_blaster_folder = os.path.join(folder_location, 'execution_data')
os.makedirs(repetitive_blaster_folder, exist_ok=True)

tic = time.perf_counter()  # Start the timer
repetitive_blaster(data_input=first_blaster_bedops,
                   genome_fasta=blastn_dict_path_out,  # path to the genome dict
                   folder_path=repetitive_blaster_folder,
                   numbering=1,
                   start_time=formatted_start_time,
                   identity_1=identity_1,
                   tic_start=tic_main,
                   word_size=word_size_param,
                   min_length=min_length_param,
                   extend_number=extend_number_param,
                   limit_len=limit_length_param)

# Move "blastoise_df.csv" to the main folder
final_data_path = os.path.join(repetitive_blaster_folder, "blastoise_df.csv")
new_final_data_path = os.path.join(folder_location, "blastoise_df.csv")
if os.path.exists(new_final_data_path):
    os.remove(new_final_data_path)
else:
    shutil.move(final_data_path, new_final_data_path)


# Create a new temporary folder 'tmpBlastoise' and move files into it
tmp_folder = os.path.join(folder_location, "tmpBlastoise")
os.makedirs(tmp_folder, exist_ok=True)

for file_or_folder in os.listdir(folder_location):
    # Skip the 'tmpBlastoise' folder itself
    if file_or_folder not in ("original_data", "tmpBlastoise", "blastoise_df.csv"):  # Except for "tmpBlastoise" and "blastoise_df.csv"
        full_path = os.path.join(folder_location, file_or_folder)
        destination_path = os.path.join(tmp_folder, os.path.basename(full_path))
        if os.path.exists(destination_path):  # If the destination exists
            if os.path.isdir(destination_path):
                shutil.rmtree(destination_path)  # Safely remove directory at destination
            else:
                os.remove(destination_path)  # Safely remove a file at destination
        shutil.move(full_path, tmp_folder)



# =============================================================================
# End time
# =============================================================================
toc_main = time.perf_counter()  # Stop the timer
end_time = datetime.now()
formatted_end_time = end_time.strftime("%Y %B %d at %H:%M")
boxymcboxface(message="END OF THE PROGRAM")
print(f"\t- Execution time: {toc_main - tic_main:0.2f} seconds\n",
      f"\t- Program started: {formatted_start_time}\n",
      f"\t- Program ended: {formatted_end_time}\n",
      f"\t- Blastoise final file saved at: {final_data_path}")

# Exit program
boxymcboxface("BLASTING OUT!")
print(f"""
                   o O       o O       o O       o O       o O
                 o | | O   o | | O   o | | O   o | | O   o | | O
               O | | | | O | | | | O | | | | O | | | | O | | | | O
              O-oO | | o   O | | o   O | | o   O | | o   O | | oO-o
             O---o O o       O o       O o       O o      O o  O---o
            O-----O                                           O-----o
            o-----O           ⣠⣴⣾⣶⣿⣿⣶⣶⣶⣿⡟⠀⠀⠀                  o-----O
             o---O          ⣠⣼⣿⣿⡿⣋⣠⣿⣿⣿⣿⣿⡶⢶⣶⣤⣤⣀⣤⣶⣿⣗⡤⠶⢦⡀⠀        o---O
              o-O          ⣤⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠀⠘⢿⣿⣿⣿⡿⢁⣾⡀⠀⢀⡷⠀         o-O
               O           ⠻⣿⠟⠛⠛⠻⠟⠛⠋⢹⣿⣿⣿⢣⡆⠀⠈⠛⠛⠋⠀⢻⣿⣿⣶⠟⣻⣦          O  
              o-O             ⣴⡿⠁   ⢸⣿⣿⢇⣿⠃⠀⣠⣤⣤⣤⣤⣀⠉⠙⢁⣴⡿⠁         o-O    
             o---O       ⢀⣾⣿⣿ ⢠⠀⠀⠀⣰⣤⣿⣿⢋⣬⡄⢀⣾⣿⣿⣿⣿⣿⣿⣧⠀⣿⣯⠀⠀        o---O    
            O-----O  ⢀⣾⣿⣿⣿⣿   ⣠⠻⠿⠿⠿⠿⣛⣵⣿⣿⣧⢸⣿⣿⣿⣿⣿⣿⣿⣿⣄⣿⣿⡆⠀       O-----o    
            O-----O⢀⣾⣿⣿⣿⣿    ⣠⣿⡀⢸⣿⣿⣿⣿⣿⣿⣿⠿⠆⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⠀       O-----o  
            O---⣿⣿⣿⣿⣿⣿⣿  ⢀⣀⣴⣾⣿⣿⡇⣬⣭⣭⣭⣭⣭⣶⣶⣿⣷⡄⢈⣻⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀       O-----o  
            o-⣿⣿⣿⣿⣿      ⠰⢾⣿⣿⣿⣿⡇⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⢐⣛⡻⣿⣿⣿⣿⣿⣿⠻⣿⣿⠀       o-----O
             o-⣿⣿         ⠁⠀⣠⣶⣿⣷⢸⣿⣿⣿⣿⣿⣿⣿⡿⠿⠛⠋⡵⠿⢿⣿⣿⣿⢟⣄⢹⡏⠀        o---O
              o-O          ⣰⣿⣿⣿⣿⣆⢲⣶⣶⣶⣶⣶⣶⣶⣿⢇⣷⣾⣿⡇⣟⣯⣶⣿⣿⡾⠀⠀         o-O  
               O           ⣿⣿⣿⣿⣿⣿⣦⠹⣿⣿⣿⣿⣿⣿⣿⡜⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡾⠋          O
              O-o          ⠘⣿⣿⣿⣿⣿⣿⣷⣬⠉⠿⣛⣻⣿⣯⣥⣹⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀         O-O
             O---o        ⣠⣶⣿⣿⣿⣿⣿⣿⠿⠿⠦⠀⠀⠀⠉⠉⠁⠀⠹⣿⣿⣿⣿⣿⣿⣿⡿⠀⠀        O---o
            O-----o                         ⠀⠹⠿⠛⠿⣿⠟⠛⠛⠀        O-----o
            o-----O                                           o-----O
             o---O o O       o O       o O       o O       o O o---O
              o-Oo | | O   o | | O   o | | O   o | | O   o | | Oo-O
               O | | | | O | | | | O | | | | O | | | | O | | | | O
                 O | | o   O | | o   O | | o   O | | o   O | | o
                   O o       O o       O o       O o       O o     
""")

# Add right to groups and users
subprocess.run(["chmod", "-R", "a+w", data_location], check=True)



