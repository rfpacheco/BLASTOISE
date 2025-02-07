from Bio import SeqIO

from blastoise.modules.bedops import bedops_main

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def chromosome_filter(path_input, name):
    """
    Reads a FASTA file and generates a list of chromosome identifiers in a formatted manner.

    With this filter we obtain the labels/titles for each chromosome of our file, e.g., in **Leishmania** case we'll obtain the labels "LinJ.01", "LinJ.02", etc.
    This filter reads all the sequences in a FASTA file. Then with a prefix ``name``, it adds a numbering in the format ``.XX``, being X the numbers in order for each sequences it finds.

    For example, for:

    - ``name = "LinJ"``
    - A fasta file of 36 sequences:

    The output would be:

    .. code-block:: bash

       ['LinJ.01', 'LinJ.02', 'LinJ.03', 'LinJ.04', 'LinJ.05', 'LinJ.06', 'LinJ.07', 'LinJ.08', 'LinJ.09', 'LinJ.10', 'LinJ.11', 'LinJ.12', 'LinJ.13', 'LinJ.14', 'LinJ.15', 'LinJ.16', 'LinJ.17', 'LinJ.18', 'LinJ.19', 'LinJ.20', 'LinJ.21', 'LinJ.22', 'LinJ.23', 'LinJ.24', 'LinJ.25', 'LinJ.26', 'LinJ.27', 'LinJ.28', 'LinJ.29', 'LinJ.30', 'LinJ.31', 'LinJ.32', 'LinJ.33', 'LinJ.34', 'LinJ.35', 'LinJ.36']

    The objective of this function is to be able to correctly name the files since programming languages do not usually admit a non-string format, numbers that start at 0, in this way we can automate their correct labeling, especially for numbers from 01 to 09 .

    .. attention::
       This IDs need to be matched **exactly** with the row[1] from the CSV to filter.

    Reads a FASTA file and generates a list of chromosome identifiers in a formatted manner.

    Args:
        path_input (str): The path to the input FASTA file.
        name (str): A prefix name to be added to each chromosome identifier.

    Returns:
        list: A list of formatted chromosome identifiers.
    """
    max_chr = len(list(SeqIO.parse(path_input, "fasta")))  # It reads the FASTA file and gets the total number of chromosomes.
    chromosome_number = []
    main_list = (list(range(1, max_chr + 1)))  # We index correctly with "+ 1" since Python starts everything in 0
    for number in main_list:
        number = str(number)
        if len(number) == 1:
            chromosome_number.append(name + ".0" + number)  # For it to be the same as the CSV
        else:
            chromosome_number.append(name + "." + number)

    return chromosome_number

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def global_filters_main(data_input, genome_fasta, writing_path, min_length):
    """
    Applies a series of filters to the input data and then processes the data using `bedops_main` function.

    Parameters:
    data_input (DataFrame): The input data to be filtered.
    genome_fasta (str): Path to the genome fasta file required for further processing.
    writing_path (str): Path where the processed data will be saved.
    min_length (int): The minimum length to filter the input data.

    Returns:
    DataFrame: The processed data after applying filters and `bedops_main` function.
    """

    data_filtered = data_input[data_input["length"].astype(int) >= min_length]  # Filter by length; here it was a 100 before

    data_filtered = data_filtered.apply(lambda x: x.replace("-", ""))  # Filter dashes

    if data_filtered.empty:  # Checks if the data is empty. If it is, it will skip the next part of the code
        return  # Skip the next part of the code

    data_bedops = bedops_main(data_input=data_filtered,
                              genome_fasta=genome_fasta,
                              writing_path_input=writing_path)
   
    return data_bedops