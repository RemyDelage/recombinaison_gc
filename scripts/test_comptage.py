# Importation of the necessary modules
import argparse
import os

def comptage(file_name:str):
    '''
    This function allows to compute the number of nucleotides contained in
    a fasta file.

    Parameters
    ----------
        file_name : str
            The name of the file on which the computations will be made.
    
    Returns
    -------
    list :
        A list of 3 intergers which corresponds to the number of nucleotides
        [A, C, G, T]
    '''
    # Initialisation of the counters
    a = 0
    c = 0
    g = 0
    t = 0
    _, file_extension = os.path.splitext(file_name) # Checks the files extensions
    if file_extension.lower() not in [".fasta", ".fa", ".fna", ".ffn", ".faa", \
     ".fap", ".fsa"]: # Checks if the file are in FASTA format
        raise ValueError(f"ERROR : File '{file_name}' is not in FASTA format.") # Error message if the files are not in FASTA format
    with open (file_name, "r") as fasta_file: # Open the FASTA file
        for line in fasta_file: # Read each line of the file
            if line.startswith(">") : 
                continue # Ignore the lines which starts with a chevron sign (>)
            for letter in line : # Read each nucleotide
                if letter == "A" or letter == "a":
                    a += 1 # Increment the a counter if the nucleotide is A (or a)
                if letter == "C" or letter == "c":
                    c += 1 # Increment the c counter if the nucleotide is C (or c)
                if letter == "G" or letter == "g":
                    g += 1 # Increment the g counter if the nucleotide is G (or g)
                if letter == "T" or letter == "t":
                    t +=1 # Increment the t counter if the nucleotide is T (or t)
        count=[a,c,g,t] # Cretate a list which contains all the counters
        fasta_file.close()
    
    with open ("nucleotides_counts.txt", "a") as results_file: # Write the results on txt file
        count2 = ",".join(map(str, count)) # Transform the lists into strings
        results_file.write(count2) # Write the file
        results_file.write("\t") # Separate each results counts with tabulation

    return count # Return the list of counters


# Script execution
# Script execution
if __name__ == "__main__":
    # Create the parser to execute this script in command line 
    # and short description of what the script does
    parser = argparse.ArgumentParser(description="Count nucleotides in a FASTA file")
    
    # Create argument of the parser (only file name required)
    parser.add_argument("-f", "--files", type=str, nargs="*", \
    help="Name(s) of the file(s) to process (FASTA format)")
    parser.add_argument("-d", "--directory", type=str,\
     help="Path of the files you wants to study (FASTA format)")
    
    args = parser.parse_args() # Add the arguments to the parser

    file_names = []  # List to store file names

    if args.files:
        file_names.extend(args.files)  # Extend list with file names specified using -f

    if args.directory:
        directory_files = [
            os.path.join(args.directory, f) # Search the directory
            for f in os.listdir(args.directory) # Open the directory and list the files
            if os.path.isfile(os.path.join(args.directory, f)) and 
            f.endswith((".fasta", ".fa", ".fna", ".ffn", ".faa", ".fap", ".fsa")) # Checking if the files are FASTA files
        ]
        file_names.extend(directory_files)  # Extend list with files from directory

    # Error message if the user dose not specify a file or a directory
    if not file_names:
        parser.error("You must specify either a file or a directory.") 

    # Execute the function
    for file_name in file_names:
        try:
            print(comptage(file_name))
        except ValueError as e:
            print(e)
