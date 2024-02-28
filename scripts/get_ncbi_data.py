####################### GET_NCBI_DATA #########################################
# This script allows to directly download the genomic fasta and gff files from 
# the NCBI Database only thanks to this accession number
###############################################################################



## Modules importation
import importlib.util
import sys
import ftplib
import os
import wget
import argparse

def is_module_installed(module_name):
    '''
    This finction checks if the necessary modules for run this script are 
    already installed

    Parameters
    ----------
        module_name : the name of the module

    Returns
    -------
    Bool
        True if the module is already installed, False otherwise
    '''
    return importlib.util.find_spec(module_name) is not None

# Function to install a module using pip
def install_module(module_name):
    '''
    Install a Python module using pip.

    Parameters
    ----------
        module_name (str): The name of the module to install.
    
    Returns
    -------
        None
    '''
    try: # Try to intsall the module
        import pip
        pip.main(['install', module_name]) # Module installation
    except Exception as e:
        print(f"Error installing module {module_name}: {e}") 
        # Error message if the module installation failed

# List of required modules
required_modules = ['ftplib', 'os', 'wget', 'argparse']

# Check and install required modules
installed_modules = []
for module in required_modules: # Check if the module is installed
    if not is_module_installed(module): # If the module is not installed
        print(f"The module {module} is not installed. Installing...")
        install_module(module) # Install the module
        installed_modules.append(module) # Add the module to the list of installed modules
    else:
        installed_modules.append(module)

# Creation of the arguments for the command line
## Create the parser
parser = argparse.ArgumentParser(description = 'Download fasta and gff files from the NCBI GeneBank database')
## Add the accession number to the parser
parser.add_argument('-i', '--input', type = str, 
                    help = 'The accession number of the data')

def download_files(accession_number: str):
    '''
    This function allows to download the genomic fasta and gff files from 
    the NCBI database.

    Parameters
    ----------
        accession_number : str
            The accession_number of the data to download

    Returns
    -------
        None
    '''

    accession_num = accession_number.replace("_","") # Remove the underscore
    extension = accession_num[-2:] # Store the extension of the accession number
    accession_no_extension = accession_num[:-2] # Remove the extension
    accession_first_letters = accession_num[0:3] # Keep the first letters
    accession_part = [accession_no_extension[i:i+3] for i in range(0, len(accession_no_extension), 3)]
    # Partitionning the accession number like : AAA/xxx/xxx/xxx/xxx

    # Define the path to the folder containing the data on FTP server
    ftp_folder_path = f"/genomes/all/{'/'.join(accession_part)}/"   
    
    # Connexion to the NCBI FTP server
    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()
    
    ftp.cwd(ftp_folder_path) # Go to the folder that contains the data
    folders = ftp.nlst() # Get the list of the folders of the path.
    
    # Get the folder that matches the accession number
    matching_folders = [folder for folder in folders if folder.startswith(accession_number)]

    # Error message if no matching folders are found
    if not matching_folders: 
        print("No files found. Check the accession number or if data exists on the NCBI GeneBank database.")
        return
        selected_folder = matching_folders[0] # Select the first folder that mathes with the accession number

    ftp.cwd(selected_folder)  # Change to the selected folder  
    files = ftp.nlst()  # Get the list of files in the selected folder

    # Download the files
    for file in files:
        if file.endswith("_genomic.fna.gz") or file.endswith("_genomic.gff.gz"):
            if not file.endswith("_from_genomic.fna.gz"): # Download only the genomic files
                file_url = f"ftp://{ftp.host}{ftp.pwd()}/{file}" # get the URL of the file
                print(f"Downloading genomic file : {file}") # Print the downoaded file name
                wget.download(file_url) # Download the file
                print(" -> Download complete") # Inform the user that the download is complete
    
    # Close the FTP server connexion and not display errors messages
    try:
        ftp.quit() 
    except ftplib.error_temp as e: 
        pass

# Parse the arguments
args = parser.parse_args()
if args.input:
    download_files(args.input)
else:
    print("No accession number provided. Please provide an accession number using the -i or --input option.")
