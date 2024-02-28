#### GET_NCBI_DATA
## This script allows to directly download the genomic fasta and gff files from the NCBI Database 
#  only thanks to this accession number



## Modules importation
import importlib.util
import sys
import ftplib
import os
import wget
import argparse

def is_module_installed(module_name):
    '''
    This finction checks if the necessary modules for run this script are already installed

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
    """
    Install a Python module using pip.

    Parameters
    ----------
        module_name (str): The name of the module to install.
    
    Returns
    -------
        None
    """
    try:
        import pip
        pip.main(['install', module_name])
    except Exception as e:
        print(f"Error installing module {module_name}: {e}")

# List of required modules
required_modules = ['ftplib', 'os', 'wget', 'argparse']

# Check and install required modules
installed_modules = []
for module in required_modules:
    if not is_module_installed(module):
        print(f"The module {module} is not installed. Installing...")
        install_module(module)
        installed_modules.append(module)
    else:
        installed_modules.append(module)

## Creation of the arguments for the command line
parser = argparse.ArgumentParser(description = 'Download fasta and gff files from the NCBI GeneBank database')
parser.add_argument('-i', '--input', type = str, 
                    help = 'The accession number of the data')

def download_files(accession_number: str):

    accession_num = accession_number.replace("_","")
    extension = accession_num[-2:]
    accession_no_extension = accession_num[:-2]
    accession_first_letters = accession_num[0:3]
    accession_part = [accession_no_extension[i:i+3] for i in range(0, len(accession_no_extension), 3)]

    ftp_folder_path = f"/genomes/all/{'/'.join(accession_part)}/"    
    
    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()
    
    ftp.cwd(ftp_folder_path)
    folders = ftp.nlst()
    
    matching_folders = [folder for folder in folders if folder.startswith(accession_number)]
  
    if not matching_folders: 
        print("No files found. Check the accession number or if data exists on the NCBI GeneBank database.")
        return
    
    selected_folder = matching_folders[0]
    ftp.cwd(selected_folder)  # Change to the selected folder
    
    files = ftp.nlst()  # Get the list of files in the selected folder

    for file in files:
        if file.endswith("_genomic.fna.gz") or file.endswith("_genomic.gff.gz"):
            if not file.endswith("_from_genomic.fna.gz"): 
                file_url = f"ftp://{ftp.host}{ftp.pwd()}/{file}"
                print(f"Downloading genomic file : {file}")
                wget.download(file_url)
                print(" -> Download complete")
    
    try:
        ftp.quit()
    except ftplib.error_temp as e:
        pass


args = parser.parse_args()
if args.input:
    download_files(args.input)
else:
    print("No accession number provided. Please provide an accession number using the -i or --input option.")