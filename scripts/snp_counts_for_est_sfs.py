# -*- coding:utf-8 -*-

import argparse

def snp_counts(input_file: str, output_file: str):
    '''
    This function allows to filter the data file obtain by the 
    nucleotides_counts_for_est_sfs.py script to only keep the SNPs data for 
    using ST-SFS Software (Keightley et al., 2016).

    Parameters
    ----------
        input_file: str
            The TXT file obtain by the nucleotides_counts_for_est_sfs.py script.

        output_file : str
            The TXT file created with only the SNPs that will be used in 
            EST-SFS software.

    Returns
    -------
        None
    '''
    # Read the TXT input file and create the TXT output file
    with open(input_file, 'r') as input, open(output_file, 'a') as output:
        # Browse the imput file
        for line in input:
            # Separate the file into columns
            col = line.strip().split()
            # Set a variable to keep th lines
            keep_line = True
            # Do the computations only for the two last columns (the outgroups)
            for column in col[-2:]:
                for i in column.split(","):
                    # Checks if it is a number
                    if not i.isdigit():  
                        continue
                    # Don't keep the line if the number is greater than 2
                    if int(i) > 2:
                        keep_line = False
                        break
                # Skip the lines where the numbers are greater than 2
                if not keep_line:
                    break
            # Write the keeping lines into the output file
            if keep_line:
                output.write(line)
    # Message for the user when the process is complete
    print(f"File created : {output_file}")

def main():
    '''
    This function allows to create the parser for automatically run this script
    in a single command line.
    '''
    parser = argparse.ArgumentParser(description="Script for only keep the SNP for using EST-SFS Software (Keightley et al., 2016).")
    parser.add_argument("-i", "--input", type=str, required=True, help="The TXT file with all the sequences.")
    parser.add_argument("-o", "--output", type=str, required=True, help="The TXT file with only the SNPs for EST-SFS.")
    args = parser.parse_args()

    snp_counts(args.input, args.output)

# Execute the script :
if __name__ == "__main__":
    main()
