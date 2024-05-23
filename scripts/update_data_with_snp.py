import pandas as pd
import argparse

def merge_process_snp_info(data_file:str, snp_file:str, output_file:str):
    '''
    This function allows to add the SNP information into the table Species_information.tsv.

    Parameters
    ----------
        - data_file : str
            The TSV file which contains the data which will be updated.
        - snp_file : str
            The TSV file which contains the information about the SNPs.
        - output_file : str
            The output file which contains the updated information with the SNP.

    Returns
    -------
        None
    '''
    # Read the data frame which contains the data of the species.
    data = pd.read_csv(data_file, sep = '\t')
    # Correct the columns names.
    data.columns = ['GENE_ID', 'TRANSCRIPT_ID', 'STRAND', 'START', 'END', 'CHROMOSOME_ID', 'CHROMOSOME_NUM', 'SNP_GENOME', 'SNP_GENE']
    # Drop the columns SNP_GENOME and SNP_GENE (no information).
    data = data.drop(columns = ['SNP_GENOME', 'SNP_GENE'])
    # Remove the lines with the chrloroplastic (Pldt) and mitochondial (MT) chromosomes.
    data = data[~data['CHROMOSOME_NUM'].isin(['Pltd', 'MT'])]
    # Update the GENE_ID to drop the spaces (replace them by underscores)
    data['GENE_ID'] = data['GENE_ID'].str.replace(' ', '_')

    # Read the data frame of the SNPs and rename the columns
    snp = pd.read_csv(snp_file, sep = '\t')
    snp.columns = ['CHROMOSOME_NUM', 'SNP_POSITION', 'REF_ALLELE', 'ALT_ALLELE']

    #List for merge the data
    merged_data = []

    # Get the datas of the chromosome num, the start and the end.
    for i, row in data.iterrows():
        chromosome_num = row['CHROMOSOME_NUM']
        start = row['START']
        end = row['END']
        
        # Filter the data to only keep the SNPs into the genes.
        filtered_snp = snp[(snp['CHROMOSOME_NUM'] == chromosome_num) & 
                           (snp['SNP_POSITION'] >= start) & 
                           (snp['SNP_POSITION'] <= end)]

        # Get all the data.
        for j, snp_row in filtered_snp.iterrows():
            merged_data.append([
                row['GENE_ID'],
                row['TRANSCRIPT_ID'],
                row['STRAND'],
                row['START'],
                row['END'],
                row['CHROMOSOME_ID'],
                row['CHROMOSOME_NUM'],
                snp_row['REF_ALLELE'],
                snp_row['ALT_ALLELE'],
                snp_row['SNP_POSITION']
              ])

    # Create the data frame with the updated data and store it.
    merged_data = pd.DataFrame(merged_data, columns = ['GENE_ID', 'TRANSCRIPT_ID', 'STRAND', 'START', 'END', 'CHROMOSOME_ID', 'CHROMOSOME_NUM', 'REF_ALLELE', 'ALT_ALLELE', 'SNP_POSITION'])
    merged_data.to_csv(output_file, sep = '\t', index = False)
    # Output message to inform the user when the output file is created.
    print(f"Merged data written to {output_file}")


def main():
    '''
    This function allows to create the parser for automatically run this script
    in a single command line.
    '''
    parser = argparse.ArgumentParser(description="Update the genes/SNP informations.")
    parser.add_argument("-i", "--input", type=str, required=True, help="The TSV file containing the first information.")
    parser.add_argument("-snp", "--snp", type=str, required=True, help="The TSV file containing the SNP information from the VCF file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="The TSV file with all the data informations")
    args = parser.parse_args()

    merge_process_snp_info(args.input, args.snp, args.output)

if __name__ == "__main__":
    main()