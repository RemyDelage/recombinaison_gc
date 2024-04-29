import argparse
import gzip
import os
import glob
import vcf
import pandas as pd

def filter_bed_from_vcf(input_path:str, vcf_file:str, output_path:str):
    '''
    This function allows to filter the BED files (compressed) produced by the
    Cactus aligner (Armstrong et al., 2020) for only keep the SNPs presents
    into the VCF file of the species.

    Parameters
    ----------
        input_path: str
            The path where the BED files are stored.

        vcf_file: str
            The species' VCF file containing the SNPs we need to conserve.

        output_path: str
          The path where the produced BED files will be stored.  
        
    Retruns
    -------
        None
            BED files (compressed) with only the sequences contained into
            the VCF file.
    '''
    vcf_pos = set()
    # Read the VCF file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    # Browse each record in the VCF file and add the SNP position to a set.
    for record in vcf_reader:
        vcf_pos.add(record.POS)

    # Unzip and read the BED file (as data frame, panda object)
    for bed_file in glob.glob(os.path.join(input_path, "*.bed.gz")):
        with gzip.open(bed_file, 'rt') as bed:
            bed_df = pd.read_csv(bed, sep='\t', header=None)
        # Filter the data from the BED to only keep the matching sequences with the VCF
        filtered_bed_df = bed_df[bed_df[1].isin(vcf_pos)]

        # Create the output file name
        filename = os.path.basename(bed_file).replace(".bed.gz", "_filtered.bed.gz")
        output_file = os.path.join(output_path, filename)

        # Write the output files
        with gzip.open(output_file, 'wt') as output:
            filtered_bed_df.to_csv(output, sep='\t', header=None, index=None)

        print("Filtered file saved :", output_file)

def main():
    '''
    This function allows to create the parser for automatically run this script
    in a single command line.
    '''
    parser = argparse.ArgumentParser(description="Filter BED files to only keep the sequences present in the VCF file")
    parser.add_argument("--input_path", type=str, help="Path where the BED files are located")
    parser.add_argument("--vcf", type=str, help="The VCF file containing the sequences")
    parser.add_argument("--output_path", type=str, help="Path where the output BED files will be generated")
    args = parser.parse_args()

    filter_bed_from_vcf(args.input_path, args.vcf, args.output_path)

# Execute the script
if __name__ == "__main__":
    main()
