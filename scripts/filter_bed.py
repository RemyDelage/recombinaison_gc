import argparse
import gzip
import os
import glob
import vcf
import pandas as pd

def filter_bed_from_vcf(input_path:str, vcf_file:str, output_path:str):
    vcf_pos = set()
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        vcf_pos.add(record.POS)

    for bed_file in glob.glob(os.path.join(input_path, "*.bed.gz")):
        with gzip.open(bed_file, 'rt') as bed:
            bed_df = pd.read_csv(bed, sep='\t', header=None)

        filtered_bed_df = bed_df[bed_df[1].isin(vcf_pos)]

        filename = os.path.basename(bed_file).replace(".bed.gz", "_filtered.bed.gz")
        output_file = os.path.join(output_path, filename)

        with gzip.open(output_file, 'wt') as output:
            filtered_bed_df.to_csv(output, sep='\t', header=None, index=None)

        print("Fichier filtré enregistré :", output_file)

def main():
    parser = argparse.ArgumentParser(description="Filter BED files to only keep the sequences present in the VCF file")
    parser.add_argument("--input_path", type=str, help="Path where the BED files are located")
    parser.add_argument("--vcf", type=str, help="The VCF file containing the sequences")
    parser.add_argument("--output_path", type=str, help="Path where the output BED files will be generated")
    args = parser.parse_args()

    filter_bed_from_vcf(args.input_path, args.vcf, args.output_path)

if __name__ == "__main__":
    main()
