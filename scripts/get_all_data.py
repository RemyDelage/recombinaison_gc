import sys
import argparse

def add_snp_info(species_info_file, snp_file, data_file):
    species_data = dict()
    with open(species_info_file, 'r') as species_file:
        next(species_file)
        for line in species_file:
            # Access to the columns of the data frame
            columns = line.strip().split('\t')
            gene_id = columns[0]
            transcript_id = columns[1]
            strand = columns[2]
            start = columns[3]
            end = columns[4]
            chromosome_id = columns[5]
            chromosome_num = columns[6]
            snp_genome = columns[7]
            snp_gene = columns[8]
            if gene_id not in species_data:
                # Create list for store the genes id
                species_data[gene_id] = []
            species_data[gene_id].append((transcript_id, int(start), int(end), [])) # empty list for SNPs

    with open(snp_file, 'r') as snp:
        next(snp)
        for line in snp:
            columns = line.strip().split('\t')
            chromosome_num = columns[0]
            pos = int(columns[1])  # Convertir la position en entier
            for gene_id, gene_info in species_data.items():
                for transcript_id, start, end, snps in gene_info:
                    if start <= pos <= end:
                        snps.append(line)  # Ajouter la ligne complète du SNP

    with open(data_file, 'w') as data:
        data.write("GENE_ID\tTRANSCRIPT_ID\tSTRAND\tSTART\tEND\tSNP_GENOME\tSNP_GENE\n")
        for gene_id, gene_info in species_data.items():
            for transcript_id, start, end, snps in gene_info:
                data.write(f"{gene_id}\t{transcript_id}\t{strand}\t{start}\t{end}\t{chromosome_id}\t{chromosome_num}\t{','.join(snps)}\t.\n")


def merge_snp_info(data_file, snp_file, output_file):
    species_data = {}
    with open(data_file, 'r') as data_table:
        next(data_table)
        for line in data_table:
            fields = line.strip().split('\t')
            chromosome_num = int(fields[6])
            start = int(fields[3])
            end = int(fields[4])
            species_data[(chromosome_num, start, end)] = fields

    with open(snp_file, 'r') as snp_data:
        next(snp_data)
        for line in snp_data:
            data = line.strip().split('\t')
            chrom = int(data[0])
            start = int(data[1])
            snp = data[2]

            for key in species_data:
                if key[0] == chrom and key[1] <= start <= key[2]:
                    species_data[key][7] = snp

    with open(output_file, 'w') as output:
        output.write("GENE_ID\tTRANSCRIPT_ID\tSTRAND\tSTART\tEND\tCHROMOSOME_ID\tCHROMOSOME_NUM\tSNP_GENOME\tSNP_GENE\n")
        for key in sorted(species_data.keys()):
            output.write('\t'.join(species_data[key]) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add and merge SNP informations')

    parser.add_argument('-i', type=str, help='Data with the informations without SNP')
    parser.add_argument('-snp', type=str, help='The SNP file')
    parser.add_argument('-o', type=str, help='The output file')

    args = parser.parse_args()

    data_file = add_snp_info(args.i, args.snp, "temporary_data.tsv")
    merge_snp_info(data_file, args.snp, args.o)


