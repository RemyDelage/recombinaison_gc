import gffutils
import glob
import os
from os.path import join

# Initialisation of the variables
data_dir = '/home/genouest/cnrs_umr6553/rdelage/data/'
results_dir = '/home/genouest/cnrs_umr6553/rdelage/results/02_genes_fasta/02a_gffutils'
gff_input_files = glob.glob(data_dir + '*.gff')
for files in gff_input_files:
	print(f"gff file found : {files}")
	
	# Create the output files names
	files = files.replace('./', '')
	separator = files.split('_')
	species = join(separator[0], separator[1], separator[2], separator[3]).replace('/', '_')
	output_filename = f"{species}_filtered.gff"
	gff_output_files = os.path.join(results_dir, output_filename) 
	
	# Create the databases
	db  = gffutils.create_db(files, dbfn = 'arabidopsis.db', force = True, keep_order = False, merge_strategy = 'merge', sort_attribute_values = False)
	# Inform the database are created
	print("Database created")

	# Load the databases
	db = gffutils.FeatureDB('arabidopsis.db', keep_order=True)
	print("Database loaded")

	# Open the output files 
	with open(gff_output_files, 'w') as filtered_files:
		for feature in db.features_of_type("gene"):
			if feature.strand != '?':
				gene_name = feature.attributes['Name'][0] if 'Name' in feature.attributes else 'Unknown'
				filtered_files.write(str(feature) + '\n')
				print(f"gene {gene_name} treated")
