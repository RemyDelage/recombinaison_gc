import gffutils

gff = '/home/genouest/cnrs_umr6553/rdelage/data/GCF_000001735.4_Arabidopsis_thaliana_genomic.gff'

#Create the database
db  = gffutils.create_db(gff, dbfn = 'arabidopsis_thaliana.db', force = True, keep_order = False, merge_stategy = 'merge', sort_attribute_values = False)

#Inform the database is created
print("Database created")

#Load the database
db = gffutils.FeatureDB('arabidopsis_thaliana.db', keep_order=True)
print("Database loaded")
