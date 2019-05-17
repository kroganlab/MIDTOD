

# ------------------------------------------ BEGIN EDITABLE REGION ------------------------------------

###################################
# CONFIGURATION:  User preferences
###################################

# FILE LOCATIONS:
#-----------------
# Source files (location of the other R script files)
source('~/github/kroganlab/fluMetabolomics/scripts/mapMasses.R')
source('~/github/kroganlab/fluMetabolomics/scripts/mz_to_kegg/searchDB4KEGG.R')
source('~/github/kroganlab/fluMetabolomics/scripts/mz_to_hmdb/searchDB4HMDB.R')
source('~/github/kroganlab/fluMetabolomics/scripts/aggregateResults.R')



# path to the results file in 'wide' format
results_file = '~/Desktop/debug/MKVWH5HTBE_results-wide.txt'
# path to the evidence file used with MSStats
evidence_file = '~/Desktop/debug/MKVWH5HTBE_evidence.txt'
# The path to the file where you want the results outputted to
out_file = '~/Desktop/debug/output.txt'

# SEARCH CONSTRAINTS
#--------------------
# log2FC condition (column name) in results to check for significant fold change in
log2FC_condition = 'H5-12hr-MockH5-12hr_log2FC'
# adjusted p-value condition (column name) in results file to check for significant p-values in
adj_pval_condition = 'H5-12hr-MockH5-12hr_adj.pvalue'

# any log2 fold change above this value is considered significant (also applies to the negative value in the oposite way)
log2FC = 1
# any p-value below this is considered significant
pvalue = 0.05

THRESH = 0.05  # this is the amount we are willing to let the masses be off for identification +/-
max_weight = 500 # The maximum weight (in Daltons) of metabolites in KEGG to be included in the search


# Database Files
#~~~~~~~~~~~~~~~~
# path to the Fluomics database file containing the significant hits
flu_file = '~/github/kroganlab/fluMetabolomics/files/AllSignificantData.txt'
# path to the KEGG database file
kegg_file = '~/github/kroganlab/fluMetabolomics/files/KEGG_EC_uniprot_mapping_20170125.txt'
# path to the HMDB database file
hmdb_file = '~/github/kroganlab/fluMetabolomics/files/HMDB_20150729.txt'

# path to the file containing the mapping between HMDB and Entrez id's (MOUSE or HUMAN?)
#     NOTE: Be sure to select the file for the correct species (mouse/human)
hmdb.entrez_file = '~/github/kroganlab/fluMetabolomics/files/HMDB2entrez_mouse.tsv'
# hmdb.entrez_file = '~/github/kroganlab/fluMetabolomics/files/HMDB2entrez_human.tsv'

# ------------------------------------------ END EDITABLE REGION ------------------------------------


# convert "-" to "." for R purposes
log2FC_condition <- gsub("-", ".", log2FC_condition)
adj_pval_condition <- gsub("-", ".", adj_pval_condition)

# Mapp the non-rounded masses to the peak names ( rounded m.z./time )
to_search = mapMasses(results_file, evidence_file, out_file)
# Keep only significant fold changes in the results
to_search = to_search_all = to_search[(abs(to_search[,log2FC_condition])>log2FC) & (to_search[,adj_pval_condition]<pvalue),]
# prep for the DB search
to_search = to_search[,c('Protein','m.z',log2FC_condition,adj_pval_condition)]


# SEARCH HMDB AGAINST OTHER DATA SETS
# -----------------------------------
# This version is designed to work with the flat files offsite (not live mysql db)
# This section will take the non-rounded masses and search for matches in HMDB. The search process includes including adducts and
#   searches for any matches within a window +/- a THRESHold
# Search the Fluomics OMICS database for genes associated with metabolites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results.hmdb = searchDB4HMDB(to_search, hmdb_file, hmdb.entrez_file, flu_file, out_file, THRESH, max_weight)


# SEARCH KEGG AGAINST OTHER DATA SETS
# -----------------------------------
# This section will take the non-rounded masses and search for matches in KEGG. The search process includes including adducts and
#   searches for any matches within a window +/- a THRESHold
# Search the Fluomics OMICS database for genes associated with metabolites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results.kegg = searchDB4KEGG(to_search, kegg_file, flu_file, out_file, THRESH, max_weight)


# combine the results and aggregate the information
results.agg <- aggregateResults(results.kegg, results.hmdb, out_file)








