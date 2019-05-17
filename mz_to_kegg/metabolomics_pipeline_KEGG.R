

# ------------------------------------------ BEGIN EDITABLE REGION ------------------------------------

###################################
# CONFIGURATION:  User preferences
###################################

# FILE LOCATIONS:
#-----------------
# Source files (location of the other R script files)
source('~/github/kroganlab/fluMetabolomics/scripts/mapMasses.R')
source('~/github/kroganlab/fluMetabolomics/scripts/mz_to_kegg/searchDB4KEGG.R')

# path to the results file in 'wide' format
results_file = '~/Box Sync/projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/results/H1N1_lung/H1N1_lung-results-wide.txt'
# path to the evidence file used with MSStats
evidence_file = '~/Box Sync/projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/data/H1N1_lung_evidence.txt'
# The path to the file where you want the results outputted to
out_file = '~/Desktop/metabolomics_debug/H1N1_test_output.txt'

# path to the KEGG database file
kegg_file = '~/Box Sync/projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite_IDs/peak_matching/KEGG_EC_uniprot_mapping.txt'
# path to the Fluomics database file containing the significant hits
flu_file = '~/Box Sync/projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite_IDs/peak_matching/AllSignificantData.txt'  #includes both human and mouse data, but the ertrez id's won't map between the species so ok to use with both


# SEARCH CONSTRAINTS
#--------------------
# log2FC condition (column name) in results to check for significant fold change in
log2FC_condition = 'LungH1.LungPBS_log2FC'
# adjusted p-value condition (column name) in results file to check for significant p-values in
adj_pval_condition = 'LungH1.LungPBS_adj.pvalue'


# any log2 fold change above this value is considered significant (also applies to the negative value in the oposite way)
log2FC = 1
# any p-value below this is considered significant
pvalue = 0.05

THRESH = 0.05  # this is the amount we are willing to let the masses be off for identification +/-
max_weight = 500 # The maximum weight (in Daltons) of metabolites in KEGG to be included in the search

# ------------------------------------------ END EDITABLE REGION ------------------------------------


#########################
# SEARCH KEGG 
#########################

# ~~~~~~~~~~~~~~~
# MAPPING MASSES
# ---------------
# This section will take an MSStats wide formatted file and map the non-rounded masses to the results for identification purposes.

# convert "-" to "." for R purposes
log2FC_condition <- gsub("-", ".", log2FC_condition)
adj_pval_condition <- gsub("-", ".", adj_pval_condition)

# Mapp the non-rounded masses to the peak names ( rounded m.z./time )
to_search = mapMasses(results_file, evidence_file, out_file)
# Keep only significant fold changes in the results
to_search = to_search_all = to_search[(abs(to_search[,log2FC_condition])>log2FC) & (to_search[,adj_pval_condition]<pvalue),]
# prep for the DB search
to_search = to_search[,c('Protein','m.z',log2FC_condition,adj_pval_condition)]


# SEARCH KEGG AGAINST OTHER DATA SETS
# -----------------------------------
# This version is designed to work with the flat files offsite (not live mysql db)
# This section will take the non-rounded masses and search for matches in KEGG. The search process includes including adducts and
#   searches for any matches within a window +/- a THRESHold
# Search the Fluomics OMICS database for genes associated with metabolites
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results = searchDB4KEGG(to_search, kegg_file, flu_file, out_file, THRESH, max_weight)

















