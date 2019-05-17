#! /usr/bin/env Rscript
#
# ------------------
script <- 'metabolomics_pipeline.R'
# ------------------

# Description
# Generate a config.yaml file required by RMSQv3

args <- commandArgs(TRUE)

# Number of arguments
n <- 4

if (length(args) != n) {
  cat("\nHey! you need to provide,",n,"arguments, in THIS EXACT ORDER:\n\n")
  cat("\t>",script,"
      \t\tresults_file-wide.txt
      \t\tevidence.txt
      \t\toutput_folder_name (e.g: output_midtod_20180508 )
      \t\tspecie [human or mouse]
      \n"
  )
  stop("bye bye!\n\n")
}

# ------------------------------------------ BEGIN EDITABLE REGION ------------------------------------

# # # DEBUG
# setwd('~/sourcecode/artms/metabolomics/results/')
# # 
# # path to the results file in 'wide' format
# results_file <- 'metabolomics-results-wide.txt'
# # path to the evidence file used with MSStats
# # evidence_file <- '/Users/shrivergroup/Documents/MSstat_Projects/Celina/Macrophage/H5/MSstats_v3/MKVWalignMDMpos_evidence.txt'
# evidence_file <- '../H5THP1align_05292018-evidence.txt'
# # The path to the file where you want the results outputted to
# # out_file <- '/Users/shrivergroup/Documents/MSstat_Projects/Celina/Macrophage/H5/MSstats_v3/MKVWalignMDMpos_18h_results.txt'
# out_folder <- 'output_midtod_debugging'
# #  Specie!
# specie <- 'human'

# resutls file in wide format
results_file <- args[1]
# path to the evidence file used with MSStats
evidence_file <- args[2]
# The path to the file where you want the results outputted to
out_folder <- args[3]
# Choose the right hmdb
specie <- args[4]

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


# SEARCH CONSTRAINTS
#--------------------

# any log2 fold change above this value is considered significant (also applies to the negative value in the oposite way)
log2FC <- 1
# any p-value below this is considered significant
pvalue <- 0.05

THRESH <- 0.05  # this is the amount we are willing to let the masses be off for identification +/-
max_weight <- 500 # The maximum weight (in Daltons) of metabolites in KEGG to be included in the search


# # Database Files
# #~~~~~~~~~~~~~~~~
# # path to the Fluomics database file containing the significant hits
flu_file <- '~/github/kroganlab/fluMetabolomics/files/AllSignificantData.txt'
# # path to the KEGG database file
kegg_file <- '~/github/kroganlab/fluMetabolomics/files/KEGG_EC_uniprot_mapping_20170125.txt'
# # path to the HMDB database file
# OLD
hmdb_file <- '~/github/kroganlab/fluMetabolomics/files/HMDB_20150729.txt'
# Updated
# hmdb_file <- '~/Downloads/HMDB_endog_20180618.txt'

# path to the file containing the mapping between HMDB and Entrez id's (MOUSE or HUMAN?)
#     NOTE: Be sure to select the file for the correct species (mouse/human)

if(specie == "human"){
  hmdb.entrez_file <- '~/github/kroganlab/fluMetabolomics/files/HMDB2entrez_human.tsv'  
}else if(specie == "mouse"){
  hmdb.entrez_file <- '~/github/kroganlab/fluMetabolomics/files/HMDB2entrez_mouse.tsv'
}else{
  stop("\n\nWRONG SPECIE: only human or mouse is supported\n")
}

# ------------------------------------------ END EDITABLE REGION ------------------------------------

#########################
# SEARCH HMDB 
#########################
cat("METABOLOMICS POST-PROCESSING PIPELINE\n\n")
cat("---+ Loading files\n")

# ~~~~~~~~~~~~~~~ 
# MAPPING MASSES
# ---------------
# This section will take an MSStats wide formatted file and map the non-rounded masses to the results for identification purposes.

# Read the results file
datadf <- read.delim(results_file, stringsAsFactors=F)

# Get the column names (comparisons)
comparisons <- names(datadf)[grepl("_log2FC", names(datadf))]
comparisons <- gsub("_log2FC","", comparisons)

cat("---+ Ready to process all the relative quantifications:\n")

for (one in comparisons){
  
  # Debugging
  # one <- 'H5THP1Inf3h.H5THP1Mo3h'
  
  cat("--------------------------------------------\n")
  cat(one," ===> ")
  
  # convert "-" to "." for R purposes
  one <- gsub("-", ".", one)
  
  log2FC_condition <- paste0(one,"_log2FC")
  adj_pval_condition <- paste0(one,"_adj.pvalue")
  
  cat(log2FC_condition," + ", adj_pval_condition,"\n")
  cat("--------------------------------------------\n")
  
  output_file_pre <- gsub(".txt", "", results_file)
  out_file <- paste0(out_folder, "/", output_file_pre, "-", one, "-midtod.txt")
  
  # Mapp the non-rounded masses to the peak names ( rounded m.z./time )
  to_search <- mapMasses(results_file, evidence_file, out_file)
  # Keep only significant fold changes in the results
  
  to_search <- to_search[(abs(to_search[,log2FC_condition])>log2FC) & (to_search[,adj_pval_condition]<pvalue),]

  # Check point for finite values
  rNames <- rownames(to_search)
  to_search <- lapply(to_search, FUN = function(x) {
    if ("numeric" %in% class(x))
      x[!is.finite(x)] <- NA
    return(value = x)
  })
  to_search <- as.data.frame(to_search, row.names = rNames)
  to_search <- to_search[complete.cases(to_search), , drop = FALSE]

  if(dim(to_search)[1] > 0){
    # prep for the DB search
    to_search <- to_search[,c('Protein','m.z',log2FC_condition,adj_pval_condition)]
    
    # SEARCH KEGG AGAINST OTHER DATA SETS
    # -----------------------------------
    # This section will take the non-rounded masses and search for matches in KEGG. The search process includes including adducts and
    #   searches for any matches within a window +/- a THRESHold
    # Search the Fluomics OMICS database for genes associated with metabolites
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    results.kegg <- searchDB4KEGG(to_search, kegg_file, flu_file, out_file, THRESH, max_weight)

    # SEARCH HMDB AGAINST OTHER DATA SETS
    # -----------------------------------
    # This version is designed to work with the flat files offsite (not live mysql db)
    # This section will take the non-rounded masses and search for matches in HMDB. The search process includes including adducts and
    #   searches for any matches within a window +/- a THRESHold
    # Search the Fluomics OMICS database for genes associated with metabolites
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    results.hmdb <- searchDB4HMDB(to_search, hmdb_file, hmdb.entrez_file, flu_file, out_file, THRESH, max_weight)
    
    # combine the results and aggregate the information
    results.agg <- aggregateResults(results.kegg, results.hmdb, out_file)
    
    cat("\n",one,"is done\n\n")
  }else{
    cat("\nNOT ENOUGH SIGNIFICANT RESULTS FOR THIS COMPARISON\n\n")
  }
}

cat(script, "is done. Have a nice day! \n\n")
