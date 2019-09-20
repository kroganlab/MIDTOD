midtod  <- function(resultsFile, evidenceFile, species, outputDir) {
###################################
# CONFIGURATION:  User preferences
###################################

# FILE LOCATIONS:
#-----------------
# Source files (location of the other R script files)
source(file = "scripts/mapMasses.R")
source(file = "scripts/mz_to_kegg/searchDB4KEGG.R")
source(file = "scripts/mz_to_hmdb/searchDB4HMDB.R")
source(file = "scripts/aggregateResults.R")


# SEARCH CONSTRAINTS
#--------------------

# any log2 fold change above this value is considered significant (also applies to the negative value in the opposite way)
log2FC <- 1
# any p-value below this is considered significant
pvalue <- 0.05

threshold <- 0.05  # this is the amount we are willing to let the masses be off for identification +/-
maxWeight <- 500 # The maximum weight (in Dalton) of metabolites in KEGG to be included in the search


# # Database Files
# #~~~~~~~~~~~~~~~~
# # path to the Fluomics database file containing the significant hits
fluFile <- "files/AllSignificantData.txt"
# # path to the KEGG database file
keggFile <- "files/KEGG_EC_uniprot_mapping_20170125.txt"
# # path to the HMDB database file
# OLD
hmdbFile <- "files/HMDB_20150729.txt"
# Updated
# hmdb_file <- '~/Downloads/HMDB_endog_20180618.txt'

# path to the file containing the mapping between HMDB and Entrez id's (MOUSE or HUMAN?)
#     NOTE: Be sure to select the file for the correct species (mouse/human)

if (species == "human") {
  hmdbEntrezFile <- "files/HMDB2entrez_human.tsv"
} else if (species == "mouse") {
  hmdbEntrezFile <- "files/HMDB2entrez_mouse.tsv"
} else {
  stop("\n\nWRONG SPECIES: only human or mouse is supported\n")
}

# ------------------------------------------ END EDITABLE REGION ------------------------------------

#########################
# SEARCH HMDB 
#########################
message("METABOLOMICS POST-PROCESSING PIPELINE\n\n")
message("---+ Loading files\n")

# ~~~~~~~~~~~~~~~ 
# MAPPING MASSES
# ---------------
# This section will take an MSStats wide formatted file and map the non-rounded masses to the results for identification purposes.

# Read the results file
dataDF <- read.delim(resultsFile, stringsAsFactors = FALSE)

# Get the column names (comparisons)
comparisons <- names(dataDF)[grepl("_log2FC", names(dataDF))]
comparisons <- gsub(pattern = "_log2FC", replacement = "", comparisons)

message("---+ Ready to process all the relative quantifications:\n")

for (one in comparisons){
  
  # Debugging
  # one <- 'H5THP1Inf3h.H5THP1Mo3h'
  
  message("--------------------------------------------\n")
  message(one," ===> ")
  
  # convert "-" to "." for R purposes
  one <- gsub(pattern = "-", replacement = ".", one)
  
  log2FCcondition <- paste0(one, "_log2FC")
  adjPvalCondition <- paste0(one, "_adj.pvalue")
  
  message(log2FCcondition, " + ", adjPvalCondition,"\n")
  message("--------------------------------------------\n")
  
  outputFilePre <- gsub(pattern = ".txt", replacement = "", basename(path = resultsFile))
  outputFile <- paste0(outputDir, "/", outputFilePre, "-", one, "-midtod.txt")

  message(resultsFile)
  message(evidenceFile)
  message(outputFile)
  
  # Map the non-rounded masses to the peak names ( rounded m.z./time )
  toSearch <- mapMasses(results_file  = resultsFile,
			evidence_file = evidenceFile,
			out_file      = outputFile)
  # Keep only significant fold changes in the results
  toSearch <- toSearch[(abs(toSearch[, log2FCcondition]) > log2FC) &
			 (toSearch[, adjPvalCondition] < pvalue), ]

  # Check point for finite values
  rowNames <- rownames(toSearch)
  toSearch <- lapply(toSearch, FUN = function(x) {
    if ("numeric" %in% class(x))
      x[!is.finite(x)] <- NA
    return(value = x)
  })
  toSearch <- as.data.frame(toSearch, row.names = rowNames)
  toSearch <- toSearch[complete.cases(toSearch), , drop = FALSE]

  if (nrow(toSearch) > 1) {
    # prep for the DB search
    toSearch <- toSearch[, c("Protein",
			     "m.z",
			     log2FCcondition,
			     adjPvalCondition)]
    
    # SEARCH KEGG AGAINST OTHER DATA SETS
    # -----------------------------------
    # This section will take the non-rounded masses and search for matches in KEGG. The search process includes including adducts and
    #   searches for any matches within a window +/- a THRESHold
    # Search the Fluomics OMICS database for genes associated with metabolites
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    resultsKEGG <- searchDB4KEGG(input_file = toSearch,
				 kegg_file  = keggFile,
				 flu_file   = fluFile,
				 out_file   = outputFile,
				 THRESH     = threshold,
				 max_weight = maxWeight)
    
    # SEARCH HMDB AGAINST OTHER DATA SETS
    # -----------------------------------
    # This version is designed to work with the flat files offsite (not live mysql db)
    # This section will take the non-rounded masses and search for matches in HMDB. The search process includes including adducts and
    #   searches for any matches within a window +/- a THRESHold
    # Search the Fluomics OMICS database for genes associated with metabolites
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    resultsHMDB <- searchDB4HMDB(input_file       = toSearch,
				 hmdb_file        = hmdbFile,
				 hmdb.entrez_file = hmdbEntrezFile,
				 flu_file         = fluFile,
				 out_file         = outputFile,
				 THRESH           = threshold,
				 max_weight       = maxWeight)
    
    # combine the results and aggregate the information
    resultsAggregate <- aggregateResults(results.kegg = resultsKEGG,
					 results.hmdb = resultsHMDB,
					 out_file     = outputFile)
    
    message("\n",one,"is done\n\n")
  } else{
    message("\nNOT ENOUGH SIGNIFICANT RESULTS FOR THIS COMPARISON\n\n")
  }
}
}
