#' Metabolite identification though orthogonal datasets 
#' @description identify metabolites based on evidences from other
#' dataset generated on the same samples.
#' @author Slim Fourati, \email{sxf279@case.edu}
#' @param resultsFile output of artMS
#' @param evidenceFile output of mapviewer
#' @param species human or mice
#' @param outputDir directory to write output files
#' @return None
#' @examples
#' midtod()
#' @export



# paths must resolve while sourcing, so convert all relative paths outside of a function:

## load subroutines ##
source(file = "scripts/mapMasses.R")
source(file = "scripts/mz_to_kegg/searchDB4KEGG.R")
source(file = "scripts/mz_to_hmdb/searchDB4HMDB.R")
source(file = "scripts/aggregateResults.R")

## orthogonal data files ##
# path to the Fluomics database file containing the significant hits
fluFile <- file.path(getwd(),"files/AllSignificantData.txt")
# path to the KEGG database file
keggFile <- file.path(getwd(),
                      list.files(path       = "files",
                                 pattern    = "KEGG_EC_uniprot_mapping",
                                 full.names = TRUE))
# # path to the HMDB database file
# OLD
hmdbFile <- file.path(getwd(), "files/HMDB_20150729.txt")
# Updated
# hmdb_file <- '~/Downloads/HMDB_endog_20180618.txt'

# path to the file containing the mapping between HMDB and Entrez id's (MOUSE or HUMAN?)
#     NOTE: Be sure to select the file for the correct species (mouse/human)


hmdbEntrezFiles <- list(mouse="files/HMDB2entrez_mouse.tsv",
                        human="files/HMDB2entrez_human.tsv")
hmdbEntrezFiles <- lapply (hmdbEntrezFiles, FUN = function(file)file.path(getwd(),file))



midtod  <- function(resultsFile, evidenceFile, species, outputDir, 
                    remove.infinites=FALSE, orthogonalDataFile = NULL,
                    filterResults = TRUE,
                    log2FC = 1,
                    pvalue = 0.05) {


  ## search constraints ##
  ## theses thresholds do double-duty: 1) they filter the orthogonal data 2) if filterResults == TRUE, they filter the metabolite MS data
  # any log2 fold change above this value is considered significant
  # (also applies to the negative value in the opposite way)
  #log2FC <- 1
  # any p-value below this is considered significant
  #pvalue <- 0.05
  
  # this is the amount we are willing to let the masses be off for
  # identification +/-
  threshold <- 0.05
  # the maximum weight (in Dalton) of metabolites in KEGG to be included
  # in the search
  maxWeight <- 500 


# ------------------------------------------ END EDITABLE REGION ------------------------------------

hmdbEntrezFile <- hmdbEntrezFiles[[species]]
if(is.null(hmdbEntrezFile)){
  stop("\n\nWRONG SPECIES: ", species,": only human or mouse is supported\n")
  
}
  
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

# pre-load flu data here instead of once per each search per loop. This should save a bit of time
# Load significant hits from other OMICS datasets
if (!is.null(orthogonalDataFile)){
  fluFile <- orthogonalDataFile
}
message ("loading omics data from ", fluFile)
flu <- read.delim(fluFile, sep='\t', stringsAsFactors=F, header=T)
# limit by significance and remove unnecessary variables
flu <- flu[!is.na(q_value) & flu$q_value < pvalue &  # check for reasonable q_values
             abs(flu$log2fc) > log2FC & abs(flu$log2fc) < 1e3,  # check for reasonable log2FC -- this removes outliers and infinites
           c('experiment_id','omics_type','condition_2','cell_line','strain','entrez_id','symbol')]   
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
  
  if (filterResults){
    # Keep only significant fold changes in the results
    toSearch <- toSearch[(abs(toSearch[, log2FCcondition]) > log2FC) &
                           (toSearch[, adjPvalCondition] < pvalue), ]
  }
  # Check point for finite values
  if (remove.infinites){
    # preserve rownames because we are about to lose them with lapply
    rowNames <- rownames(toSearch)
    toSearch <- lapply(toSearch,
                       FUN = function(x) {
                         if ("numeric" %in% class(x))
                           x[!is.finite(x)] <- NA
                         return(value = x)
                       })
    #convert back to data.frame and assign row names
    toSearch <- as.data.frame(toSearch, row.names = rowNames)
  }
  # prep table for the DB search and eventual output of important columns
  toSearch$contrast = gsub ("_log2FC", "", log2FCcondition)
  names(toSearch)[which(names(toSearch)==log2FCcondition)] <- "log2FC"
  names(toSearch)[which(names(toSearch)==adjPvalCondition)] <- "adj.pvalue"

  #limit to only columns of interest before checking for complete cases
  toSearch <- toSearch[, c("Protein",
                           "m.z",
                           "contrast",
                           "log2FC",
                           "adj.pvalue")]
  toSearch <- toSearch[complete.cases(toSearch), , drop = FALSE]

  if (nrow(toSearch) > 0) {
    
    # SEARCH KEGG AGAINST OTHER DATA SETS
    # -----------------------------------
    # This section will take the non-rounded masses and search for matches in KEGG. The search process includes including adducts and
    #   searches for any matches within a window +/- a THRESHold
    # Search the Fluomics OMICS database for genes associated with metabolites
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    resultsKEGG <- searchDB4KEGG(input_file = toSearch,
				 kegg_file  = keggFile,
				 flu_file   = flu,#fluFile,
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
				 flu_file         = flu,#fluFile,
				 out_file         = outputFile,
				 THRESH           = threshold,
				 max_weight       = maxWeight)
    
    # combine the results and aggregate the information
    resultsAggregate <- aggregateResults(results.kegg = resultsKEGG,
					 results.hmdb = resultsHMDB,
					 out_file     = outputFile)
    
    message("\n",one," is done\n\n")
  } else{
    message("\nNOT ENOUGH SIGNIFICANT RESULTS FOR THIS COMPARISON\n\n")
  }
}
}
