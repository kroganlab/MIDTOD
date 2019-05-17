# https://bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html

# Purpose: Search the KEGG Exact Masses for matches with the Metabolomics Masses, 
#     and return matching EC Codes. Then map these to uniprot id's

library(reshape2)
library(KEGGREST)
listDatabases()

# get a list of all the coupounds associated with the different weights
compounds1 <- keggFind("compound", 0:300, "exact_mass")
length(compounds1)
compounds2 <- keggFind("compound", 300:600, "exact_mass")
length(compounds2)
compounds <- c(compounds1, compounds2)

# convert to data frame for better searching
compounds <- data.frame(names=names(compounds), compound_id=compounds, stringsAsFactors = F)


# retrieve all the info about each of the compounds. 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# bin all the coupounds into groups of size 10 or fewer. We can pull only 10 at a time from the server, so must stagger the searches.
compounds.bins <- as.numeric( cut(1:dim(compounds)[1], ceiling(dim(compounds)[1]/10) ) )
# tmp <- data.frame(table(compounds.bins), stringsAsFactors = F)
compounds <- data.frame(compounds, bin = compounds.bins, stringsAsFactors = F)

# Loop through each bin
cat("Searching for Compounds from KEGG")
compound.list <- list()
pb <- txtProgressBar(min=1, max=length(unique(compounds$bin)), style=3)
for(bin in unique(compounds$bin)){
  idx <- which(compounds$bin == bin)
  tmp = keggGet(compounds$names[idx])
  
  # pull out all the EC enzyme codes so we can search Uniprot table for matches
  tmp <- lapply(tmp, function(y) {
      if( !is.null(y$ENZYME) ){
          # Warnings are thrown because data.frame is copying the same values in Entry, mass, and weight through to all Enzymes
          return( suppressWarnings( data.frame( y$ENTRY, y$EXACT_MASS, y$MOL_WEIGHT, y$ENZYME, name= paste(y$NAME, collapse=""), PATHWAY=paste(y$PATHWAY, collapse=""), BRITE=paste(y$BRITE, collapse=""), stringsAsFactors = F) ) )
      }
    } 
  ) 
  # combine all the different enzyme lists together
  compound.list[[bin]] <- do.call(rbind, tmp)
  
  setTxtProgressBar(pb, bin)
}
close(pb)

compounds <- do.call(rbind, compound.list)
names(compounds) = gsub("y\\.", "", names(compounds))
names(compounds)[1] = "kegg_id"

# write the results out
# write.table(compounds, '~/Box Sync/db/ontology/fluomics_KEGG_compound_MASS_EC_20170113.txt', sep='\t', quote=F, row.names=F)


# Mapp to Uniprot
#~~~~~~~~~~~~~~~~~~~~
# read in uniprot to EC mapping
ec <- read.delim("~/Box Sync/db/mappings/uniprot2EC_Mouse_Human_20170113.txt", sep='\t', stringsAsFactors = F)
# remove all the cases where there is no EC number
ec <- ec[-which(ec$EC.number == ""),]

# merge the two together
tmp <- merge(compounds, ec[c("Entry","EC.number","Cross.reference..GeneID.")], by.x="ENZYME", by.y="EC.number")
names(tmp)[grep("Cross.reference", names(tmp))] <- "entrez_id"  # rename "Cross.reference..GeneID."


# !!!!!!!!! THERE ARE CASES WHERE THERE ARE MULTIPLE ENTREZ ID'S FOR A SINGLE UNIPROT! WHAT TO DO?
# x <- str_count(tmp$Entrez,";")
# idx <- which(x>1)
# doubs <- tmp[idx,]


# remove the ; at the end of the Entrez id
tmp$entrez_id <- gsub("(^.*)(;$)", "\\1", tmp$entrez_id)


# Remove any Exogenous coupounds via looking at the pathways
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read in list of exogenous KEGG pathways
kegg.exogenous <- read.delim('~/github/kroganlab/fluMetabolomics/files/KEGG_exogenous_pathways_20170125.txt', stringsAsFactors=F, header=F)
query <- paste0('"', paste0(unlist(kegg.exogenous), collapse ="|"), '"')
idx <- grep(query, tmp$PATHWAY)  # ~ 43915/72991

tmp <- tmp[-idx,]

# write the results out
write.table(tmp, '~/github/kroganlab/fluMetabolomics/files/KEGG_EC_uniprot_mapping.txt', sep='\t', quote=F, row.names=F)






































