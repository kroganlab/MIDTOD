# https://bioconductor.org/packages/release/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html

# Purpose: Search the KEGG Exact Masses for matches with the Metabolomics Masses, and return matching uniprots

library(KEGGREST)
listDatabases()


org <- keggList("organism")
head(org)

queryables <- c(listDatabases(), org[,1], org[,2])

keggList("hsa")


query <- keggGet(c("hsa:10458", "ece:Z5100"))

names(query[[1]])
names(query[[2]])

# get a list of all the coupounds associated with the different weights
compounds1 <- keggFind("compound", 0:300, "mol_weight")
length(compounds1)
compounds2 <- keggFind("compound", 300:600, "mol_weight")
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

compound.list <- list()
pb <- txtProgressBar(min=1, max=length(unique(compounds$bin)), style=3)
for(bin in unique(compounds$bin)){
  idx <- which(compounds$bin == bin)
  tmp = keggGet(compounds$names[idx])
  
  # pull out all the EC enzyme codes so we can search Uniprot table for matches
  tmp <- lapply(tmp, function(y) {
      if( !is.null(y$ENZYME) ){
          # Warnings are thrown because data.frame is copying the same values in Entry, mass, and weight through to all Enzymes
          return( suppressWarnings( data.frame( y$ENTRY, y$EXACT_MASS, y$MOL_WEIGHT, y$ENZYME, stringsAsFactors = F) ) )
      }
    } 
  ) 
  # combine all the different enzyme lists together
  compound.list[[bin]] <- do.call(rbind, tmp)
  
  setTxtProgressBar(pb, bin)
}
close(pb)

tmp <- do.call(rbind, compound.list)










temp = data.frame( y$ENTRY, y$EXACT_MASS, y$MOL_WEIGHT, y$ENZYME, stringsAsFactors = F)







unlist(lapply(x, function(y) return(y$ENTRY)))
#










# COMPOUNDS
x = keggGet("C00008")
x[[1]]$ENZYME



# REACTIONS
x = keggGet("R00002")

keggGet("3.6.1.10")


x = keggGet(compounds1)
length(x)
str(x)

head(keggConv("compound", "ncbi-geneid"))
keggConv("hsa", "uniprot:Q06187") #WORKS
keggConv("uniprot", "compound:C00008") #DOESNT




y = names(tmp)
y = gsub("cpd","compound",y)
keggConv("uniprot", names(tmp))


#~~~~~~~~~~~~~~~~

x = org.Hs.egUNIPROT
mapped_genes = mappedkeys(x)
y = unique(as.data.frame(x[mapped_genes])[,2])

y = x[mapped_genes]

y=keys(org.Hs.egUNIPROT)








x = org.Hs.egPATH
mapped_genes = mappedkeys(x)
y = unique(as.data.frame(x[mapped_genes]))












































