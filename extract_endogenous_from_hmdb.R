#! /usr/bin/env Rscript
#
# ------------------
script <- 'extract_endogenous_from_hmdb.R'
# ------------------

# Description

# This script generates a flat file of metabolites obtain from HMDB 
# and filtered by the status "quantified" or "detected".
# To get the info from HMDB requires some work. It is tricky due to the huge size of the 
# HMDB file. Therefore, due to time constraints, I followed a "divide and conquer" strategy:
# 1. Download the "All Metabolites" in XML format from the download page (http://www.hmdb.ca/downloads)
# 2. Create a folder and move the file in
# 3. Run this command: `xml_split -s 100Mb -n 3 -l 2 hmdb_metabolites.xml`. 
# Comment: I am not sure # why it requires the `-n` option, but if it is not added, 
# then it does not split the files properly.As a result, thousands of files are obtain with the 
# info for each individual metabolite in xml format.
# 4. Then REMOVE the `hmdb_metabolites.xml` file from the folder and run

args <- commandArgs(TRUE)

# Number of arguments
n = 1

if (length(args) != n) {
  cat("\nHey! you need to provide,",n,"argument:\n\n")
  cat("\t>",script," split_hmdb_folder
      \n"
  )
  stop("Number of Argument(s) is not right. Bye bye!\n\n")
}

cat("\n",script," is running\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WARNING: make sure that the original <hmdb_metabolites.xml> file is not in the folder!\n\n")

# xml_location <- '~/Downloads/split_hmdb/'
xml_location <- args[1]

if(!grepl("\\/", xml_location)){
  xml_location <- paste0(xml_location,"/")
}

# Generate output filename
t <- Sys.time()
timeStamp <-  strftime(t,"%Y%m%d")
out_file = paste0("HMDB_endog_",timeStamp,".txt")

library(XML)

cat("--+ Ready to process the HMDB split files in the folder\n")
# REmove the file 000.xml which contain the info to do the merge back
file_list = dir(xml_location)[-1]

removeindex <- dir(xml_location)[1] # should be "hmdb_metabolites-000.xml". Just double check
if(grepl("0.xml", removeindex)){
  cat("----- Index file <",removeindex,"> (should have been) ignored\n")  
}else{
  stop("The folder should contain an index file (000.xml) but could not be found. 
       This file was going to be ignored, but if it is not there, it might be a major issue\n\n")
}


all_hmdb = list()
write.table(t(c("id", "name", "weight", "status", "direct_parent", "kegg_id", "chemspider_id", "pubchem_compound_id", "drugbank_id")), out_file, sep='\t', quote=F, append=F, row.names=F, col.names=F)

cat("--+ Time to process each individual file\n\n")
pb <- txtProgressBar(min=0, max=length(file_list), style=3)
a <- character(0)
for(i in 1:length(file_list)){

  # Progress bar
  setTxtProgressBar(pb, i)
  
  # Parse the XML file
  x = xmlTreeParse(paste0(xml_location,file_list[i]))
  xmltop = xmlRoot(x)
  
  # names(xmltop)
  
  # First check whether is predicted
  status = xmlValue(xmltop[["status"]])
  
  # Ready to filter by status
  # if (status == "quantified" | status == "detected" ){
  # if (status == "quantified" | status == "detected" | status == "expected" ){
    # Get the essential values
    hmdb_accession = xmlValue(xmltop[["accession"]])
    metabolite_name = xmlValue(xmltop[["name"]])
    if(length(metabolite_name)<1) metabolite_name = NA
    
    direct_parent = xmlValue(xmltop[["taxonomy"]][["direct_parent"]])
    if(length(direct_parent)<1){
      stop("Wait a minute: there is no DIRECT_PARENT for this entry: ",file_list[i],"\n\n")
    } 
    
    kegg_id = xmlValue(xmltop[["kegg_id"]])
    if(length(kegg_id)<1) kegg_id = NA
    monisotopic_molecular_weight = xmlValue(xmltop[["monisotopic_molecular_weight"]])
    if(length(monisotopic_molecular_weight)<1) monisotopic_molecular_weight = NA
    chemspider_id = xmlValue(xmltop[["chemspider_id"]])
    if(length(chemspider_id)<1) chemspider_id = NA
    pubchem_compound_id = xmlValue(xmltop[["pubchem_compound_id"]])
    if(length(pubchem_compound_id)<1) pubchem_compound_id = NA
    drugbank_id  = xmlValue(xmltop[["drugbank_id"]])
    if(length(drugbank_id)<1) drugbank_id = NA
    
    # Here is where the exact info about Engoneous should be... Don't know how to get there by name
    # origin = xmlValue(xmltop[["ontology"]][[2]][[6]][[1]][[6]][[1]][[1]])
    
    # Check first if the ontology attribute is there, and if is, check whether is Endogenous
    if(any(grepl("ontology", names(xmltop)))){
      if( !identical(a, xmlValue(xmltop[["ontology"]]) ) ){
        if(grepl("Endogenous", xmlValue(xmltop[["ontology"]]))){
          origin = "Endogenous"
        }else{
          origin = NA
        }
        # keep only Endogenous
        if( grepl("Endogenous", origin) ){
          # Old version
          # tmp = data.frame(id = hmdb_accession, name=metabolite_name, kegg_id=kegg_id, origin=origin, weight=monisotopic_molecular_weight)
          # New version
          tmp = data.frame(id = hmdb_accession, name=metabolite_name, weight=monisotopic_molecular_weight, status=status, direct_parent = direct_parent, kegg_id=kegg_id, chemspider_id=chemspider_id, pubchem_compound_id=pubchem_compound_id, drugbank_id=drugbank_id)
          
          write.table(tmp, out_file, sep='\t', quote=F, append=T, row.names=F, col.names=F)
        }
      }
    }
  # } # Option to filter by status
}
close(pb)

cat("----+ ",out_file,"should be ready in the folder\n\n")


