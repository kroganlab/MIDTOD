

# Aggregate the results from the KEGG and HMDB search together
aggregateResults <- function(results.kegg, results.hmdb, out_file){
  cat("Combining HMDB and KEGG serch results...\n")
  # fix some of the names for easier data wrangling
  names(results.kegg)[grep("EXACT_MASS", names(results.kegg))] = "weight"
  
  # add columns in order to combine HMDB and KEGG searches
  results.hmdb$PATHWAY = results.hmdb$BRITE = ""
  
  # combine all the HMDB and KEGG results together
  features <- c("Protein", 'entrez_id', 'symbol', "kegg_id","name",'m.z', 'omics_type', "experiment_id", "weight", 'PATHWAY', 'BRITE')
  results <- rbind(results.kegg[,features], results.hmdb[,features] )
  
  cat("Writing results out in long form...\n")
  write.table(results, gsub(".txt", "_ALL_RESULTS_long.txt", out_file), quote=F, row.names=F, sep='\t')
  
  # Aggregate data  c("Protein", 'entrez_id', 'symbol', "kegg_id","name",'m.z', 'omics_type', "experiment_id", "weight", 'PATHWAY', 'BRITE')
  #~~~~~~~~~~~~~~~~~
  cat("Aggregating
      \tKegg ID
      \tSymbol
      \tOmics Type
      \tExperiment ID
      \tMass
      \tPathway
      \tBrite ID\n")
  results.agg <- aggregate(data=results[,c("Protein", 'm.z', "weight", 'entrez_id', 'symbol', "kegg_id","name", 'omics_type', "experiment_id", 'PATHWAY', 'BRITE')], .~Protein+m.z+weight,  function(y){ paste(unique(y) ,collapse=";") } )
  
  # write out the results
  cat("Writing aggregated results...\n")
  write.table(results.agg, gsub(".txt", "_ALL_RESULTS_Aggregated.txt", out_file), quote=F, row.names=F, sep='\t')
  
  return(unique(results.agg))
}