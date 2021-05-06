

# Aggregate the results from the KEGG and HMDB search together
aggregateResults <- function(results.kegg, results.hmdb, out_file){
  cat("Combining HMDB and KEGG serch results...\n")
  
  keggMissing <- (nrow(results.kegg) == 0)
  hmdbMissing <- (nrow(results.hmdb) == 0)

  if (keggMissing & hmdbMissing){
    return (data.frame())
  }

  features <- c("Protein", 'entrez_id', 'symbol', "kegg_id","name",'m.z', 'omics_type', "experiment_id", "weight", 'PATHWAY', 'BRITE', 'matchSource', 'contrast', 'log2FC', 'adj.pvalue')

    # fix some of the names for easier data wrangling
  # and add columns so we can track source of match in long.txt
  if (!keggMissing){
    names(results.kegg)[grep("EXACT_MASS", names(results.kegg))] = "weight"
    results.kegg$matchSource <- "kegg"
    results.kegg <- results.kegg[,features]
  }
  if (!hmdbMissing){
    results.hmdb$PATHWAY = results.hmdb$BRITE = ""
    results.hmdb$matchSource <- "hmdb"
    results.hmdb <- results.hmdb[,features]
  }
  
  # combine all the HMDB and KEGG results together
  results <- rbind(results.kegg, results.hmdb)
  
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