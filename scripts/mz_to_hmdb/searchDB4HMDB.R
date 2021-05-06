# Read in a file with a column of metabolite peaks and m/z, search HMDB for similar metabolites based on their weight, adduct, and +/- some threshold.
#   Then compare the associated entrez id's to the metabolites with the significant genes found in the other experiments (RNA seq, PH, UB, etc)
#   Return any significant genes and experiments found to match the metabolites.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  LOAD FUNCTIONS INTO MEMORY FIRST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


searchHMDB <- function(x, hmdb, THRESH, modMasses){
  cat('>> Using fuzzy search to find matching masses in HMDB...\n')

  metabolites <- unique(x$m.z)
  metabolite_matches = list()
  pb <- txtProgressBar(min=0, max=length(metabolites), style=3)
  for(i in 1:length(metabolites)){
    idx <- matchMZAgainstMetaboliteDB(queryMZ = metabolites[i], targetMasses=hmdb$weight,
                                      modMasses = modMasses, THRESH=THRESH)
    if( sum(idx)>0 ){
      idx = unique(idx)
      metabolite_matches[[i]] = cbind(hmdb[idx,], "m.z"=metabolites[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  message ("Of ", length(metabolites), " metabolites, ", length(metabolite_matches), " matches were found in DB")
  tmp = do.call(rbind,metabolite_matches)
  return(tmp)
}


#  Main function to search handle search
searchDB4HMDB <- function(input_file, hmdb_file, hmdb.entrez_file, flu_file, out_file, THRESH, max_weight=Inf, modMasses = modMasses.positive){
  cat(">> Loading files...\n")
  # Load metabolite list
  dat <- input_file
  # Load hmdb
  hmdb <- read.delim(hmdb_file, sep='\t', stringsAsFactors=F, header=T)
  
  # remove metabolites over the "max-weight" due to them probably not being useful (lipids, etc)
  if(max_weight != Inf){
    hmdb <- hmdb[!hmdb$weight > max_weight,]
  }
  
  # load hmdb file with corresponding entrez id's
  hmdb.entrez <- read.delim(hmdb.entrez_file, sep='\t', stringsAsFactors=F, header=T)
  names(hmdb.entrez)[grep('entrez_gene_id_mouse',names(hmdb.entrez))] = 'entrez_id'
  
  
  # Load significant hits from other OMICS datasets
  if ("data.frame" %in% class(flu_file)){
    flu = flu_file
  }else{
    flu <- read.delim(flu_file, sep='\t', stringsAsFactors=F, header=T)
    flu <- flu[,c('experiment_id','omics_type','condition_2','cell_line','strain','entrez_id','symbol')]   # remove unnecessary variables
  }
  
  # search for masses in HMDB 
  hmdb_hits = searchHMDB(dat, hmdb, THRESH, modMasses)

  if (is.null(hmdb_hits)){
    return (data.frame())
  }

  cat(">> Matching Peaks to significant hits in other data sets...\n")
  # match up all Peak names with the hits
  hmdb_hits = merge(hmdb_hits, dat, by='m.z')
  # get entrez id's corresponding to HMDB id's
  hmdb_hits = merge(hmdb_hits, hmdb.entrez[,c('hmdb', 'entrez_id')], by.x='id', by.y='hmdb')    ##### for MOUSE
  # check the rest of the flu datasets for significant matches
  hmdb_hits.long = merge(hmdb_hits, flu, by.x='entrez_id', by.y='entrez_id')
  
  # Remove duplicate metabolite ID's
  hmdb_hits.long <- unique(hmdb_hits.long[,-grep("condition_2|cell_line|strain", names(hmdb_hits.long))])
  
  # aggregate all the gene names and experiment_id's
  cat(">> Adding gene names and experiment ID's...\n")
  tmp = aggregate(data=hmdb_hits.long[c('id','m.z','Protein','weight','symbol')], symbol~., function(y){ paste(unique(y) ,collapse=",") } )
  tmp2 = aggregate(data=hmdb_hits.long[c('id','m.z','Protein','weight','omics_type')], omics_type~., function(y){ paste(unique(y) ,collapse=",") } )
  # combine the data with the aggregate omics type and symbol
  tmp = merge(tmp, tmp2, by=c('id','m.z','Protein','weight'))
  
  # map on the p-values and fold changes
  tmp <- merge(tmp, input_file, by=c('Protein','m.z'))
  
  cat(">> Writing results to file...\n")
  write.table(unique(hmdb_hits[, c("m.z", "Protein", "name","contrast", "log2FC", "adj.pvalue")]), gsub (".txt", "_HMDB_mass_matches.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  write.table(tmp ,gsub(".txt", "_HMDB_Aggregated.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  write.table(hmdb_hits.long ,gsub(".txt", "_HMDB_long.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  
  cat(">> Search Complete!!\n")
  return(hmdb_hits.long)
}

#~~~~~~~~~~~~~~~~~~~~~~~~  END FUNCTIONS SECTION  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









