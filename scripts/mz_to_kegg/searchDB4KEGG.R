# Read in a file with a column of metabolite peaks and m/z, search KEGG for similar metabolites based on their weight, adduct, and +/- some threshold.
#   Then compare the associated entrez id's to the metabolites with the significant genes found in the other experiments (RNA seq, PH, UB, etc)
#   Return any significant genes and experiments found to match the metabolites.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  LOAD FUNCTIONS INTO MEMORY FIRST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

searchKEGG <- function(x, kegg, THRESH){
  
  # Debugging
  # x <- dat
  
  
  cat('>> Using fuzzy search to find matching masses in KEGG...\n')
  H= 1.007276
  NH4=18.033823
  Na=22.989218
  K=38.963158
  Cl=34.969402
  
  metabolites <- unique(x$m.z)
  metabolite_matches = list()
  pb <- txtProgressBar(min=0, max=length(metabolites), style=3)
  for(i in 1:length(metabolites)){
    # 
    # # Debugging
    # i <- 1
    # 
    mHp = metabolites[i] - H
    mNH4p = metabolites[i] - NH4
    mNap = metabolites[i] - Na
    mKp = metabolites[i] - K
    #mHn = metabolites[i] + H
    #mCln = metabolites[i] + Cl
    
    mHp = which( (kegg$EXACT_MASS >= (mHp-THRESH)) & (kegg$EXACT_MASS <= (mHp+THRESH)) )
    mNH4p = which( (kegg$EXACT_MASS >= (mNH4p-THRESH)) & (kegg$EXACT_MASS <= (mNH4p+THRESH)) )
    mNap = which( (kegg$EXACT_MASS >= (mNap-THRESH)) & (kegg$EXACT_MASS <= (mNap+THRESH)) )
    mKp = which( (kegg$EXACT_MASS >= (mKp-THRESH)) & (kegg$EXACT_MASS <= (mKp+THRESH)) )
    #mHn = which( (kegg$EXACT_MASS >= (mHn-THRESH)) & (kegg$EXACT_MASS <= (mHn+THRESH)) )
    #mCln = which( (kegg$EXACT_MASS >= (mCln-THRESH)) & (kegg$EXACT_MASS <= (mCln+THRESH)) )
    
    idx = c(mHp, mNH4p, mNap, mKp) #, mHn, mCln)
    if( sum(idx)>0 ){
      idx = unique(idx)
      metabolite_matches[[i]] = cbind(kegg[idx,], "m.z"=metabolites[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  message ("Of ", length(metabolites), " metabolites, ", length(metabolite_matches), " matches were found in KEGG")
  tmp = do.call(rbind,metabolite_matches)
  return(tmp)
}


#  Main function to handle search
searchDB4KEGG <- function(input_file, kegg_file, flu_file, out_file, THRESH, max_weight=Inf){
  
  # # Debugging: 
  # input_file <- to_search
  # max_weight=Inf
  # 
  
  cat(">> Loading files...\n")
  # Load metabolite list
  dat <- input_file
  # Load hmdb
  kegg <- read.delim(kegg_file, sep='\t', stringsAsFactors=F, header=T)
  
  # remove metabolites over the "max-weight" due to them probably not being useful (lipids, etc)
  if(max_weight != Inf){
    kegg <- kegg[!kegg$EXACT_MASS > max_weight,]
  }
  
  # Load significant hits from other OMICS datasets
  if ("data.frame" %in% class(flu_file)){
    flu = flu_file
  }else{
    flu <- read.delim(flu_file, sep='\t', stringsAsFactors=F, header=T)
    flu <- flu[,c('experiment_id','omics_type','condition_2','cell_line','strain','entrez_id','symbol')]   # remove unnecessary variables
  }
  # search for masses in KEGG 
  kegg_hits = searchKEGG(dat, kegg, THRESH)
 
  cat(">> Matching Peaks to significant hits in other data sets...\n")
  # match up all Peak names with the hits
  if (is.null(kegg_hits)){
    return (data.frame())
  }
  kegg_hits = merge(kegg_hits, dat, by='m.z')
  
  # get entrez id's corresponding to HMDB id's
  # kegg_hits = merge(kegg_hits, kegg.entrez[,c('hmdb', 'entrez_id')], by.x='id', by.y='hmdb')    ##### for MOUSE
  
  # check the rest of the flu datasets for significant matches
  kegg_hits.long = merge(kegg_hits, flu, by='entrez_id')
  
  # Remove duplicate metabolite ID's
  kegg_hits.long <- unique(kegg_hits.long[,-grep("condition_2|cell_line|strain", names(kegg_hits.long))])
  
  # convert symbols to uppercase  -- why?  
  #kegg_hits.long$symbol = toupper(kegg_hits.long$symbol)
  
  
  # aggregate all the gene names and experiment_id's
  cat(">> Adding gene names and experiment ID's...\n")
  tmp = aggregate(data=kegg_hits.long[c('kegg_id','name','m.z','Protein','EXACT_MASS','symbol', 'PATHWAY', 'BRITE')], symbol~., function(y){ paste(unique(y) ,collapse=",") } )
  tmp2 = aggregate(data=kegg_hits.long[c('kegg_id','name','m.z','Protein','EXACT_MASS','omics_type', 'PATHWAY', 'BRITE')], omics_type~., function(y){ paste(unique(y) ,collapse=",") } )
  # combine the data with the aggregate omics type and symbol
  tmp = merge(tmp, tmp2, by=c('kegg_id','name','m.z','Protein','EXACT_MASS', 'PATHWAY', 'BRITE'))
  
  # map on the p-values and fold changes
  kegg_hits.agg <- merge(tmp, input_file, by=c('Protein','m.z'))
  
  cat(">> Writing results to file...\n")
  write.table(unique(kegg_hits[, c("m.z", "Protein", "name","contrast", "log2FC", "adj.pvalue")]), gsub (".txt", "_KEGG_mass_matches.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  write.table(kegg_hits.agg ,gsub(".txt", "_KEGG_aggregated.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  write.table(kegg_hits.long ,gsub(".txt", "_KEGG_long.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  
  cat(">> Search Complete!!\n")
  return(kegg_hits.long)
}

#~~~~~~~~~~~~~~~~~~~~~~~~  END FUNCTIONS SECTION  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









