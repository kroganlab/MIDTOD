#!/usr/bin/env Rscript
# convert peak aligned metabolomic data from Markview RT into MSStats compatible format

suppressMessages(library(reshape2))
suppressMessages(library(data.table))
suppressMessages(library(optparse))


main <- function(data_file, id_file=NULL, out_file){ 
  cat(paste('> Reading in data from: ', data_file, '...\n', sep=""))
  dat=fread(data_file)
  tmp = melt(data=dat, id=c(1:7), variable.name="RawFile", value.name="Intensity")
  
  tmp[,'Row'] = NULL
  tmp[,'Index'] = NULL
  setnames(tmp, c(1,3), c('Sequence', 'Retention time'))
  tmp$Sequence <- gsub(' \\([0-9]+\\)','',tmp$Sequence)
  tmp$Proteins <- tmp$Sequence
  tmp$Charge = 1
  
  tmp = as.data.frame(tmp, stringsAsFactors=F)
  
  # annotate with  known id's
  if(!is.null(id_file)){
    ids = read.delim(id_file, stringsAsFactors=F, sep='\t')
    ids$KEGG = gsub('\\\xca','',ids$KEGG)
    for( i in 1:length(ids$KEGG)){  #### This could be optimized!!!!!
      idx = tmp$'m/z' %in% ids$m.z[i]
      tmp$Proteins[idx] = ids$KEGG[i]
    }
  }
  
  cat(paste('> Writing out data to: ', out_file, '...\n', sep=""))
  write.table(as.data.frame(tmp, stringsAsFactors=F), out_file, row.names=F, col.names=T, sep='\t', quote=F)
  cat('CONVERSION COMPLETE!\n\n')
}


option_list <- list(
  make_option(c("-f", "--data_file"), 
              help="data file containing Prospector output"),
  make_option(c("-i", "--id_file"), 
              help="file containing id's to swap out"),
  make_option(c("-o", "--out_file"), 
              help="output file for summary matrix")        
)
parsedArgs = parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=T))
main(parsedArgs$data_file, parsedArgs$id_file, parsedArgs$out_file)


# data_file='~/Box Sync/Projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Markerview_data/Markerview RT alignment_H1N1 lung samples 05152015.txt'
# id_file = '~/Box Sync/Projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite IDs/KEGG metabolites_HMDB.txt'
# out_file='~/projects/Fluomics/data/20150527_H1N1_lung_evidence.txt'
# main(data_file, id_file, out_file)
