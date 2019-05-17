
# this script is to convert the random formatting of the processed data to 
#   match the formatting in the original design.


fixFormatting <- function(dat_file){
  dat <- read.delim(dat_file, stringsAsFactors = F)
  
  # alter names to match db
  dat$experiment_id = paste0(dat$OMICS, "-", dat$Model)           # experiment
  names(dat)[grep("Experiment", names(dat))] <- "omics_type"          # omics_type
names(dat)[grep("Comparison", names(dat))] <- "condition_2"         # condition_2
  names(dat)[grep("Model", names(dat))] <- "cell_line"                # cell_line
  names(dat)[grep("FLU", names(dat))] <- "strain"                     # strain
  dat$status <- "FINAL"                                               # status
                     # id    <-- Need to use entrez id's

  names(dat)[grep("Gene", names(dat))] <- "symbol"                    # symbol
  if(length(grep("Uniprot_PTM", names(dat)))>0){
    dat$id <- paste(dat$Uniprot_PTM, dat$symbol, sep='|')
    names(dat)[grep("Uniprot_PTM", names(dat))] <- "uniprot_ac"       # uniprot_ac
    dat$Protein <- NULL
                    # id    <-- Need to use entrez id's  
  }else{
    dat$id <- paste(dat$Protein, dat$symbol, sep='|')                   # id    <-- Need to use entrez id's
    names(dat)[grep("Protein", names(dat))] <- "uniprot_ac"
  }
  dat$kegg_id <- ""                                                   # kegg_id
  names(dat)[grep("iLog2FC", names(dat))] <- "log2fc"                 # log2fc
  names(dat)[grep("iPvalue", names(dat))] <- "q_value"                # q_value
  
  if(length(grep("Significance", names(dat)))>0) dat[,grep("Significance", names(dat))] = NULL
  
  write.table(dat, gsub(".txt", "_newHeaders.txt", dat_file), quote=F, row.names=F, sep='\t')
  return(dat)
}






ptm.mouse = fixFormatting('~/Box Sync/projects/FluomicsModeling/results/proteomics2017/data.fluomics.PTMsites.mouse.txt')
ptm.htbe = fixFormatting('~/Box Sync/projects/FluomicsModeling/results/proteomics2017/data.fluomics.PTMsites.htbe.txt')

data.mouse = fixFormatting('~/Box Sync/projects/FluomicsModeling/results/proteomics2017/fluomics.de.v5-mouse-total.txt')
data.htbe = fixFormatting('~/Box Sync/projects/FluomicsModeling/results/proteomics2017/fluomics.de.v5-htbe-total.txt')

dat.all <- rbind(ptm.mouse, ptm.htbe, data.mouse, data.htbe)




# Add ENTREZ ID
#~~~~~~~~~~~~~~~
mappingUniprotEntrezHuman <- read.delim('~/Box Sync/db/uniprot/human-uniprot2entrez_20170104.txt', header = T, sep = "\t", stringsAsFactors = F)
mappingUniprotEntrezMouse <- read.delim('~/Box Sync/db/uniprot/mouse-uniprot2entrez_20170105.txt', header = T, sep = "\t", stringsAsFactors = F)

entrez_id <- rbind( mappingUniprotEntrezHuman, mappingUniprotEntrezMouse )

tmp <- merge(dat.all, entrez_id, by.x="uniprot_ac", by.y="Uniprot")
names(tmp)[grep("entrez", names(tmp))] <- "entrez_id"

write.table(tmp, '~/Box Sync/projects/FluomicsModeling/results/proteomics2017/data.fluomics.ptm.v5.txt', quote=F, row.names=F, sep='\t')




















