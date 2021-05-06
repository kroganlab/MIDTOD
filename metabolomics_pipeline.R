#! /usr/bin/env Rscript
source(file = "midtod.R")
library(package = "optparse")

optionList <- list(make_option(opt_str = c("-f", "--file"),
			       type    = "character",
			       default = NULL,
			       help    = "result file in wide format",
			       metavar = "character"),
		   make_option(opt_str = c("-e", "--evidence"),
			       type    = "character",
			       default = NULL,
			       help    = "evidence file",
			       metavar = "character"),
		   make_option(opt_str = c("-o", "--output-dir"),
			       type    = "character",
			       default = getwd(),
			       help    = "output directory",
			       metavar = "character"),
		   make_option(opt_str = c("-s", "--species"),
			       type    = "character",
			       default = "human",
			       help    = "species",
			       metavar = "character"),
		   make_option(opt_str = c("-I", "--remove-infinites"),
		               type    = "logical",
		               default = FALSE,
		               help    = "remove infinite values in log2FC column (default is to leave them in)",
		               action="store_true",),
		   make_option(opt_str = c("-x", "--orthogonal-data-file"),
		               type    = "character",
		               default = NULL,
		               help    = "file of other OMICS results to find significantly regulated metabolite-associated-genes",
		               metavar = "character"),
		   make_option(opt_str = c("-m", "--mode"),
		               type = "character",
		               default = "positive",
		               help ="look for matching masses as if ions have a 'positive' (default) or 'negative' charge")
		   
)
optionParser <- OptionParser(option_list = optionList)

option <- parse_args(optionParser)

if (is.null(option$file)) {
  print_help(optionParser)
  stop("At least two arguments must be supplied (result and evidence files)",
       call. = FALSE)
}

midtod(resultsFile  = option$file,
       evidenceFile = option$evidence,
       outputDir    = option$"output-dir",
       species      = option$species,
       remove.infinites = option$"remove-infinites",
       orthogonalDataFile = option$"orthogonal-data-file",
       mode         = option$mode)
