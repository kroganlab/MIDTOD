## MIDTOD
R code for MIDTOD

## Dependencies
- R (>3.6.0)  

## Usage
the MIDTOD package will identify metabolites. For full options see

```bash
Rscript metabolomics_pipeline.R --help
```

Minimum command:
```bash
Rscript metabolomics_pipeline.R -f (results_file) -e (evidence_file)

arguments:
f=results [f]ile
e=[e]vidence file
o=[o]utput directory
s=[s]pecies
  accepted values: human, mouse
x=orthogonal data file
```
  
Launch GUI
```bash
Rscript -e "shiny::runApp(launch.browser = TRUE)"
```


For full control, run within R:

```R
source ("parent_path_of_midtod/MIDTOD/midtod.R", chdir=TRUE)

outDir <- "midtod/output/"
dir.create (outDir, recursive = TRUE, showWarnings = FALSE)

midtod(resultsFile  = "H1HTBERerun_06302019_results-wide.txt",
       evidenceFile = "H1HTBERerun_06302019-evidence.txt",
       outputDir    = outDir,
       species      = "human",
       remove.infinites = FALSE,
       orthogonalDataFile = "data/2020_09_03foldChangesForMIODTOD.htbe.thp1.csv.gz",
       log2FC = 1,
       pvalue = 0.05,
       filterResult = FALSE
       )

```