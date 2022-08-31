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
m=[m]ode: positive or negative.  The expected charge on the metabolites in MS
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
       orthogonalDataFile = "data/2020_09_03foldChangesForMIODTOD.htbe.thp1.tsv.gz",
       log2FC = 1,
       pvalue = 0.05,
       filterResult = FALSE,
       mode = "positive"
       )

```

Orthogonal data set is expected to have the following column names

```
  experiment_id
  omics_type
  condition_2
  cell_line
  strain
  entrez_id
  symbol
  q_value
  log2fc
```

Orthogonal data set example. Formatted here for ease of reading.  Actual inputted file should be tab separated text values.
```
              experiment_id omics_type condition_2 cell_line strain entrez_id  symbol      q_value    log2fc
    1: PH-site_htbe_H3N2_18    PH-site    H3N2_18H      htbe   H3N2     10528   NOP56 0.0041462891  2.798370
    2: PH-site_htbe_H3N2_12    PH-site    H3N2_12H      htbe   H3N2     26037 SIPA1L1 0.0291694566 -1.703825
    3: PH-site_htbe_H3N2_12    PH-site    H3N2_12H      htbe   H3N2      6449    SGTA 0.0157542912  1.960473
    4: PH-site_htbe_H3N2_12    PH-site    H3N2_12H      htbe   H3N2      1736    DKC1 0.0264957502  2.323113
    5: PH-site_htbe_H3N2_12    PH-site    H3N2_12H      htbe   H3N2     23122  CLASP2 0.0297024129  2.118161
  ---                                                                                                      
38280:   RNAseq_thp1_H5N1_6     RNAseq        H5N1      thp1   H5N1    168374   ZNF92 0.0019812232 -1.141245
38281:   RNAseq_thp1_H5N1_6     RNAseq        H5N1      thp1   H5N1      7589 ZSCAN21 0.0202302982 -1.804801
38282:   RNAseq_thp1_H5N1_6     RNAseq        H5N1      thp1   H5N1     90204  ZSWIM1 0.0060407318 -1.123682
38283:   RNAseq_thp1_H5N1_6     RNAseq        H5N1      thp1   H5N1    140831  ZSWIM3 0.0006518764 -1.451436
38284:   RNAseq_thp1_H5N1_6     RNAseq        H5N1      thp1   H5N1    158586    ZXDB 0.0400759167 -1.048902
```


