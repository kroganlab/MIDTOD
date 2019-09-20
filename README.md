## MIDTOD
R code for MIDTOD

## Dependencies
- R (>3.6.0)  

## Usage
the MIDTOD package will identify metabolites  
```bash
Rscript metabolomics_pipeline.R -f (results_file) -e (evidence_file)

arguments:
f=results [f]ile
e=[e]vidence file
o=[o]utput directory
s=[s]pecies
  accepted values: human, mouse
```
  
Launch GUI
```bash
Rscript -e "shiny::runApp()"
```
