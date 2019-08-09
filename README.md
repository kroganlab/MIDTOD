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

## Note
Evidence files can be created using the R package artMS from MarkerView peak table
files
```{r}
> library(package = "artMS")
> artmsConvertMetabolomics(input_file = "markerview.output.txt",
                           out_file   = "markerview.evidence.txt")
```

```bash
> head markerview.output.txt
Row Index Peak Name    m/z     Ret. Time Group Use  Sample016
1   1	  70.0599_9.56 70.0599 9.56            True 1.0559
2   2	  72.0739_9.93 72.0739 9.93            True 0.11863
3   3	  86.0879_9.56 86.0879 9.56            True 1.6989

> head markerview.evidence.txt
Modified.sequence m/z     Retention time Group Use  RawFile   Intensity Proteins     Charge
70.0599_9.56      70.0599 9.56                 TRUE Sample016 1.0559    70.0599_9.56 1
72.0739_9.93      72.0739 9.93                 TRUE Sample016 0.11863	72.0739_9.93 1
86.0879_9.56      86.0879 9.56                 TRUE Sample016 1.6989    86.0879_9.56 1
```
