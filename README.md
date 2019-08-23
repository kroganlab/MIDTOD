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
Evidence files can be created using the R package artMS using as input MarkerView peak 
table files
```{r}
> library(package = "artMS")
> artmsConvertMetabolomics(input_file = "fluomics.markerview.output.txt",
                           out_file   = "fluomics.artms.evidence.txt")
```

```bash
> head fluomics.markerview.output.txt
Row Index Peak Name    m/z     Ret. Time Group Use  Sample016
1   1	  70.0599_9.56 70.0599 9.56            True 1.0559
2   2	  72.0739_9.93 72.0739 9.93            True 0.11863
3   3	  86.0879_9.56 86.0879 9.56            True 1.6989

> head fluomics.artms.evidence.txt
Modified.sequence m/z     Retention time Group Use  RawFile   Intensity Proteins     Charge
70.0599_9.56      70.0599 9.56                 TRUE Sample016 1.0559    70.0599_9.56 1
72.0739_9.93      72.0739 9.93                 TRUE Sample016 0.11863	72.0739_9.93 1
86.0879_9.56      86.0879 9.56                 TRUE Sample016 1.6989    86.0879_9.56 1
```

Result files can be created using the R function artMS::artmsQuantification
```{r}
library(package = "artMS")
# create copy of template
artmsWriteConfigYamlFile(config_file_name = "metab_config.yaml")
# edit template
artmsQuantification(yaml_config_file = "fluomics.metab_config.yaml")
```

```bash
diff metab_config.yaml fluomics.metab_config.yaml
files:
<   evidence: /path/to/the/evidence.txt
<   keys: /path/to/the/keys.txt
<   contrasts: /path/to/the/contrast.txt
<   summary: /path/to/the/summary.txt
<   output: /path/to/the/output/results/results.txt
---
files:
>   evidence: fluomics.artms.evidence.txt
>   keys: fluomics.keys.txt
>   contrasts: fluomics.contrasts.txt
>   output: fluomics.artms.result.txt
8,9c7,8
qc:
<   basic: 1
<   extended: 1
---
qc:
>   basic: 0
>   extended: 0
18,19c17,18
  filters:
<     enabled: 1
<     contaminants: 1
---
  filters:
>     enabled: 0
>     contaminants: 0
35c34
output_extras:
<   enabled: 1
---
output_extras:
>   enabled: 0
37c36
  annotate:
<     enabled: 1
---
  annotate:
>     enabled: 0
```
