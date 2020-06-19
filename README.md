# Versatile population genetic analyses using more than 10,000 individuals world-wide from 12 data sets

**Author:** Inken Wohlers  

**Summary:** This is a population genetic workflow that allows to perform analyses by selecting individuals (from n>10,000) from >100 populations compiled from 12 data sets. See the instructions in folder analysis_config for details.  
It has been used to generate population genetics results presented in the Egyptref paper.
Overall, individuals used for analyses can be selected from 12 public data sets.  
In order to run the analyses, the genotype data of the public data sets used has to be obtained from the relevant sources, which are listed in Supplementary Table 12 of the Egyptref Paper. 

**Instructions:**  
Raw base files need to be in place and Bioconda installed (https://bioconda.github.io).  
Please note: Some analyses need a considerable amount of memory (up to 180 GB RAM) and/or CPU (up to 48 cores).  
Run (parts of) the workflow to generate a target file via  

```
snakemake --use-conda TARGETFILENAME  
```  
