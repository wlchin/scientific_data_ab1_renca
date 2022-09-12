# Introduction 

This snakemake repository provides the codebase to reproduce the analyses in the manuscript:

**Time-course bulk RNAseq and single-cell RNAseq of murine AB1 mesothelioma and Renca renal cancer following immune checkpoint therapy**

These workflows further describe the datasets used in [this repository](https://github.com/wlchin/IFNsignalling).  

# Notes

Raw data preprocessing (alignment and QC) is not performed in this workflow, although the code to reproduce these steps is available as individual subworkflows. 

Instead, the workflow leverages preprocessed data downstream of these two steps, downloading intermediate files from Cloudstor. 

To reproduce these workflows from scratch (raw sequencing files), please download raw from the GEO repository [GSE153943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153942).


