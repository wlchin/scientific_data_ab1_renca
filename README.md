# Introduction 

This snakemake repository provides the codebase to reproduce the analyses in the manuscript:

**Time-course bulk RNAseq and single-cell RNAseq of murine AB1 mesothelioma and Renca renal cancer following immune checkpoint therapy**

These workflows further describe the datasets used in [this repository](https://github.com/wlchin/IFNsignalling).  

# Notes

Raw data preprocessing (alignment and QC) is not performed in this workflow, although the code to reproduce the alignement steps is available as individual subworkflows. [MutliQC html reports](https://multiqc.info/) describing QC metrics are provided for bulk RNA-seq data and [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) per-sample diagnostics for single cell data are provided in this repo. 

The workflow leverages preprocessed data downstream of these two steps, downloading intermediate files from Cloudstor. 

To reproduce these workflows from scratch (raw sequencing files), please download raw files from the GEO repository [GSE153943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153942).

The workflows run on a linux machine. Where required, steps in the workflow are containerized. Therefore, invoke the --run-singularity flag when running snakemake. 

Label transfer is used to identify tumour cells. For the Renca (renal cell cancer dataset), the label transfer step uses external single cell data from Park et al. (2018) [^1]

# References

[^1]: Park, J., Shrestha, R., Qiu, C., Kondo, A., Huang, S., Werth, M., Li, M., Barasch, J., & Suszták, K. (2018). Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease. Science (New York, N.Y.), 360(6390), 758–763. https://doi.org/10.1126/science.aar2131

