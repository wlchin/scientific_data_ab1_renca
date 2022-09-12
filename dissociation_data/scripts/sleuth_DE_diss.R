
library(sleuth)

path <- list.files("data/kallisto", full.names = T)
sample <- list.files("data/kallisto", full.names = F)
x <- strsplit(sample, "_")
sample_name <- unlist(lapply(x, `[[`, 3))
pheno <- substr(sample_name,1,nchar(sample_name)-2)
metadata <- data.frame(sample, sample_name, pheno, path, stringsAsFactors = F)

tot <- readRDS("resources/gene_to_transcript_map_gencode.rds")
design <- ~ pheno
so <- sleuth_prep(metadata, 
                 full_model = design, 
                 target_mapping = tot, 
                 read_bootstrap_tpm = TRUE,
                 extra_bootstrap_summary = TRUE,
                 #aggregation_column = "gene_name",
                 transformation_function = function(x) log2(x + 0.5))

so <- sleuth_fit(so)
oe <- sleuth_wt(so, which_beta = "phenoSingleCell")
sleuth_results_oe <- sleuth_results(oe, test = "phenoSingleCell", show_all = TRUE)

filt <- sleuth_results_oe[sleuth_results_oe$qval < 0.05 & sleuth_results_oe$b > 1,]
write(unique(filt$gene_name), file = snakemake@output[[1]])

filt <- sleuth_results_oe[sleuth_results_oe$qval < 0.05 & sleuth_results_oe$b < -1,]
write(unique(filt$gene_name), file = snakemake@output[[2]])

plot_pca(so, units = "tpm", color_by = "pheno")
ggplot2::ggsave(snakemake@output[[3]])

saveRDS(sleuth_results_oe, file = snakemake@output[[4]])
saveRDS(so, file = snakemake@output[[5]])
saveRDS(metadata, file = snakemake@output[[6]])


sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
saveRDS(sleuth_matrix, file = snakemake@output[[7]])