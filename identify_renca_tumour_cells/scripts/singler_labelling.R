
library(Seurat)
library(SingleR)
library(cowplot)

x <- readRDS(snakemake@input[[1]])  # load the seurat object
ref_mouse <- readRDS("data/mouse.rds")

DefaultAssay(x) <- "RNA" # from integrated
cluster_ids <- x$seurat_clusters  # Idents(x) singleR requires clusters
raw_counts <- x@assays$RNA@counts #documentation says that raw counts acceptable - get from seurat object

singler_mouse = SingleR(method = "cluster",
                        sc_data = raw_counts,
                        ref_data = ref_mouse$data,
                        types = ref_mouse$main_types,
                        clusters = cluster_ids,
  genes = "de", quantile.use = 0.8, p.threshold = 0.05,
  fine.tune = TRUE, fine.tune.thres = 0.05, sd.thres = 1,
  do.pvals = T, numCores = 2)

#x[["cluster_input"]] <- cluster_ids
#Idents(x) <- "cluster_input"

new.cluster.ids <- singler_mouse$labels1
names(new.cluster.ids) <- levels(x)
x <- RenameIdents(x, new.cluster.ids)
singlerlabels <- as.character(Idents(x))
saveRDS(singlerlabels, file = snakemake@output[[1]])

p1 <- DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
p2 <- DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot_grid(p1, p2)
ggplot2::ggsave("results/renca_labels.jpeg", height = 7, width = 14)






