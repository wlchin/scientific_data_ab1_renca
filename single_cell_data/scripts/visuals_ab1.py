import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

ab1 = sc.read(snakemake.input[0])

sc.set_figure_params(fontsize = 15, dpi_save = 600)
ab1.obs["hums"] = ab1.obs["analysis_ident"].astype("string")
sc.pl.umap(ab1, color = ["hums"], legend_loc = None, title = "UMAP of all AB1 samples", save = "ab1.png")

sc.set_figure_params(fontsize = 20, figsize = (6,6), dpi_save = 600)
sc.pl.umap(ab1, color = ["nCount_RNA", "nFeature_RNA", "percent.mt"], 
           title = ["RNA counts per cell", "Number of RNA features per cell", "% Mitochondrial genes per cell"], wspace = 0.3, save = "ab1_cellcomp.png")



sc.set_figure_params(fontsize = 10, dpi_save = 600)
ab1.obs["hums"] = ab1.obs["seurat_clusters"].astype("string")
sc.pl.umap(ab1, color = ["hums"], legend_loc = "on data", title = "UMAP of all AB1 samples", save = "ab1cluster.png")


sc.set_figure_params(fontsize = 15, dpi_save = 600)
ab1.obs["hums"] = ab1.obs["sample"].astype("category")
sc.pl.umap(ab1, color = ["hums"], 
           title = "UMAP of all AB1 samples", 
           size=10, alpha = 0.5, palette = "Dark2", sort_order=False, save = "ab1sample.png")


df_ab1 = ab1.obs
tu = df_ab1.groupby(["sample", "analysis_ident"]).size()
df = tu.unstack(-1)
df.index = df.index.astype("string")
df.columns = df.columns.astype("string")
df = df.reset_index()

df_total = df.sum(axis=1)
df_rel = df[df.columns[1:]].div(df_total, 0)*100
df_rel["samples"] = df["sample"]
df_rel.columns = df_rel.columns.astype("category")

sc.set_figure_params(fontsize = 15, figsize = (10,6), dpi_save = 600)
ax = df_rel.plot(
    x = 'samples',
    kind = 'barh',
    stacked = True,
    mark_right = True)

ax.legend(bbox_to_anchor=(1.0, 1.0))
ax.spines["left"]


ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.set_xlabel('percentage', fontsize = "x-large")
ax.set_ylabel('samples', fontsize = "x-large", labelpad = 20)
ax.set_title('Proportion of Cell Types by Sample (AB1)', fontsize = "x-large")
ax.grid(False)
ax.tick_params(axis='y', which='major', labelsize=15)
plt.tight_layout()
plt.savefig("figures/ab1_cell_proportions.png")




