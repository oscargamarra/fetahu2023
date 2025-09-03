library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(dplyr)
library(ggplot2)
library(pheatmap)

set.seed(1000)
CFG <- list()
CFG$work_dir <- "/Users/oscargamarrasiapo/Desktop/Fetahu_2023/analysis"

integrated <- readRDS(file = file.path(CFG$work_dir,
                                       "integrated_cellsidentified.rds"))
remove <- readRDS(file = file.path(CFG$work_dir,
                                   "remove.rds"))

neurons <- subset(integrated, subset = singler == "Neurons")

neurons_list <- SplitObject(neurons, split.by = "lib")

neurons_list <- lapply(neurons_list, function(so){
  
  DefaultAssay(so) <- "RNA"
  # SCTransform (no regression at this stage)
  so <- SCTransform(so, vst.flavor = "v2", verbose = FALSE,
                    variable.features.n = 3000)
  return(so)
})

neurons_preint <- merge(neurons_list[[1]], neurons_list[2:length(neurons_list)], merge.data = TRUE)
neurons_preint$lib <- as.factor(neurons_preint$lib)

neurons_preint <- FindVariableFeatures(neurons_preint, nfeatures = 3000)
features <- VariableFeatures(neurons_preint)
# PCA
neurons_preint <- RunPCA(neurons_preint, npcs = 50, verbose = FALSE)

neurons_preint <- RunUMAP(neurons_preint, dims = 1:50)

features_clean <- setdiff(features, remove)


neurons_integrated <- IntegrateLayers(
  object = neurons_preint,            # list of Seurat objects
  method = CCAIntegration,
  normalization.method = "SCT",
  features = features_clean,         # use the cleaned feature list
  k.weight = 80,
  verbose = FALSE
)
remove(neurons_list)


neurons_cds <- as.cell_data_set(neurons_integrated)
rowData(neurons_cds)$gene_short_name <- rownames(neurons_cds)

neurons_cds <- preprocess_cds(
  neurons_cds,
  num_dim = 50,
  method = "PCA",
  use_genes = rownames(neurons_cds), 
  scaling = FALSE
)

reducedDims(neurons_cds)$PCA <- Embeddings(neurons_integrated, "pca")

neurons_cds <- cluster_cells(neurons_cds, reduction_method = "PCA")
neurons_cds <- cluster_cells(neurons_cds, reduction_method = "UMAP")

plot_cells(neurons_cds, show_trajectory_graph = FALSE)

neurons_cds <- learn_graph(neurons_cds, use_partition = FALSE)

plot_cells(neurons_cds,
           genes=c("NCAM1", "TOP2A", "PHOX2B", "MKI67"),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)

neurons_cds <- order_cells(neurons_cds)

plot_cells(neurons_cds, color_cells_by = "pseudotime")


dynamic_genes <- graph_test(neurons_cds, neighbor_graph = "principal_graph", cores = 2)

dynamic_genes %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>%
  head()

dynamic_genes_sig <- subset(dynamic_genes, q_value < 0.05)


deg_ids <- rownames(
  subset(
    dynamic_genes[order(dynamic_genes$morans_I, decreasing = TRUE), ],
    q_value < 0.05)
)
    

# 10. Plot heatmap of dynamic genes ordered by pseudotime

gene_modules <- find_gene_modules(neurons_cds,
                                  resolution=c(10^seq(-6,-1)))

table(gene_modules$module)

# Bin cells along pseudotime (e.g., 50 bins)
pseudotime_bins <- cut_number(neurons_cds@principal_graph_aux[["UMAP"]]$pseudotime, 15)

cell_groups <- data.frame(
  cell = row.names(colData(neurons_cds)),
  pseudotime_bin = pseudotime_bins
)

agg_mat <- aggregate_gene_expression(
  neurons_cds,
  gene_group_df = gene_modules,
  cell_group_df = cell_groups
)

row.names(agg_mat) <- paste0("Module ", row.names(agg_mat))

# Heatmap: modules (rows) vs pseudotime bins (columns)
pheatmap::pheatmap(
  agg_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  clustering_method = "ward.D2"
)

module_sizes <- table(gene_modules$module)
print(module_sizes)

#Extract top 5 genes per module
# Select top 5 genes by score for each module
top5_genes <- gene_modules %>%
  group_by(module) %>%
  arrange(desc(score)) %>%  # Sort by score in descending order
  slice_head(n = 5) %>%     # Take top 5 genes
  ungroup() %>%
  select(module, gene, score)  # Select relevant columns

# Display top 5 genes per module
message("Top 5 genes per module (by score):")
print(top5_genes)
