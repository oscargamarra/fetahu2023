library(Seurat)
library(ggplot2)
library(org.Hs.eg.db)
library(Polychrome)
library(SingleR)
library(celldex)

set.seed(1000)
CFG <- list()
CFG$work_dir <- "/Users/oscargamarrasiapo/Desktop/Fetahu_2023/analysis"

integrated <- readRDS(file = file.path(CFG$work_dir, "integrated_filtered.rds"))

integrated <- FindClusters(integrated, resolution = seq(0.4, 2, 0.2), algorithm = 1)

columns <- grep("SCT_snn_res", names(integrated@meta.data), value = TRUE)

## show numbers of cells per cluster (top row shows the cluster ids)
## See changes in number of clusters and the cells at each cluster
for (col in columns) {
  cat("\n====\n", col, ":")
  print(table(integrated@meta.data[[col]]))
}


p_clusbyreso <- list()  # set up the list to hold the intermediate plots

for (reso in c("0.8", "1", "1.2")) {
  
  clus.column <- paste0("SCT_snn_res.", reso)
  clusters <- sort(unique(integrated@meta.data[[clus.column]]))
  
  ## switching color scheme to make it look different:
  clus.colors <- Polychrome::palette36.colors(length(clusters))
  ## clus.colors <-
  ## Polychrome::kelly.colors(length(clusters))
  names(clus.colors) <- as.character(clusters)
  
  p <- DimPlot(integrated, pt.size = 0.3, group.by = clus.column, cols = clus.colors) +
    labs(title = reso)
  p_clusbyreso[[reso]] <- p
}

## use patchwork to plot them side-by-side:
p_clusbyreso[[1]] | p_clusbyreso[[2]] | p_clusbyreso[[3]]

#After seeing this, i think that .8 is adequate.
clus.column <- paste0("SCT_snn_res.", 0.8)
Idents(integrated) <- clus.column


#Now let's do cell typing
# ---- Cell identification with SingleR ----

## Load reference
hpca <- HumanPrimaryCellAtlasData()

## Run SingleR on your integrated object
singler_hpca <- SingleR(
  test = GetAssayData(integrated, assay = "RNA", slot = "data"),
  ref = hpca,
  labels = hpca$label.main
)

## Plot the heatmap
plotScoreHeatmap(singler_hpca, 
                 show.labels = TRUE, max.labels = 100,
                 show.pruned = FALSE, order.by = "clusters",
                 clusters = integrated@meta.data[[clus.column]])

## Add annotation to metadata
integrated@meta.data$singler <- singler_hpca$labels

# ---- Explore composition ----
cbind(table(integrated@meta.data$singler))

# ---- UMAP plots ----
# Assign colors
typenames <- unique(integrated@meta.data$singler)
singler.colors <- Polychrome::palette36.colors(length(typenames))
names(singler.colors) <- typenames

# UMAP colored by SingleR cell type
p_singler <- DimPlot(integrated, reduction = "umap", pt.size = 0.5,
                     group.by = "singler", cols = singler.colors)

# UMAP colored by sample ("lib")
p_sample <- DimPlot(integrated, reduction = "umap", pt.size = 0.5,
                    group.by = "lib")

# UMAP colored by patient
p_patient <- DimPlot(integrated, reduction = "umap", pt.size = 0.5,
                     group.by = "patient")

# Combine plots
p_singler | p_sample
p_singler | p_patient


# Build contingency table
tab <- table(integrated$singler, integrated$SCT_snn_res.0.8)

# Find predominant cell type per cluster
majority_celltype <- apply(tab, 2, function(x) {
  names(which.max(x))
})

# Check mapping: cluster → predominant cell type
majority_celltype

# Add to metadata
integrated$cluster_celltype <- plyr::mapvalues(
  x = Idents(integrated), 
  from = names(majority_celltype),
  to   = majority_celltype
)

#See new simplified annotation
DimPlot(
  integrated, 
  reduction = "umap", 
  pt.size = 0.5,
  group.by = "cluster_celltype"
)


# ---- Save ----
saveRDS(integrated, file = file.path(CFG$work_dir,
                                     "integrated_cellsidentified.rds"))
# integrated <- readRDS(file = file)

