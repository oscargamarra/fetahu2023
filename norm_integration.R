# 1. Load packages ----
library(Seurat)        # Core single-cell RNA-seq analysis
library(patchwork)     # Combine multiple plots
library(dplyr)         # Data manipulation
library(ggplot2)       # Plotting
library(viridis)       # Color scales for visualization
library(org.Hs.eg.db)  # Gene annotation database

# 2. Global configuration ----
set.seed(1000)
CFG <- list()
CFG$work_dir <- "/Users/oscargamarrasiapo/Desktop/Fetahu_2023/analysis"
CFG$ndims <- 15
CFG$gradient.colors <- viridis(101, direction = -1)

# 3. Load post-QC Seurat objects ----
files <- list.files(CFG$work_dir, pattern = "_1postQC\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No *_1postQC.rds files found in work_dir.")

slist <- lapply(files, function(f) {
  message("Reading: ", f)
  so <- readRDS(f)
  samp <- tools::file_path_sans_ext(basename(f))
  
  # Debugging: Print sample name and number of cells
  message("Sample: ", samp)
  message("Number of cells: ", ncol(so))
  message("First few cell barcodes: ", paste(head(rownames(so@meta.data)), collapse = ", "))
  
  # Create metadata vector for sample
  sample_metadata <- rep(samp, ncol(so))
  names(sample_metadata) <- rownames(so@meta.data)  # Explicitly set names to match cell barcodes
  message("Sample metadata vector (first few): ", paste(head(sample_metadata), collapse = ", "))
  
  # Add sample metadata
  so <- AddMetaData(so, metadata = sample_metadata, col.name = "sample")
  
  # Map libraries to patients
  patient_map <- c("M1" = "P1", "M2a" = "P2", "M2b" = "P2", "M3" = "P3", "M4" = "P4")
  lib_id <- sub("_1postQC$", "", samp)
  
  # Debugging: Check library ID and patient mapping
  message("Library ID: ", lib_id)
  if (!lib_id %in% names(patient_map)) {
    stop(sprintf("Library ID '%s' from file '%s' not found in patient_map", lib_id, f))
  }
  
  # Create metadata vector for patient
  patient_metadata <- rep(as.character(patient_map[lib_id]), ncol(so))  # Convert to plain character
  names(patient_metadata) <- rownames(so@meta.data)  # Explicitly set names to match cell barcodes
  message("Patient metadata vector (first few): ", paste(head(patient_metadata), collapse = ", "))
  
  # Add patient metadata
  so <- AddMetaData(so, metadata = patient_metadata, col.name = "patient")
  
  # Debugging: Inspect metadata
  message("Metadata preview: ")
  print(head(so@meta.data))
  
  return(so)
})
names(slist) <- sapply(slist, function(x) unique(x$patient))

# 4. Per-sample SCTransform normalization ----
slist <- lapply(slist, function(so) {
  DefaultAssay(so) <- "RNA"
  # SCTransform: Normalizes data using regularized negative binomial regression
  # vst.flavor="v2" uses the v2 implementation; variable.features.n=3000 selects top 3000 variable genes
  so <- SCTransform(so, vst.flavor = "v2", verbose = FALSE, variable.features.n = 3000)
  return(so)
})

# 5. Merge Seurat objects ----
merged_preint <- merge(slist[[1]], slist[2:length(slist)], merge.data = TRUE)
merged_preint$lib <- as.factor(merged_preint$lib)
remove(slist)  # Clean up memory

# 6. Pre-integration analysis ----
DefaultAssay(merged_preint) <- "SCT"
merged_preint <- FindVariableFeatures(merged_preint, nfeatures = 3000)
merged_preint <- RunPCA(merged_preint, npcs = 50, verbose = FALSE)

# Visualize PCA and UMAP before integration
plots <- list()
plots$preint_pca <- DimPlot(merged_preint, reduction = "pca", pt.size = 1) +
  ggtitle("PCA before integration")
plots$preint_elbow <- ElbowPlot(merged_preint, reduction = "pca", ndims = 50) +
  ggtitle("Elbow plot for PCA")
plots$preint_heatmap <- DimHeatmap(merged_preint, dims = 1:18, cells = 100)

# Run UMAP on pre-integrated data
merged_preint <- RunUMAP(merged_preint, dims = 1:CFG$ndims)
plots$preint_umap <- DimPlot(merged_preint, reduction = "umap", group.by = "lib") +
  ggtitle("UMAP before integration")
saveRDS(merged_preint, file.path(CFG$work_dir, "merged_preint.rds"))

# 7. Integration with CCA ----
integrated <- IntegrateLayers(
  object = merged_preint,
  method = CCAIntegration,  # Canonical Correlation Analysis for integration
  normalization.method = "SCT",
  verbose = FALSE
)

# 8. Post-integration analysis ----
integrated <- FindNeighbors(integrated, reduction = "integrated.dr", dims = 1:CFG$ndims)
integrated <- FindClusters(integrated, resolution = 0.6)
integrated <- RunUMAP(integrated, reduction = "integrated.dr", dims = 1:CFG$ndims)
plots$int_umap <- DimPlot(integrated, reduction = "umap", group.by = "lib") +
  ggtitle("UMAP after integration")
saveRDS(integrated, file.path(CFG$work_dir, "integrated.rds"))

# 9. Cell cycle scoring ----
data(refdata_cellranger_GRCh38_3.0.0)
data(cc.genes)
integrated <- CellCycleScoring(
  object = integrated,
  s.features = cc.genes$s.genes,    # S-phase genes
  g2m.features = cc.genes$g2m.genes, # G2/M-phase genes
  assay = "RNA"
)
plots$cell_cycle <- DimPlot(integrated, reduction = "umap", pt.size = 0.5, group.by = "Phase") +
  ggtitle("Cell cycle phase")

# 10. Identify confounding genes ----
# Calculate correlations of genes with cell cycle scores
Scor <- metadataCorrelations(integrated, "S.Score")  # Correlation with S-phase score
G2Mcor <- metadataCorrelations(integrated, "G2M.Score")  # Correlation with G2/M-phase score

# Identify additional cell cycle-related genes
# derivedCellcycleGenes2: Identifies non-canonical genes correlated with cell cycle
additional <- derivedCellcycleGenes2(
  Scor = Scor,
  Sgenes = genelists$s.genes,
  G2Mcor = G2Mcor,
  G2Mgenes = genelists$g2m.genes
)

# Add module scores for confounding factors
integrated <- AddModuleScore(
  integrated,
  features = list(
    stress = genelists$stress,
    hb = genelists$hemo,
    fem = genelists$female,
    male = genelists$male
  ),
  name = c("stress", "hb", "fem", "male"),
  assay = "RNA"
)

# Visualize feature plots for QC metrics and stress score
plots$features <- FeaturePlot(
  integrated,
  pt.size = 0.5,
  features = c("nCount_RNA", "percent_mito", "nFeature_RNA", "stress1"),
  order = TRUE
)

# Compile genes to remove (cell cycle, stress, hemoglobin)
remove <- unique(c(
  cc.genes$s.genes,
  cc.genes$g2m.genes,
  additional$S.derived,
  additional$G2M.derived,
  genelists$stress,
  genelists$hemo
))
cat("Number of variable features to be removed:", length(intersect(VariableFeatures(integrated), remove)), "\n")

saveRDS(integrated, file.path(CFG$work_dir, "integrated_predeconfounding.rds"))
saveRDS(remove, file.path(CFG$work_dir, "remove_genes.rds"))

# 11. Re-run integration without confounding genes ----
merged_preint <- readRDS(file.path(CFG$work_dir, "merged_preint.rds"))
features_clean <- setdiff(rownames(merged_preint), remove)

integrated <- IntegrateLayers(
  object = merged_preint,
  method = CCAIntegration,
  normalization.method = "SCT",
  features = features_clean,
  verbose = FALSE
)
saveRDS(integrated, file.path(CFG$work_dir, "integrated_postdeconfounding.rds"))

# 12. Final dimensionality reduction and clustering ----
integrated <- RunPCA(integrated, assay = "SCT", npcs = 50, verbose = FALSE)
plots$final_pca <- DimPlot(integrated, reduction = "pca", pt.size = 1) +
  ggtitle("PCA after deconfounding")
plots$final_elbow <- ElbowPlot(integrated, reduction = "pca", ndims = 50) +
  ggtitle("Elbow plot after deconfounding")
plots$final_heatmap <- DimHeatmap(integrated, dims = 1:18, cells = 100) +
  ggtitle("Heatmap of top 18 PCs after deconfounding")

integrated <- FindNeighbors(integrated, reduction = "integrated.dr", dims = 1:CFG$ndims)
integrated <- FindClusters(integrated, resolution = 0.8)
integrated <- RunUMAP(integrated, reduction = "integrated.dr", dims = 1:CFG$ndims)

# Visualize final UMAPs
plots$final_umap_lib <- DimPlot(integrated, reduction = "umap", group.by = "lib") +
  ggtitle("UMAP after deconfounding (by library)")
plots$final_umap_phase <- DimPlot(integrated, reduction = "umap", group.by = "Phase") +
  ggtitle("UMAP after deconfounding (by cell cycle phase)")

# Save final results
saveRDS(integrated, file.path(CFG$work_dir, "integrated_filtered.rds"))

# 13. Display plots ----
print(plots$preint_pca)
print(plots$preint_elbow)
print(plots$preint_heatmap)
print(plots$preint_umap)
print(plots$int_umap)
print(plots$cell_cycle)
print(plots$features)
print(plots$final_pca)
print(plots$final_elbow)
print(plots$final_heatmap)
print(plots$final_umap_lib)
print(plots$final_umap_phase)

# Clean up memory
gc()