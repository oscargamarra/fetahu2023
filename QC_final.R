# 1. Load packages ----
library(Seurat)        # Core single-cell RNA-seq analysis
library(Matrix)        # Sparse matrix operations
library(ggplot2)       # Plotting
library(patchwork)     # Combine multiple plots
library(org.Hs.eg.db)  # Gene annotation database
library(scDblFinder)   # Doublet detection
library(viridis)       # Color scales for visualization

# 2. Global configuration ----
set.seed(1000)
CFG <- list()

# Paths
CFG$data_dir   <- "~/Desktop/Fetahu_2023/data/"
CFG$output_dir <- "~/Desktop/Fetahu_2023/analysis"

# Parameters
CFG$ndims       <- 15
CFG$random_seed <- 2033

# QC thresholds per library
CFG$thresholds <- list(
  M1  = list(min_txpts = 1000, max_txpts = 15000, max_pctmito = 10, max_pct_hemo = 2),
  M2a = list(min_txpts = 1000, max_txpts = 40000, max_pctmito = 20, max_pct_hemo = 2),
  M2b = list(min_txpts = 1000, max_txpts = 40000, max_pctmito = 10, max_pct_hemo = 2),
  M3  = list(min_txpts =  500, max_txpts = 30000, max_pctmito = 20, max_pct_hemo = 2),
  M4  = list(min_txpts =  500, max_txpts = 30000, max_pctmito = 20, max_pct_hemo = 2)
)

CFG$gradient.colors <- viridis(101, direction = -1)

# 3. Storage objects ----
srat_list <- list()  # Store Seurat objects for each library
qc_plots  <- list()  # Store QC plots (density, scatter, violin)

# 4. Load data and initial QC ----
filt_files <- list.files(path = CFG$data_dir,
                         pattern = "_filtered_feature_bc_matrix.h5$",
                         full.names = TRUE)

for (f in filt_files) {
  # Extract library name from file path
  lib <- basename(f)
  lib <- sub("_filtered_feature_bc_matrix.h5$", "", lib)
  lib <- sub("^GSM[0-9]+_", "", lib)
  
  thr <- CFG$thresholds[[lib]]
  
  # Load counts from 10X h5 file
  counts <- Read10X_h5(f)
  cnts_per_cell <- Matrix::colSums(counts)
  
  cat("\n###", lib, "###\n")
  print(summary(cnts_per_cell))
  
  # Density plot of counts per cell
  df <- data.frame(counts = cnts_per_cell)
  qc_plots[[lib]]$density <- ggplot(df, aes(x = counts)) +
    geom_density() + geom_rug() +
    scale_x_continuous(trans = "log2") +
    ggtitle(paste0("Density of counts per cell:"))
  
  # Rename cells with library prefix
  colnames(counts) <- paste0(lib, "_", colnames(counts))
  
  # Calculate mitochondrial fraction
  mitos <- grep("^mt-", rownames(counts), value = TRUE, ignore.case = TRUE)
  cat("Mitochondrial genes found:\n", paste(mitos, collapse = ", "), "\n")
  
  percent_mito <- 100 * Matrix::colSums(counts[mitos,]) / cnts_per_cell
  nuclear <- setdiff(rownames(counts), mitos)
  
  # Create metadata for QC
  meta <- data.frame(
    percent_mito = percent_mito,
    log2_counts  = log2(1 + Matrix::colSums(counts[nuclear,])),
    log2_features = log2(1 + Matrix::colSums(counts[nuclear,] > 0)),
    lib = rep(lib, ncol(counts))
  )
  
  # QC scatter plot: log2 counts vs features
  qc_plots[[lib]]$scatter <- ggplot(meta, aes(x = log2_counts, y = log2_features)) +
    geom_point(color = "blue", alpha = 0.5) +
    ggtitle(paste0("QC plot: ", lib))
  
  # Create Seurat object with nuclear genes only
  counts <- counts[nuclear,]
  srat <- CreateSeuratObject(counts = counts, project = lib, meta.data = meta)
  
  # QC plots
  # Violin plot of transcript counts
  v <- VlnPlot(srat, "nCount_RNA", group.by = "lib")
  qc_plots[[lib]]$vln_linear <- v + geom_hline(yintercept = seq(0, 1e5, 2000), col = "grey", lwd = 0.1)
  qc_plots[[lib]]$vln_log <- v + scale_y_continuous(trans = "log2") +
    geom_hline(yintercept = seq(0, 1e4, 250), col = "grey", lwd = 0.1)
  
  # Scatter plot: nCount vs nFeature
  qc_plots[[lib]]$scatter_linear <- FeatureScatter(srat, "nCount_RNA", "nFeature_RNA", pt.size = 0.5) +
    geom_vline(xintercept = c(thr$min_txpts, thr$max_txpts), linetype = 2)
  qc_plots[[lib]]$scatter_log <- qc_plots[[lib]]$scatter_linear +
    scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")
  
  # Violin plot of mitochondrial content
  qc_plots[[lib]]$mito_vln <- VlnPlot(srat, "percent_mito", group.by = "lib") +
    geom_hline(yintercept = seq(0, 100, 5), col = "grey", lwd = 0.1)
  
  # Scatter plot: nCount vs percent mito
  qc_plots[[lib]]$mito_scatter <- FeatureScatter(srat, "nCount_RNA", "percent_mito", pt.size = 0.5) +
    geom_hline(yintercept = thr$max_pctmito, linetype = 2) +
    geom_vline(xintercept = c(thr$min_txpts, thr$max_txpts), linetype = 2)
  qc_plots[[lib]]$mito_scatter_log <- qc_plots[[lib]]$mito_scatter +
    scale_x_continuous(trans = "log2")
  
  # Store Seurat object
  srat_list[[lib]] <- srat
}

# 5. Visual QC summaries ----
# Display all QC plots in grids
grid.arrange(grobs = lapply(qc_plots, `[[`, "density"), ncol = 3, nrow = 2,
             top = "Density of counts per cell")
grid.arrange(grobs = lapply(qc_plots, `[[`, "scatter"), ncol = 3, nrow = 2,
             top = "Counts vs Features (QC)")
grid.arrange(grobs = lapply(qc_plots, `[[`, "vln_linear"), ncol = 3, nrow = 2,
             top = "Transcript counts per cell (linear)")
grid.arrange(grobs = lapply(qc_plots, `[[`, "vln_log"), ncol = 3, nrow = 2,
             top = "Transcript counts per cell (log)")
grid.arrange(grobs = lapply(qc_plots, `[[`, "scatter_linear"), ncol = 3, nrow = 2,
             top = "nCount vs nFeature (linear)")
grid.arrange(grobs = lapply(qc_plots, `[[`, "scatter_log"), ncol = 3, nrow = 2,
             top = "nCount vs nFeature (log-log)")
grid.arrange(grobs = lapply(qc_plots, `[[`, "mito_vln"), ncol = 3, nrow = 2,
             top = "Mitochondrial content per cell")
grid.arrange(grobs = lapply(qc_plots, `[[`, "mito_scatter"), ncol = 3, nrow = 2,
             top = "nCount RNA vs percent mito (linear)")
grid.arrange(grobs = lapply(qc_plots, `[[`, "mito_scatter_log"), ncol = 3, nrow = 2,
             top = "nCount RNA vs percent mito (log-log)")

# 6. Hemoglobin annotation and diagnostic plots ----
data(refdata_cellranger_GRCh38_3.0.0)
hb_genes_all <- genelists$hemo
hemo_plots <- list()

for (lib in names(srat_list)) {
  srat <- srat_list[[lib]]
  thr <- CFG$thresholds[[lib]]
  
  # Identify hemoglobin genes present in the library
  hb_genes <- intersect(hb_genes_all, rownames(srat))
  if (length(hb_genes) == 0) {
    warning(paste("No hemoglobin genes found for", lib))
    next
  }
  
  # Calculate hemoglobin percentage and add to metadata
  hb_counts <- Matrix::colSums(srat@assays$RNA$counts[hb_genes,])
  srat <- AddMetaData(srat, col.name = "pct_hemo", metadata = 100 * hb_counts / srat@meta.data$nCount_RNA)
  srat_list[[lib]] <- srat
  
  # Violin plot of hemoglobin content
  hemo_plots[[lib]]$violin <- VlnPlot(srat, features = "pct_hemo", pt.size = 0.1) +
    ggtitle(paste(lib, "- Hemoglobin Content"))
  
  # Scatter plot: log2 counts vs features colored by hemoglobin
  df_diag <- srat@meta.data
  df_diag$log2_counts <- log2(df_diag$nCount_RNA + 1)
  df_diag$log2_features <- log2(df_diag$nFeature_RNA + 1)
  hemo_plots[[lib]]$scatter <- ggplot(df_diag, aes(x = log2_counts, y = log2_features, color = pct_hemo)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 5), oob = scales::squish) +
    ggtitle(paste(lib, "- Hemoglobin signal"))
  
  # Before/after hemoglobin filtering scatter plots
  p_before <- ggplot(srat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = pct_hemo)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 5), oob = scales::squish) +
    ggtitle(paste(lib, "- Before Hb filtering")) +
    xlab("RNA counts") + ylab("Features")
  
  srat_filtered <- subset(srat, pct_hemo <= thr$max_pct_hemo)
  p_after <- ggplot(srat_filtered@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = pct_hemo)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 5), oob = scales::squish) +
    ggtitle(paste(lib, "- After Hb filtering")) +
    xlab("RNA counts") + ylab("Features")
  
  hemo_plots[[lib]]$before_after <- p_before / p_after
}

# Display hemoglobin plots
grid.arrange(grobs = lapply(hemo_plots, `[[`, "violin"), ncol = 3, nrow = 2,
             top = "Hemoglobin Content Violin Plots")
grid.arrange(grobs = lapply(hemo_plots, `[[`, "scatter"), ncol = 3, nrow = 2,
             top = "Features vs Counts Scatter Plots (Hemoglobin Signal)")

# Display before/after hemoglobin filtering plots
wrap_plots(lapply(hemo_plots[1:3], `[[`, "before_after"), ncol = 3) +
  plot_annotation(title = "Before vs After Hb filtering (Libraries 1–3)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
wrap_plots(lapply(hemo_plots[4:5], `[[`, "before_after"), ncol = 2) +
  plot_annotation(title = "Before vs After Hb filtering (Libraries 4–5)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

# 7. Cell filtering summaries (including Hb) ----
# Summarize cells discarded based on QC thresholds
discard_summary <- list()
for (lib in names(srat_list)) {
  srat <- srat_list[[lib]]
  thr <- CFG$thresholds[[lib]]
  
  # Function to count cells safely (returns 0 if subset is empty)
  safe_count <- function(obj) if (length(colnames(obj)) == 0) 0 else ncol(obj)
  
  discard_summary[[lib]] <- data.frame(
    library = lib,
    total = ncol(srat),
    low_counts = safe_count(subset(srat, nCount_RNA < thr$min_txpts)),
    high_counts = safe_count(subset(srat, nCount_RNA > thr$max_txpts)),
    high_mito = safe_count(subset(srat, percent_mito > thr$max_pctmito)),
    high_hemo = safe_count(subset(srat, pct_hemo > thr$max_pct_hemo)),
    discarded = safe_count(subset(srat, nCount_RNA < thr$min_txpts |
                                    nCount_RNA > thr$max_txpts |
                                    percent_mito > thr$max_pctmito |
                                    pct_hemo > thr$max_pct_hemo))
  )
}
discard_summary_df <- do.call(rbind, discard_summary)
print(discard_summary_df)

# 8. Final filtering ----
for (lib in names(srat_list)) {
  srat <- srat_list[[lib]]
  thr <- CFG$thresholds[[lib]]
  
  # Filter cells based on QC thresholds
  srat <- subset(srat,
                 subset = nCount_RNA >= thr$min_txpts &
                   nCount_RNA <= thr$max_txpts &
                   percent_mito <= thr$max_pctmito &
                   pct_hemo <= thr$max_pct_hemo)
  
  srat_list[[lib]] <- srat
}

# 9. Doublet detection (scDblFinder) ----
dbl_plots <- list()

for (lib in names(srat_list)) {
  srat <- srat_list[[lib]]
  
  # Run scDblFinder to detect doublets
  # Parameters: iter=10 for 10 iterations, includePCs=14 for PCA components, dims=20 for dimensionality
  sce <- scDblFinder(GetAssayData(srat, assay = "RNA", slot = "counts"),
                     iter = 10, includePCs = 14, dims = 20)
  
  # Add doublet scores and classifications to metadata
  srat$scDblFinder.score <- sce$scDblFinder.score
  srat$scDblFinder.class <- sce$scDblFinder.class
  
  # Scatter plot before doublet removal
  p_before <- ggplot(srat@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = scDblFinder.score)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 1), oob = scales::squish) +
    ggtitle(paste(lib, "- Before doublet removal")) +
    xlab("RNA counts") + ylab("Features")
  
  # Filter singlets and create scatter plot
  srat_singlets <- subset(srat, scDblFinder.class == "singlet")
  p_after <- ggplot(srat_singlets@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = scDblFinder.score)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red", limits = c(0, 1), oob = scales::squish) +
    ggtitle(paste(lib, "- After doublet removal")) +
    xlab("RNA counts") + ylab("Features")
  
  # Store plots
  dbl_plots[[lib]] <- p_before / p_after
  
  # Update Seurat object with singlets only
  srat_list[[lib]] <- srat_singlets
}

# Display before/after doublet scatter plots
wrap_plots(dbl_plots[1:3], ncol = 3) +
  plot_annotation(title = "Before vs After Doublet Filtering (Libraries 1–3)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
wrap_plots(dbl_plots[4:5], ncol = 2) +
  plot_annotation(title = "Before vs After Doublet Filtering (Libraries 4–5)",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))

# 10. Save filtered Seurat objects ----
for (lib in names(srat_list)) {
  srat <- srat_list[[lib]]
  saveRDS(srat, file = file.path(CFG$output_dir, paste0(lib, "_1postQC.rds")))
}

# Clean up memory
gc()