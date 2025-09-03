library(Seurat)
library(ggplot2)
library(tidyverse)
library(paletteer)

set.seed(1000)
CFG <- list()
CFG$work_dir <- "/Users/oscargamarrasiapo/Desktop/Fetahu_2023/analysis"

integrated <- readRDS(file = file.path(CFG$work_dir,
                                       "integrated_cellsidentified.rds"))

# Find markers
de_genes <- FindAllMarkers(integrated, assay = "RNA", only.pos = TRUE,
                           min.pct = 0.1, logfc.threshold = log(1.5))

# 2. Get top 5 per cluster
top5 <- de_genes %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# 3. Dot plot of clusters × genes
DotPlot(integrated, features = unique(top5$gene)) +
  RotatedAxis() +
  labs(y = "Clusters", x = "Genes") + 
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  scale_color_gradient(
    low = "#3A9BFF",
    high = "#FF5E3A"
  )


#To see altogether both ID clust and cell type:
integrated$cluster_id_celltype <- paste0(
  Idents(integrated), " | ", integrated$cluster_celltype
)

DotPlot(integrated, 
        features = unique(top5$gene), 
        group.by = "cluster_id_celltype") +
  RotatedAxis() +
  labs(y = "Cluster | Cell type", x = "Genes") +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  scale_color_gradient(
    low = "#3A9BFF", 
    high = "#FF5E3A"
  )
