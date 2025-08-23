###############################################################
# CODEX Neighborhood Analysis
# USAGE NOTE
#
# This routine is designed for merged Seurat objects generated
# from CODEX datasets
#
# Requirements:
#   • Seurat version ≥ 3.0 (tested with v4 & v5)
#   • Input must be a merged Seurat object containing one row per cell
#   • The following metadata fields are required in obj@meta.data:
#       - EventID (unique row ID, usually rownames)
#       - orig.ident (sample ID)
#       - cell_id (unique cell label per sample)
#       - region, tile_num (image segmentation info)
#       - x, y (numeric spatial coordinates in pixels)
#       - x_tile, y_tile (tile-based coordinates if available)
#       - size (cell size, in pixels or µm²)
#       - merge.annot (cell type annotation, categorical)
#
# Notes:
#   • The analysis is metadata-driven (does not use assays)
#   • Not compatible with Visium Seurat objects
###############################################################

library(Seurat)        
library(RANN)          # Fast nearest-neighbor search
library(ClusterR)      # MiniBatch K-means clustering
library(pheatmap)      
library(RColorBrewer)  
library(scales)        
library(ggplot2)      
library(ggrepel)
library(corrplot)
library(openxlsx)
library(cowplot)

###############################################################
# 1. Setup project directories & parameters
###############################################################

project <- "ProjectName"

sample.names <- list.files("path/to/single_sample_folders")

# Define folder structure
base.path <- "path/to/project_root"
data.dir <- file.path(base.path, project, "results", "merged_samples")
main.result.dir <- file.path(base.path, project, "results", "neighborhood_analysis")
dir.create(main.result.dir, recursive = TRUE, showWarnings = FALSE)

# Marker file with cell type marker definitions
marker.file <- "markerList.txt"
marker.path <- file.path(base.path, project, "data", marker.file)

# Analysis parameters
cell.neigh <- 14    # number of CN clusters (k in k-means)
cell.window <- 10   # number of neighbors considered per cell
result.dir <- file.path(main.result.dir, paste0("CN", cell.neigh, "_window", cell.window))
dir.create(result.dir, showWarnings = FALSE)

# Distance settings
check.dist <- FALSE           # whether to filter neighbors by distance
pixel.to.micron <- 0.377474   # conversion factor (pixel -> micron)
min.dist <- 50                # minimum distance threshold (used if check.dist=TRUE)

###############################################################
# 2. Load merged Seurat object & prepare metadata
###############################################################

# Load merged Seurat object (contains all samples merged together)
obj.merged <- readRDS(file.path(data.dir, "merged.rds"))

# Extract metadata
metadata <- obj.merged[[]]
metadata$EventID <- rownames(metadata)

# Select key metadata columns to keep
keep.cols <- c("EventID","orig.ident","cell_id","region","tile_num",
               "x","y","x_tile","y_tile","size","merge.annot")

# Create a "dummy matrix" for neighborhood frequency calculation
# - One column per unique cell type (merge.annot categories)
# - All initialized to zero
dummy.matrix <- matrix(0, nrow = nrow(metadata), ncol = length(unique(metadata$merge.annot)))
colnames(dummy.matrix) <- unique(metadata$merge.annot)

# Bind metadata with dummy CN frequency matrix
tmp.CN.input <- cbind(metadata[, keep.cols], dummy.matrix)

# Load marker list for annotation (if needed later)
markers <- read.table(marker.path, sep = '\t', header = TRUE)

###############################################################
# 3. Compute neighborhood matrices per sample
###############################################################

CN.input <- do.call(rbind, lapply(sample.names, function(s) {
  # Subset metadata for one sample
  sub.tmp <- tmp.CN.input[tmp.CN.input$orig.ident == s, ]
  
  # Find nearest neighbors (using RANN::nn2)
  nn <- nn2(data = sub.tmp[, c("x", "y")], k = cell.window)
  
  # For each cell, count the frequency of cell types in its neighborhood
  for (j in 1:nrow(nn$nn.idx)) {
    neighbors <- if (check.dist) {
      # Optionally filter by physical distance (microns)
      valid <- which(nn$nn.dists[j, ] * pixel.to.micron <= min.dist)
      nn$nn.idx[j, valid]
    } else {
      nn$nn.idx[j, ]
    }
    # Count frequencies of neighboring cell annotations
    freq <- table(sub.tmp[neighbors, "merge.annot"])
    sub.tmp[j, names(freq)] <- freq
  }
  sub.tmp
}))

# Extract only the CN frequency matrix (drop metadata columns)
CN.input.kmean <- CN.input[, unique(metadata$merge.annot)]

###############################################################
# 4. Neighborhood clustering using MiniBatch K-means
###############################################################

MbatchKm <- MiniBatchKmeans(CN.input.kmean,
                            clusters = cell.neigh,   # number of CNs
                            batch_size = 1024,
                            num_init = 3,            # how many random initializations
                            max_iters = 100,
                            initializer = "kmeans++",
                            early_stop_iter = 10,
                            seed = 1)

###############################################################
# 5. Enrichment heatmap of cell types per CN
###############################################################

# Get cluster centroids (CN profiles)
clus.centroids <- MbatchKm$centroids

# Background tissue frequencies of each cell type
tissue.fract <- table(metadata$merge.annot) / nrow(metadata)

# Normalize CN profiles by tissue abundance
scale.centroids <- t(apply(t(clus.centroids) + tissue.fract, 2, function(x) x / sum(x)))

# Compute enrichment score: relative abundance vs tissue background
data.to.heatmap <- t(apply(scale.centroids, 1, function(x) x / tissue.fract))
rownames(data.to.heatmap) <- as.character(1:cell.neigh)

# Save enrichment heatmap as PDF
suffix <- ifelse(check.dist, "_filter", "_NOfilter")
pdf(file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, "_identification", suffix, ".pdf")),
    width = 10, height = 6)

# Define color scale (log2 fold enrichment)
colors <- c(-7, seq(-6, 6, by = 0.01), 7)
my_palette <- c("blue", colorRampPalette(c("blue", "white", "red"))(length(colors) - 3), "red")

# Plot heatmap
pheatmap(log2(data.to.heatmap), scale = "none", color = my_palette, breaks = colors,
         cluster_rows = TRUE, cluster_cols = TRUE, border_color = "grey90",
         main = "Cell type enrichment in CNs")
dev.off()

###############################################################
# 6. Assign CN cluster labels to individual cells
###############################################################

# Predict cluster assignment for each cell
cell.CN <- predict_MBatchKMeans(CN.input.kmean, MbatchKm$centroids)

# Add CN assignment to Seurat object metadata
obj.merged <- AddMetaData(obj.merged, cell.CN, col.name = paste0("neighborhood", cell.neigh))

###############################################################
# 7. Per-sample CN spatial visualization
###############################################################

for (sample in sample.names) {
  df <- obj.merged@meta.data[obj.merged$orig.ident == sample, ]
  
  # Plot cell positions colored by CN cluster
  p <- ggplot(df, aes(x = x, y = y, color = as.factor(df[[paste0("neighborhood", cell.neigh)]]))) +
    geom_point(size = 0.3) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(sample)
  
  # Save per-sample CN plot
  ggsave(filename = file.path(result.dir, paste0(sample, "_CN_plot", suffix, ".pdf")), 
         plot = p, width = 4, height = 4)
}

###############################################################
# 8. Save results (RDS, CSV, Excel)
###############################################################

# Save Seurat object with CN assignments
saveRDS(obj.merged, file = file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, suffix, ".rds")))

# Export full CN assignment table
out.CN <- data.frame(metadata[, keep.cols], CN.input.kmean, cell.CN)
colnames(out.CN) <- c(keep.cols, unique(metadata$merge.annot), paste0("neighborhood", cell.neigh))
write.csv(out.CN, file = file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, suffix, ".csv")), row.names = FALSE)

# Export results to Excel
wb <- createWorkbook(title = "results")
modifyBaseFont(wb, fontSize = 12)
addWorksheet(wb, "CN assignment")
writeData(wb, "CN assignment", out.CN)
saveWorkbook(wb, file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, suffix, ".xlsx")), overwrite = TRUE)
