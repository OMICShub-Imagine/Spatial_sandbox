library(Seurat)
library(RANN)
library(ClusterR)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(openxlsx)
library(cowplot)

project <- "ProjectName"
sample.names <- list.files("path/to/single_sample_folders")  # OR manually: c("s1", "s2", ...)

# Folder structure
base.path <- "path/to/project_root"
data.dir <- file.path(base.path, project, "results", "merged_samples")
main.result.dir <- file.path(base.path, project, "results", "neighborhood_analysis")
dir.create(main.result.dir, recursive = TRUE, showWarnings = FALSE)

# Marker file
marker.file <- "markerList.txt"
marker.path <- file.path(base.path, project, "data", marker.file)

# Output subfolder
cell.neigh <- 14
cell.window <- 10
result.dir <- file.path(main.result.dir, paste0("CN", cell.neigh, "_window", cell.window))
dir.create(result.dir, showWarnings = FALSE)

# Distance settings
check.dist <- FALSE
pixel.to.micron <- 0.377474
min.dist <- 50

### Load and prepare input data

obj.merged <- readRDS(file.path(data.dir, "merged.rds"))
metadata <- obj.merged[[]]
metadata$EventID <- rownames(metadata)

# Keep relevant metadata
keep.cols <- c("EventID","orig.ident","cell_id","region","tile_num","x","y","x_tile","y_tile","size","merge.annot")
dummy.matrix <- matrix(0, nrow = nrow(metadata), ncol = length(unique(metadata$merge.annot)))
colnames(dummy.matrix) <- unique(metadata$merge.annot)
tmp.CN.input <- cbind(metadata[, keep.cols], dummy.matrix)

# Read marker list
markers <- read.table(marker.path, sep = '\t', header = TRUE)

### Compute neighborhood matrix

CN.input <- do.call(rbind, lapply(sample.names, function(s) {
  sub.tmp <- tmp.CN.input[tmp.CN.input$orig.ident == s, ]
  nn <- nn2(data = sub.tmp[, c("x", "y")], k = cell.window)
  for (j in 1:nrow(nn$nn.idx)) {
    neighbors <- if (check.dist) {
      valid <- which(nn$nn.dists[j, ] * pixel.to.micron <= min.dist)
      nn$nn.idx[j, valid]
    } else {
      nn$nn.idx[j, ]
    }
    freq <- table(sub.tmp[neighbors, "merge.annot"])
    sub.tmp[j, names(freq)] <- freq
  }
  sub.tmp
}))

CN.input.kmean <- CN.input[, unique(metadata$merge.annot)]

### K-means clustering (MiniBatch)

MbatchKm <- MiniBatchKmeans(CN.input.kmean,
                            clusters = cell.neigh,
                            batch_size = 1024,
                            num_init = 3,
                            max_iters = 100,
                            initializer = "kmeans++",
                            early_stop_iter = 10,
                            seed = 1)

### Enrichment heatmap

clus.centroids <- MbatchKm$centroids
tissue.fract <- table(metadata$merge.annot) / nrow(metadata)
scale.centroids <- t(apply(t(clus.centroids) + tissue.fract, 2, function(x) x / sum(x)))
data.to.heatmap <- t(apply(scale.centroids, 1, function(x) x / tissue.fract))
rownames(data.to.heatmap) <- as.character(1:cell.neigh)

# Plot
suffix <- ifelse(check.dist, "_filter", "_NOfilter")
pdf(file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, "_identification", suffix, ".pdf")),
    width = 10, height = 6)
colors <- c(-7, seq(-6, 6, by = 0.01), 7)
my_palette <- c("blue", colorRampPalette(c("blue", "white", "red"))(length(colors) - 3), "red")
pheatmap(log2(data.to.heatmap), scale = "none", color = my_palette, breaks = colors,
         cluster_rows = TRUE, cluster_cols = TRUE, border_color = "grey90",
         main = "Cell type enrichment in CNs")
dev.off()

### Assign CNs to metadata

cell.CN <- predict_MBatchKMeans(CN.input.kmean, MbatchKm$centroids)
obj.merged <- AddMetaData(obj.merged, cell.CN, col.name = paste0("neighborhood", cell.neigh))

### Per-sample CN plots

for (sample in sample.names) {
  df <- obj.merged@meta.data[obj.merged$orig.ident == sample, ]
  p <- ggplot(df, aes(x = x, y = y, color = as.factor(df[[paste0("neighborhood", cell.neigh)]]))) +
    geom_point(size = 0.3) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(sample)
  ggsave(filename = file.path(result.dir, paste0(sample, "_CN_plot", suffix, ".pdf")), plot = p, width = 4, height = 4)
}

### Outs

saveRDS(obj.merged, file = file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, suffix, ".rds")))
out.CN <- data.frame(metadata[, keep.cols], CN.input.kmean, cell.CN)
colnames(out.CN) <- c(keep.cols, unique(metadata$merge.annot), paste0("neighborhood", cell.neigh))
write.csv(out.CN, file = file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, suffix, ".csv")), row.names = FALSE)

wb <- createWorkbook(title = "results")
modifyBaseFont(wb, fontSize = 12)
addWorksheet(wb, "CN assignment")
writeData(wb, "CN assignment", out.CN)
saveWorkbook(wb, file.path(result.dir, paste0("CN_k", cell.neigh, "_w", cell.window, suffix, ".xlsx")), overwrite = TRUE)
