library(Seurat)
library(sva)        # ComBat-seq
library(harmony)
library(mclust)     # ARI
library(cluster)    # silhouette
library(pheatmap)
library(ggplot2)
library(patchwork)
library(dplyr)

# ---------------------------
# Helper function silhouette on subsample

silhouette_subsample <- function(embedding_matrix, labels, subsample_size = 1000, seed = 42) {
  set.seed(seed)
  n_cells <- nrow(embedding_matrix)
  if (n_cells <= subsample_size) {
    warning("Subsample size >= number of cells; using full data")
    subsample_idx <- seq_len(n_cells)
  } else {
    subsample_idx <- sample(n_cells, subsample_size)
  }
  subsample_embeddings <- embedding_matrix[subsample_idx, , drop = FALSE]
  subsample_labels <- labels[subsample_idx]
  dist_mat <- dist(subsample_embeddings)
  sil <- silhouette(as.numeric(as.factor(subsample_labels)), dist_mat)
  avg_sil <- summary(sil)$avg.width
  list(silhouette = sil, average = avg_sil, indices = subsample_idx)
}

# ---------------------------
# Load

files <- list.files(pattern = "M[12]_L[1-4]_B[1-6]\\.rds")
ob.list <- lapply(files, function(file) {
  obj <- readRDS(file)
  spatial_assay <- GetAssay(obj, "Spatial")
  
  # Capitalize marker names
  features_cap <- toupper(rownames(spatial_assay))
  rownames(spatial_assay) <- features_cap
  obj[["Spatial"]] <- spatial_assay
  DefaultAssay(obj) <- "Spatial"
  
  # Add metadata
  obj$Batch <- ifelse(grepl("^M1", obj$orig.ident), "M1",
                      ifelse(grepl("^M2", obj$orig.ident), "M2", NA))
  if (is.null(obj$Biopsy)) {
    # Optional: define Biopsy from filename or ident
    obj$Biopsy <- gsub("\\.rds", "", file)
  }
  
  obj
})

# ---------------------------
# Prepare raw counts and metadata

DefaultAssay(ob.list[[1]]) <- "Spatial"
raw_counts_list <- lapply(ob.list, function(obj) {
  DefaultAssay(obj) <- "Spatial"
  GetAssayData(obj, assay = "Spatial", layer = "counts")
})

# Get common features
common_features <- Reduce(intersect, lapply(raw_counts_list, rownames))
raw_counts_list <- lapply(raw_counts_list, function(m) m[common_features, , drop = FALSE])
combat_seq_counts <- do.call(cbind, raw_counts_list)

# Merge metadata
meta_list <- lapply(ob.list, function(obj) obj@meta.data)
combat_seq_meta <- do.call(rbind, meta_list)
combat_seq_meta$Batch <- as.factor(combat_seq_meta$Batch)
combat_seq_meta$Biopsy <- as.factor(combat_seq_meta$Biopsy)

# ---------------------------
# ComBat-seq

combat_seq_corrected <- ComBat_seq(
  counts = as.matrix(combat_seq_counts),
  batch = combat_seq_meta$Batch,
  group = NULL
)

# Build corrected Seurat object
seurat_combat_seq <- CreateSeuratObject(counts = combat_seq_corrected, assay = "Spatial", meta.data = combat_seq_meta)
LayerData(seurat_combat_seq[["Spatial"]], layer = "data") <- LayerData(seurat_combat_seq[["Spatial"]], layer = "counts")

# Dimensionality Reduction
seurat_combat_seq <- FindVariableFeatures(seurat_combat_seq, assay = "Spatial", selection.method = "vst")
#seurat_combat_seq <- ScaleData(seurat_combat_seq, assay = "Spatial", features = VariableFeatures(seurat_combat_seq))
seurat_combat_seq <- ScaleData(seurat_combat_seq, assay = "Spatial", features = VariableFeatures(seurat_combat_seq),
							   vars.to.regress = c("Nucleus.Cell.area.ratio","Cell..Solidity"))
seurat_combat_seq <- RunPCA(seurat_combat_seq, assay = "Spatial", npcs = 10)
seurat_combat_seq <- RunUMAP(seurat_combat_seq, dims = 1:10, reduction.name = "umap_after_combatseq")

seurat_combat_seq <- FindNeighbors(seurat_combat_seq, dims = 1:10)
seurat_combat_seq <- FindClusters(seurat_combat_seq, resolution = 0.1)

saveRDS(seurat_combat_seq, file = "CODEX_combat_seq.rds")

# UMAP plot by Batch and Biopsy
p_batch <- DimPlot(seurat_combat_seq, reduction = "umap_after_combatseq", group.by = "Batch") +
  ggtitle("ComBat-seq: UMAP by Batch")

p_biopsy <- DimPlot(seurat_combat_seq, reduction = "umap_after_combatseq", group.by = "Biopsy") +
  ggtitle("ComBat-seq: UMAP by Biopsy")

print(p_batch + p_biopsy + plot_layout(ncol = 2))

# Extract UMAP embeddings
embedding_combat_seq <- Embeddings(seurat_combat_seq, "umap_after_combatseq")

# Compute silhouette on UMAP using Batch
sil_combat_seq <- silhouette_subsample(
  embedding_matrix = embedding_combat_seq,
  labels = seurat_combat_seq$Batch,
  subsample_size = 10000
)

# Adjusted Rand Index (ARI) clustering vs. known cell types (Classification)
ari_combat_seq <- adjustedRandIndex(seurat_combat_seq$seurat_clusters, seurat_combat_seq$Classification)

# Contingency table and purity
contingency_combat_seq <- table(seurat_combat_seq$seurat_clusters, seurat_combat_seq$Classification)
purity_combat_seq <- sum(apply(contingency_combat_seq, 1, max)) / sum(contingency_combat_seq)

# Plot cluster vs cell type heatmap
pheatmap(
  contingency_combat_seq,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "ComBat-seq: Clusters vs Cell Types"
)

# ---------------------------
# asinh + Harmony 

ob.list_asinh <- lapply(ob.list, function(obj) {
  DefaultAssay(obj) <- "Spatial"
  counts <- GetAssayData(obj, assay = "Spatial", slot = "counts")
  transformed <- asinh(counts / 5)
  obj[["Spatial"]] <- CreateAssayObject(counts = transformed)
  obj
})

seurat_harmony <- merge(x = ob.list_asinh[[1]], y = ob.list_asinh[-1])
seurat_harmony <- FindVariableFeatures(seurat_harmony, assay = "Spatial", selection.method = "vst")
seurat_harmony <- ScaleData(seurat_harmony, assay = "Spatial", features = VariableFeatures(seurat_harmony),
							vars.to.regress = c("Nucleus.Cell.area.ratio","Cell..Solidity"))
seurat_harmony <- RunPCA(seurat_harmony, assay = "Spatial", features = VariableFeatures(seurat_harmony), npcs = 10)
seurat_harmony <- RunUMAP(seurat_harmony, dims = 1:10, reduction.name = "umap_before_harmony")

p_before_harmony <- DimPlot(seurat_harmony, reduction = "umap_before_harmony", group.by = "Batch") +
  ggtitle("Harmony: UMAP Before Integration")
print(p_before_harmony)

seurat_harmony <- RunHarmony(seurat_harmony, group.by.vars = "Batch", reduction.use = "pca", assay.use = "Spatial")
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:10, reduction.name = "umap_after_harmony")
seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:10)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.2)

saveRDS(seurat_harmony, file = "CODEX_harmony.rds")

p_after_harmony <- DimPlot(seurat_harmony, reduction = "umap_after_harmony", group.by = "Batch") +
  ggtitle("Harmony: UMAP After Integration")
print(p_before_harmony + p_after_harmony + plot_layout(ncol = 2))

embedding_harmony <- Embeddings(seurat_harmony, "umap_after_harmony")
sil_harmony <- silhouette_subsample(embedding_harmony, seurat_harmony$Batch, subsample_size = 10000)
ari_harmony <- adjustedRandIndex(seurat_harmony$seurat_clusters, seurat_harmony$Classification)
contingency_harmony <- table(seurat_harmony$seurat_clusters, seurat_harmony$Classification)
purity_harmony <- sum(apply(contingency_harmony, 1, max)) / sum(contingency_harmony)

pheatmap(contingency_harmony, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Harmony: Clusters vs Cell Types")

# ---------------------------
# CLR normalization + RPCA Integration 

all_features <- Reduce(intersect, lapply(ob.list, function(obj) rownames(GetAssay(obj, "Spatial"))))

ob.list_clr <- lapply(ob.list, function(obj) {
  DefaultAssay(obj) <- "Spatial"
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 1)
  obj <- ScaleData(obj, features = all_features,
                   vars.to.regress = c("Nucleus.Cell.area.ratio", "Cell..Solidity"))
  obj <- RunPCA(obj, features = all_features, npcs = 20)
  obj
})

anchors <- FindIntegrationAnchors(object.list = ob.list_clr, anchor.features = all_features,
                                  reduction = "rpca", dims = 1:9)
integrated <- IntegrateData(anchorset = anchors, dims = 1:9)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 30)
integrated <- RunUMAP(integrated, dims = 1:9)
integrated <- FindNeighbors(integrated, dims = 1:9)
integrated <- FindClusters(integrated, resolution = 0.2)

saveRDS(integrated, file = "CODEX_rpca_integrated.rds")

clr_merged <- merge(x = ob.list_clr[[1]], y = ob.list_clr[-1])
clr_merged <- FindVariableFeatures(clr_merged, assay = "Spatial", selection.method = "vst")
clr_merged <- ScaleData(clr_merged, assay = "Spatial", features = VariableFeatures(clr_merged))
clr_merged <- RunPCA(clr_merged, assay = "Spatial", features = VariableFeatures(clr_merged), npcs = 10)
clr_merged <- RunUMAP(clr_merged, dims = 1:10, reduction.name = "umap_before_rpca")

p_before_rpca <- DimPlot(clr_merged, reduction = "umap_before_rpca", group.by = "Batch") +
  ggtitle("RPCA: UMAP Before Integration")

p_after_rpca <- DimPlot(integrated, reduction = "umap", group.by = "Batch") +
  ggtitle("RPCA: UMAP After Integration")
print(p_before_rpca + p_after_rpca + plot_layout(ncol = 2))

embedding_rpca <- Embeddings(integrated, "umap")
sil_rpca <- silhouette_subsample(embedding_rpca, integrated$Batch, subsample_size = 10000)
ari_rpca <- adjustedRandIndex(integrated$seurat_clusters, integrated$Classification)
contingency_rpca <- table(integrated$seurat_clusters, integrated$Classification)
purity_rpca <- sum(apply(contingency_rpca, 1, max)) / sum(contingency_rpca)

pheatmap(contingency_rpca, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "RPCA: Clusters vs Cell Types")

# ---------------------------
# Summary metrics

summary_df <- data.frame(
  Method = c("ComBat", "Harmony", "RPCA"),
  ARI = c(ari_combat_seq, ari_harmony, ari_rpca),
  Purity = c(purity_combat_seq, purity_harmony, purity_rpca),
  Silhouette_Batch = c(sil_combat_seq$average, sil_harmony$average, sil_rpca$average)
)

print(summary_df)
write.csv(summary_df, "integration_summary_metrics.csv", row.names = FALSE)


# ---------------------------
# Add kBET per method

compute_kBET <- function(embedding_matrix, batch_labels, k0 = NULL, subsample = 5000, seed = 42) {
  set.seed(seed)
  if (nrow(embedding_matrix) > subsample) {
    idx <- sample(seq_len(nrow(embedding_matrix)), subsample)
    embedding_matrix <- embedding_matrix[idx, , drop = FALSE]
    batch_labels <- batch_labels[idx]
  }
  embedding_matrix <- as.matrix(embedding_matrix)
  batch_labels <- as.factor(batch_labels)
  if (is.null(k0)) {
    k0 <- floor(sqrt(nrow(embedding_matrix)))
  }
  kbet_result <- kBET(embedding_matrix, batch_labels, k0 = k0)
  if (!is.null(kbet_result$stats$kBET.signif)) {
    rejection_rate <- mean(kbet_result$stats$kBET.signif, na.rm = TRUE)
  } else {
    rejection_rate <- NA
    warning("kBET result did not contain 'kBET.signif'.")
  }
  return(rejection_rate)
}

# Compute kBET for all methods
batch_labels_combat <- seurat_combat_seq$Batch
names(batch_labels_combat) <- colnames(seurat_combat_seq)
batch_labels_combat <- batch_labels_combat[rownames(embedding_combat_seq)]

kbet_combat_seq <- compute_kBET(embedding_combat_seq, batch_labels_combat)
kbet_harmony <- compute_kBET(embedding_harmony, seurat_harmony$Batch)
kbet_rpca <- compute_kBET(embedding_rpca, integrated$Batch)

# Add to summary
summary_df$kBET <- c(kbet_combat_seq, kbet_harmony, kbet_rpca)

# Save summary
print(summary_df)
write.csv(summary_df, "integration_summary_metrics.csv", row.names = FALSE)

# ---------------------------
# RidgePlots Before/After Integration

# ComBat-seq RidgePlots
RidgePlot(seurat_combat_seq, features = VariableFeatures(seurat_combat_seq)[1:6], group.by = "Batch") +
  ggtitle("ComBat-seq: After Correction")
RidgePlot(ob.list[[1]], features = VariableFeatures(ob.list[[1]])[1:6], group.by = "Batch") +
  ggtitle("ComBat-seq: Before Correction")

# Harmony RidgePlots
RidgePlot(seurat_harmony, features = VariableFeatures(seurat_harmony)[1:6], group.by = "Batch") +
  ggtitle("Harmony: After Correction")
RidgePlot(seurat_harmony, reduction = "umap_before_harmony", features = VariableFeatures(seurat_harmony)[1:6], group.by = "Batch") +
  ggtitle("Harmony: Before Correction")

# RPCA RidgePlots
RidgePlot(integrated, features = VariableFeatures(integrated)[1:6], group.by = "Batch") +
  ggtitle("RPCA: After Correction")
RidgePlot(clr_merged, features = VariableFeatures(clr_merged)[1:6], group.by = "Batch") +
  ggtitle("RPCA: Before Correction")

# ---------------------------
# Top markers

plot_top_markers_density <- function(seurat_obj, reduction_name, method_name, assay = "Spatial") {
  library(ComplexHeatmap)
  library(circlize)
  library(Seurat)
  library(ggplot2)
  library(glue)
  library(patchwork)
  library(dplyr)
  library(grid)

  DefaultAssay(seurat_obj) <- assay
  Idents(seurat_obj) <- "seurat_clusters"

  message(glue("Finding markers for {method_name}..."))
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, logfc.threshold = 0)

  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 5) %>%
    ungroup()

  top_genes <- unique(top_markers$gene)

  # Get average expression (genes x clusters)
  avg_expr_list <- AverageExpression(
    seurat_obj,
    features = top_genes,
    assays = assay,
    return.seurat = FALSE
  )

  if (!(assay %in% names(avg_expr_list))) {
    stop(glue("Assay '{assay}' not found in AverageExpression output. Available: {paste(names(avg_expr_list), collapse = ', ')}"))
  }

  expr_mat <- as.matrix(avg_expr_list[[assay]])
  expr_mat <- apply(expr_mat, 2, as.numeric)
  rownames(expr_mat) <- rownames(avg_expr_list[[assay]])

  # Rename cluster columns from g0,g1,... to actual cluster IDs if needed
  cluster_ids <- levels(Idents(seurat_obj))
  if (all(grepl("^g[0-9]+$", colnames(expr_mat))) && length(cluster_ids) == ncol(expr_mat)) {
    colnames(expr_mat) <- cluster_ids
  }

  # Scale expression per gene across clusters (z-score)
  expr_scaled <- t(scale(t(expr_mat)))
  expr_scaled <- expr_scaled[complete.cases(expr_scaled), , drop = FALSE]

  # Cluster colors
  library(scales)
  cluster_colors <- setNames(hue_pal()(length(cluster_ids)), cluster_ids)

  # Column annotation for clusters
  col_anno <- HeatmapAnnotation(
    Cluster = colnames(expr_scaled),
    col = list(Cluster = cluster_colors),
    annotation_name_side = "left"
  )

  # Color function for heatmap
  col_fun <- colorRamp2(
    c(min(expr_scaled, na.rm = TRUE), 0, max(expr_scaled, na.rm = TRUE)),
    c("navy", "white", "firebrick3")
  )

  heatmap_plot <- Heatmap(
    expr_scaled,
    name = "Scaled Expr",
    top_annotation = col_anno,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Z-score Expression"),
    col = col_fun,
    column_title = paste0(method_name, " - Top Marker Scaled Expression per Cluster"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )

  draw(heatmap_plot, heatmap_legend_side = "right")

  # Plot density for top 6 markers
  top_features <- top_genes[1:min(6, length(top_genes))]

  p_list <- lapply(top_features, function(gene) {
    Plot_Density_Custom(
      seurat_object = seurat_obj,
      features = gene,
      reduction = reduction_name,
      viridis_palette = "magma",
      pt.size = 0.5,
      combine = TRUE
    ) +
      ggtitle(glue("{method_name}: {gene}")) +
      theme(plot.title = element_text(size = 10))
  })

  wrap_plots(p_list, ncol = 3) &
    plot_annotation(title = glue("{method_name} - Top Marker Density Expression"))
}


# ComBat
plot_top_markers_density(seurat_combat_seq, reduction_name = "umap_after_combat", method_name = "ComBat", assay = "Spatial")

# Harmony
plot_top_markers_density(seurat_harmony, 
                         reduction_name = "umap_after_harmony", 
                         method_name = "Harmony", 
                         assay = "Spatial")
                         
# RPCA
plot_top_markers_density(integrated, reduction_name = "umap", method_name = "RPCA", assay = "integrated")
