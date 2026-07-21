#!/usr/bin/env Rscript
#
# Usage:
#   Rscript Xenium_QC.R --xenium-dir /path/to/output-XETG... --sample P1_L
#   Rscript Xenium_QC.R --xenium-dir /path/to/output --sample ROI_B --mode fixed
#   Rscript Xenium_QC.R --xenium-dir /path/to/output --sample P1_NL --mode auto --nmads 3
#
# Modes:
#   auto  - data-driven thresholds (MAD / percentile), clipped to safe bounds (default)
#   fixed - legacy fixed cutoffs from skin 5K samples

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(viridis)
  library(tibble)
  library(readr)
  library(patchwork)
})

# --- CLI ---
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i >= length(args)) return(default)
  args[i + 1]
}

xenium_dir <- get_arg("--xenium-dir")
sample_name <- get_arg("--sample", "Xenium_sample")
out_dir <- get_arg("--out-dir", ".")
mode <- get_arg("--mode", "auto")
nmads <- as.numeric(get_arg("--nmads", "3"))
cluster_res <- as.numeric(get_arg("--cluster-res", "0.2"))
skip_clustering <- "--skip-clustering" %in% args

if (is.null(xenium_dir) || !dir.exists(xenium_dir)) {
  stop(
    "Provide a valid Xenium output directory:\n",
    "  Rscript Xenium_QC.R --xenium-dir /path/to/output-XETG... --sample SAMPLE_ID\n",
    call. = FALSE
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sample_name)

# --- helpers ---
robust_lower <- function(x, nmads = 3, floor_val = 0, ceiling_val = Inf) {
  x <- x[is.finite(x)]
  if (!length(x)) return(floor_val)
  lx <- log10(x + 1)
  bound <- median(lx) - nmads * mad(lx)
  val <- 10^bound - 1
  val <- max(floor_val, val)
  min(val, ceiling_val)
}

robust_upper <- function(x, nmads = 3, floor_val = 0, ceiling_val = Inf) {
  x <- x[is.finite(x)]
  if (!length(x)) return(ceiling_val)
  bound <- median(x) + nmads * mad(x)
  val <- min(quantile(x, 0.975, na.rm = TRUE), bound)
  val <- max(floor_val, val)
  min(val, ceiling_val)
}

compute_qc_metrics <- function(obj) {
  obj$frac_blank <- obj$nCount_BlankCodeword / pmax(obj$nCount_Xenium, 1)
  obj$frac_control <- obj$nCount_ControlProbe / pmax(obj$nCount_Xenium, 1)
  obj$frac_genomic <- obj$nCount_GenomicControl / pmax(obj$nCount_Xenium, 1)
  obj$frac_deprecated <- if ("nCount_DeprecatedCodeword" %in% colnames(obj@meta.data)) {
    obj$nCount_DeprecatedCodeword / pmax(obj$nCount_Xenium, 1)
  } else {
    NA_real_
  }

  mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
  if (length(mt_genes) > 0) {
    obj$percent.mt <- 100 * Matrix::colSums(obj[mt_genes, ]) / pmax(obj$nCount_Xenium, 1)
  }
  obj
}

compute_thresholds <- function(meta, mode = c("auto", "fixed"), nmads = 3) {
  mode <- match.arg(mode)
  fixed <- list(
    min_counts = 100L,
    min_features = 100L,
    max_frac_blank = 0.2,
    max_frac_control = 0.2,
    max_frac_genomic = 0.2,
    max_percent_mt = 25
  )
  if (mode == "fixed") {
    return(c(fixed, list(mode = "fixed", nmads = nmads)))
  }

  thr <- list(
    min_counts = round(robust_lower(meta$nCount_Xenium, nmads = nmads, floor_val = 50, ceiling_val = 200)),
    min_features = round(robust_lower(meta$nFeature_Xenium, nmads = nmads, floor_val = 50, ceiling_val = 200)),
    max_frac_blank = robust_upper(meta$frac_blank, nmads = nmads, floor_val = 0.05, ceiling_val = 0.25),
    max_frac_control = robust_upper(meta$frac_control, nmads = nmads, floor_val = 0.05, ceiling_val = 0.25),
    max_frac_genomic = robust_upper(meta$frac_genomic, nmads = nmads, floor_val = 0.05, ceiling_val = 0.25),
    max_percent_mt = if ("percent.mt" %in% colnames(meta)) {
      robust_upper(meta$percent.mt, nmads = nmads, floor_val = 10, ceiling_val = 25)
    } else {
      NA_real_
    },
    mode = "auto",
    nmads = nmads
  )
  thr
}

apply_filters <- function(meta, thr) {
  keep <- meta$nCount_Xenium > thr$min_counts &
    meta$nFeature_Xenium > thr$min_features &
    meta$frac_blank < thr$max_frac_blank &
    meta$frac_control < thr$max_frac_control &
    meta$frac_genomic < thr$max_frac_genomic
  if ("percent.mt" %in% colnames(meta) && is.finite(thr$max_percent_mt)) {
    keep <- keep & meta$percent.mt < thr$max_percent_mt
  }
  keep
}

save_plot <- function(p, path, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
}

thresholds_to_df <- function(thr) {
  tibble(
    parameter = c(
      "mode", "nmads", "min_counts", "min_features",
      "max_frac_blank", "max_frac_control", "max_frac_genomic", "max_percent_mt"
    ),
    value = c(
      thr$mode, thr$nmads, thr$min_counts, thr$min_features,
      thr$max_frac_blank, thr$max_frac_control, thr$max_frac_genomic,
      ifelse(is.null(thr$max_percent_mt) || !is.finite(thr$max_percent_mt), NA, thr$max_percent_mt)
    )
  )
}

cat("Loading Xenium:", xenium_dir, "\n")
xenium.obj <- LoadXenium(xenium_dir)
DefaultAssay(xenium.obj) <- "Xenium"
pre_ncells <- ncol(xenium.obj)

xenium.obj <- compute_qc_metrics(xenium.obj)
meta_pre <- xenium.obj@meta.data
thr <- compute_thresholds(meta_pre, mode = mode, nmads = nmads)
keep <- apply_filters(meta_pre, thr)

cat("Threshold mode:", thr$mode, "\n")
print(thresholds_to_df(thr))

# --- QC plots (pre-filter) ---
p_count <- ImageFeaturePlot(xenium.obj, features = "nCount_Xenium", size = 0.5) +
  scale_fill_viridis_c(option = "viridis") +
  ggtitle(paste0(sample_name, " — transcripts per cell (pre-QC)")) +
  theme(legend.position = "right")
save_plot(p_count, paste0(prefix, "_QC_spatial_nCount.png"), width = 8, height = 7)

p_feat <- ImageFeaturePlot(xenium.obj, features = "nFeature_Xenium", size = 0.5) +
  scale_fill_viridis_c(option = "plasma") +
  ggtitle(paste0(sample_name, " — features per cell (pre-QC)")) +
  theme(legend.position = "right")
save_plot(p_feat, paste0(prefix, "_QC_spatial_nFeature.png"), width = 8, height = 7)

qc_features <- c("nFeature_Xenium", "nCount_Xenium", "frac_blank", "frac_control", "frac_genomic")
if ("percent.mt" %in% colnames(meta_pre)) qc_features <- c(qc_features, "percent.mt")

p_vln <- VlnPlot(xenium.obj, features = qc_features, ncol = 3, pt.size = 0, combine = FALSE)
for (i in seq_along(p_vln)) {
  p_vln[[i]] <- p_vln[[i]] + ggtitle(qc_features[i])
}
save_plot(
  patchwork::wrap_plots(p_vln, ncol = 3),
  paste0(prefix, "_QC_violin_metrics.png"),
  width = 12, height = 3.5 * ceiling(length(qc_features) / 3)
)

# --- filter ---
xenium.obj$qc_pass <- keep
xenium.filt <- subset(xenium.obj, subset = qc_pass)
post_ncells <- ncol(xenium.filt)
cat(sprintf("Cells kept: %d / %d (%.1f%%)\n", post_ncells, pre_ncells, 100 * post_ncells / pre_ncells))

if (post_ncells < 100) {
  warning("Very few cells after QC (n < 100). Review thresholds and plots.", call. = FALSE)
}

# --- optional clustering / markers ---
all.markers <- tibble()
if (!skip_clustering && post_ncells >= 50) {
  xenium.filt <- NormalizeData(xenium.filt, verbose = FALSE)
  xenium.filt <- FindVariableFeatures(xenium.filt, nfeatures = 2000, verbose = FALSE)
  xenium.filt <- ScaleData(xenium.filt, verbose = FALSE)
  xenium.filt <- RunPCA(xenium.filt, verbose = FALSE)
  xenium.filt <- FindNeighbors(xenium.filt, dims = 1:20, verbose = FALSE)
  xenium.filt <- FindClusters(xenium.filt, resolution = cluster_res, verbose = FALSE)
  xenium.filt <- RunUMAP(xenium.filt, dims = 1:20, verbose = FALSE)

  p_umap <- DimPlot(xenium.filt, label = TRUE) + ggtitle(paste0(sample_name, " — UMAP"))
  save_plot(p_umap, paste0(prefix, "_QC_umap.png"), width = 7, height = 6)

  p_img <- ImageDimPlot(xenium.filt, cols = "polychrome", size = 0.75) +
    ggtitle(paste0(sample_name, " — clusters on tissue"))
  save_plot(p_img, paste0(prefix, "_QC_image_clusters.png"), width = 8, height = 7)

  all.markers <- FindAllMarkers(
    xenium.filt,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  ) %>% as_tibble()
}

# --- summaries ---
qc_summary <- tibble(
  sample = sample_name,
  xenium_dir = normalizePath(xenium_dir),
  pre_qc_cells = pre_ncells,
  post_qc_cells = post_ncells,
  percent_kept = round(100 * post_ncells / pre_ncells, 2),
  threshold_mode = thr$mode,
  nmads = thr$nmads
)

cell_qc <- meta_pre %>%
  rownames_to_column("cell_id") %>%
  mutate(qc_pass = keep) %>%
  dplyr::select(cell_id, nCount_Xenium, nFeature_Xenium, frac_blank, frac_control, frac_genomic,
                any_of("percent.mt"), qc_pass)

thr_df <- thresholds_to_df(thr)

wb <- createWorkbook()
addWorksheet(wb, "QC_summary")
writeData(wb, "QC_summary", qc_summary)
addWorksheet(wb, "Thresholds")
writeData(wb, "Thresholds", thr_df)
addWorksheet(wb, "Cell_QC")
writeData(wb, "Cell_QC", cell_qc)
if (nrow(all.markers) > 0) {
  addWorksheet(wb, "Cluster_Markers")
  writeData(wb, "Cluster_Markers", all.markers)
}
saveWorkbook(wb, paste0(prefix, "_QC_report.xlsx"), overwrite = TRUE)

write_csv(qc_summary, paste0(prefix, "_QC_summary.csv"))
write_csv(thr_df, paste0(prefix, "_QC_thresholds.csv"))
write_csv(cell_qc, paste0(prefix, "_QC_per_cell.csv"))
saveRDS(xenium.filt, paste0(prefix, "_QC_filtered.rds"))

cat("Outputs written with prefix:", prefix, "\n")
