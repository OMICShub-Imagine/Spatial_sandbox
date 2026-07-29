library(Seurat)
library(dplyr)
library(pheatmap)

# ob.list <- list of Seurat objects

#Define morphology feature names
morpho_feats <- c(
  "Nucleus..Area.µm.2", "Nucleus..Length.µm", "Nucleus..Circularity",
  "Nucleus..Solidity", "Nucleus..Max.diameter.µm", "Nucleus..Min.diameter.µm",
  "Nucleus.Cell.area.ratio", "Cell..Area.µm.2", "Cell..Length.µm",
  "Cell..Circularity", "Cell..Solidity", "Cell..Max.diameter.µm",
  "Cell..Min.diameter.µm"
)

#Extract combined dataframe (morphology + markers) from a single sample
extract_combined_df <- function(obj, morpho_feats) {
  # Extract morphology from metadata
  morpho_df <- obj@meta.data %>% select(all_of(morpho_feats))
  
  # Extract normalized marker intensities (cells x markers)
  markers <- GetAssayData(obj, assay = "Spatial", slot = "data") %>% as.matrix()
  
  # Check rownames matching
  if(!all(rownames(morpho_df) == colnames(markers))) {
    stop("Row names of morphology and marker columns don't match!")
  }
  
  # Combine morphology and transposed markers (cells x features)
  combined_df <- cbind(morpho_df, t(markers)) %>% as.data.frame()
  
  return(combined_df)
}

#Correlation matrix between morpho and markers per sample
compute_cor_matrix <- function(sample_df, morpho_feats) {
  morpho_numeric <- sample_df[, morpho_feats, drop=FALSE]
  marker_feats <- setdiff(colnames(sample_df), morpho_feats)
  marker_numeric <- sample_df[, marker_feats, drop=FALSE]
  
  # Convert all to numeric (safe conversion)
  morpho_numeric <- as.data.frame(lapply(morpho_numeric, function(x) as.numeric(as.character(x))))
  marker_numeric <- as.data.frame(lapply(marker_numeric, function(x) as.numeric(as.character(x))))
  
  # Remove any columns all NA
  morpho_numeric <- morpho_numeric[, colSums(is.na(morpho_numeric)) < nrow(morpho_numeric), drop=FALSE]
  marker_numeric <- marker_numeric[, colSums(is.na(marker_numeric)) < nrow(marker_numeric), drop=FALSE]
  
  # Spearman correlation
  cor_mat <- cor(morpho_numeric, marker_numeric, method = "spearman", use = "pairwise.complete.obs")
  
  return(cor_mat)
}

#Extract combined dfs for all samples
samples_list <- lapply(ob.list, extract_combined_df, morpho_feats = morpho_feats)

#Calculate correlation matrices for each sample
cor_matrices <- lapply(samples_list, compute_cor_matrix, morpho_feats = morpho_feats)

#Align correlation matrices by common features (morpho and marker)
common_morpho <- Reduce(intersect, lapply(cor_matrices, rownames))
common_marker <- Reduce(intersect, lapply(cor_matrices, colnames))

cor_matrices_aligned <- lapply(cor_matrices, function(mat) mat[common_morpho, common_marker])

#Average correlation matrices across samples
avg_cor_matrix <- Reduce("+", cor_matrices_aligned) / length(cor_matrices_aligned)

#Heatmap output
basic_redblue <- colorRampPalette(c("blue", "white", "red"))(100)

pheatmap(avg_cor_matrix,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Average Spearman Correlation Across Samples\nMorphology vs Markers",
         fontsize_row = 10,
         fontsize_col = 10,
         color = basic_redblue,
         breaks = seq(-1, 1, length.out = 101))
         
         
