##########################################################################################
# Script to reproduce manually merged cell population labels from Nowicka et al. (2017)
# for 'BCR-XL' data set
# 
# See main script 'prepare_data_BCR_XL_sim_main.R' for more details.
# 
# Lukas Weber, October 2017
##########################################################################################


DIR_BENCHMARK <- "../../../../../benchmark_data/BCR_XL_sim"
DIR_OUTPUT <- file.path(DIR_BENCHMARK, "population_IDs")
DIR_TMP <- file.path(DIR_BENCHMARK, "population_IDs/tmp")

DIR_CURRENT <- getwd()




####################################################################
# Run code from first few sections in Nowicka et al. (2017) workflow
####################################################################

setwd(DIR_TMP)


# -----------
# Data import
# -----------

library(readxl)

## Load metadata
url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
metadata_filename <- "PBMC8_metadata.xlsx"
download.file(file.path(url, metadata_filename), destfile = metadata_filename)

md <- read_excel(metadata_filename)

## Make sure condition variables are factors with the right levels
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
head(data.frame(md))

## Define colors for conditions
color_conditions <- c("#6A3D9A", "#FF7F00")
names(color_conditions) <- levels(md$condition)


## Load .fcs files
fcs_filename <- "PBMC8_fcs_files.zip"
download.file(file.path(url, fcs_filename), destfile = fcs_filename)
unzip(fcs_filename)

library(flowCore)

fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)
fcs_raw


## Load panel file
panel_filename <- "PBMC8_panel.xlsx"
download.file(file.path(url, panel_filename), destfile = panel_filename)
panel <- read_excel(panel_filename)

head(data.frame(panel))

# Replace problematic characters
panel$Antigen <- gsub("-", "_", panel$Antigen)
panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)

# Replace problematic characters
panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)

# Lineage markers
(lineage_markers <- panel$Antigen[panel$Lineage == 1])

# Functional markers
(functional_markers <- panel$Antigen[panel$Functional == 1])

# Spot checks
all(lineage_markers %in% panel_fcs$desc)
all(functional_markers %in% panel_fcs$desc)


# -------------------
# Data transformation
# -------------------

## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor=5){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[, c(lineage_markers, functional_markers)] / cofactor)
  exprs(x) <- expr
  x
})
fcs

## Extract expression
expr <- fsApply(fcs, exprs)
dim(expr)

## Normalization (for visualization only)
library(matrixStats)

rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1


# ----------------
# Diagnostic plots
# ----------------

## Per-sample marker expression distributions

## Generate sample IDs corresponding to each cell in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

library(ggplot2)
library(reshape2)

ggdf <- data.frame(sample_id = sample_ids, expr)
ggdf <- melt(ggdf, id.var = "sample_id", value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = expression, color = condition, group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_conditions)


## Number of cells per sample

cell_table <- table(sample_ids)

ggdf <- data.frame(sample_id = names(cell_table), cell_counts = as.numeric(cell_table))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = condition)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 2.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_conditions, drop = FALSE) +
  scale_x_discrete(drop = FALSE)


## Multi-dimensional scaling (MDS) plot

# Get the median marker expression per sample
library(dplyr)

expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>%
  summarize_each(funs(median))
expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

library(limma)

mds <- plotMDS(expr_median_sample, plot = FALSE)

library(ggrepel)

ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, sample_id = colnames(expr_median_sample))
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  scale_color_manual(values = color_conditions)


## Heatmap of median marker intensities

library(RColorBrewer)
library(pheatmap)

# Column annotation for the heatmap
mm <- match(colnames(expr_median_sample), md$sample_id)
annotation_col <- data.frame(condition = md$condition[mm],
                             row.names = colnames(expr_median_sample))
annotation_colors <- list(condition = color_conditions)

# Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)

pheatmap(expr_median_sample, color = color, display_numbers = TRUE,
         number_color = "black", fontsize_number = 5, annotation_col = annotation_col,
         annotation_colors = annotation_colors, clustering_method = "average")


# --------------------------------------------------
# Marker ranking based on non-redundancy score (NRS)
# --------------------------------------------------

## Define a function that calculates the NRS per sample
NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  score <- rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
  return(score)
}

## Calculate the score
nrs_sample <- fsApply(fcs[, lineage_markers], NRS, use.exprs = TRUE)
rownames(nrs_sample) <- md$sample_id
nrs <- colMeans(nrs_sample, na.rm = TRUE)

## Plot the NRS for ordered markers
lineage_markers_ord <- names(sort(nrs, decreasing = TRUE))
nrs_sample <- data.frame(nrs_sample)
nrs_sample$sample_id <- rownames(nrs_sample)

ggdf <- melt(nrs_sample, id.var = "sample_id",
             value.name = "nrs", variable.name = "antigen")

ggdf$antigen <- factor(ggdf$antigen, levels = lineage_markers_ord)
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$condition[mm]

ggplot(ggdf, aes(x = antigen, y = nrs)) +
  geom_point(aes(color = condition), alpha = 0.9,
             position = position_jitter(width = 0.3, height = 0)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = color_conditions)


# --------------------------------------------------------------------
# Cell population identification with FlowSOM and ConsensusClusterPlus
# --------------------------------------------------------------------

## FlowSOM clustering
library(FlowSOM)

fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)

## Metaclustering into 20 clusters with ConsensusClusterPlus
library(ConsensusClusterPlus)

codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 20

mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png",
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                           distance = "euclidean", seed = 1234)

## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]


## Heatmap

color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

plot_clustering_heatmap_wrapper <- function(expr, expr01,
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  pheatmap(expr_heat, color = color,
           cluster_cols = FALSE, cluster_rows = cluster_rows,
           labels_col = labels_col, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           annotation_legend = annotation_legend)
}

plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord],
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)


## Density plots

plot_clustering_distr_wrapper <- function(expr, cell_clustering){
  
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  cell_clustering <- factor(cell_clustering,
                            levels = levels(cell_clustering)[cluster_rows$order])
  
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = paste0(levels(cell_clustering), " \n(", freq_clust, "%)"))
  
  ggd <- melt(data.frame(cluster = cell_clustering, expr),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  
  ggplot(data = ggd, aes(x = expression, y = ..scaled..)) +
    geom_density(data = transform(ggd, cluster = NULL),
                 color = "darkgrey", fill = "black", adjust = 1, alpha = 0.3) +
    geom_density(color = "blue", fill = "blue", adjust = 1, alpha = 0.3) +
    facet_grid(cluster ~ antigen, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text = element_text(size = 5),
          strip.text = element_text(size = 5))
}

# plot_clustering_distr_wrapper(expr = expr[, lineage_markers_ord],
#                               cell_clustering = cell_clustering1)


## Visual representation with t-SNE

## Find and skip duplicates
dups <- which(!duplicated(expr[, lineage_markers]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)

## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 2000)

## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)

tsne_expr <- expr[tsne_inds, lineage_markers]

## Run t-SNE
library(Rtsne)

set.seed(1234)
tsne_out <- Rtsne(tsne_expr, check_duplicates = FALSE, pca = FALSE)

## Plot t-SNE colored by CD4 expression
dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                 expr[tsne_inds, lineage_markers])

ggplot(dr, aes(x = tSNE1, y = tSNE2, color = CD4)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_gradientn("CD4",
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

dr$sample_id <- sample_ids[tsne_inds]
mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]
dr$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

## Plot t-SNE colored by clusters
ggp <- ggplot(dr, aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggp

## Facet per sample
ggp + facet_wrap(~ sample_id)

## Facet per condition
ggp + facet_wrap(~ condition)


## Meta-clustering

## Get code sizes; sometimes not all the codes have mapped cells so they will have size 0
code_sizes <- table(factor(som$map$mapping[, 1], levels = 1:nrow(codes)))
code_sizes <- as.numeric(code_sizes)

## Run t-SNE on codes
set.seed(1234)
tsne_out <- Rtsne(codes, pca = FALSE)

## Run PCA on codes
pca_out <- prcomp(codes, center = TRUE, scale. = FALSE)

codes_dr <- data.frame(tSNE1 = tsne_out$Y[, 1], tSNE2 = tsne_out$Y[, 2],
                       PCA1 = pca_out$x[, 1], PCA2 = pca_out$x[, 2])
codes_dr$code_clustering1 <- factor(code_clustering1)
codes_dr$size <- code_sizes

## Plot t-SNE on codes
gg_tsne_codes <- ggplot(codes_dr, aes(x = tSNE1, y = tSNE2,
                                      color = code_clustering1, size = size)) +
  geom_point(alpha = 0.9) +
  theme_bw() +
  scale_color_manual(values = color_clusters, guide = FALSE) +
  theme(legend.position = "bottom")

## Plot PCA on codes
gg_pca_codes <- ggplot(codes_dr, aes(x = PCA1, y = PCA2,
                                     color = code_clustering1, size = size)) +
  geom_point(alpha = 0.9) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), nrow = 2)) +
  scale_size(guide = FALSE) +
  theme(legend.position = "top")

library(cowplot)

plot_grid(gg_tsne_codes, gg_pca_codes, ncol = 1, labels = c('A', 'B'))


# ------------------------
# Cluster merging (manual)
# ------------------------

## Download manual merging scheme from Nowicka et al. (2017)
cluster_merging1_filename <- "PBMC8_cluster_merging1.xlsx"
download.file(file.path(url, cluster_merging1_filename),
              destfile = cluster_merging1_filename)
cluster_merging1 <- read_excel(cluster_merging1_filename)
data.frame(cluster_merging1)

## Convert to factor with merged clusters in correct order
levels_merged <- c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells", "CD8 T-cells", 
                   "DC", "NK cells", "monocytes", "surface-")
cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster, 
                                       levels = levels_merged)

## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

mm <- match(code_clustering1, cluster_merging1$original_cluster)
code_clustering1m <- cluster_merging1$new_cluster[mm]


## t-SNE plot: merged clusters

dr$cell_clustering1m <- cell_clustering1m[tsne_inds]

ggplot(dr, aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

## Heatmap: merging scheme
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord],
                                expr01 = expr01[, lineage_markers_ord], cell_clustering = cell_clustering1,
                                color_clusters = color_clusters, cluster_merging = cluster_merging1)

## Heatmap: merged clusters
plot_clustering_heatmap_wrapper(expr = expr[, lineage_markers_ord], 
                                expr01 = expr01[, lineage_markers_ord], 
                                cell_clustering = cell_clustering1m, 
                                color_clusters = color_clusters)


# -------------------------------------------------------------------
# Stop workflow here (for complete workflow, see Nowicka et al. 2017)
# -------------------------------------------------------------------

setwd(DIR_CURRENT)




###################################
# Save population IDs for each cell
###################################

# note: using manually merged population IDs from Nowicka et al. (2017) workflow


# population labels
head(cell_clustering1m)
str(cell_clustering1m)
length(cell_clustering1m)

# corresponding sample IDs
head(sample_ids)
str(sample_ids)
length(sample_ids)

length(cell_clustering1m) == length(sample_ids)

# note order of sample IDs
unique(sample_ids)


# split by sample
df_labels <- data.frame(sample = sample_ids, population = cell_clustering1m)
head(df_labels)

labels <- split(df_labels$population, df_labels$sample)

# note sample IDs are automatically rearranged alphabetically
names(labels)

# convert sample IDs to same format as .fcs filenames
names(labels) <- gsub("BCRXL([0-9]+)", "patient\\1_BCR-XL", names(labels))
names(labels) <- gsub("Ref([0-9]+)", "patient\\1_Reference", names(labels))

names(labels)


# save as .csv files (one file per sample)
for (i in 1:length(labels)) {
  fn <- file.path(DIR_OUTPUT, paste0("population_IDs_", names(labels)[i], ".csv"))
  write.csv(data.frame(population = labels[[i]]), file = fn, row.names = FALSE)
}



###################################
# Save timestamp file for Makefiles
###################################

file_timestamp <- file.path(DIR_OUTPUT, "timestamp.txt")

sink(file_timestamp)
Sys.time()
sink()



