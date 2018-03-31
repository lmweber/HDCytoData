##########################################################################################
# Script to load CyTOF example data set: BCR-XL and reproduce manually merged
# cell population labels from Nowicka et al. (2017)
# 
# The original 'BCR-XL' data set is sourced from Bodenmiller et al. (2012), and has
# previously been used for benchmark evaluations by Bruggner et al. (2014) and Nowicka et
# al. (2017).
#
# Raw data downloaded from Cytobank (experiment 15713)
# - see Citrus wiki (section 'PBMC Example 1'):
# https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - direct link to Cytobank repository:
# https://community.cytobank.org/cytobank/experiments/15713/download_files
#
# Code in this script is adapted from the workflow paper by Nowicka et al. (2017).
# Cell population labels are reproduced from Nowicka et al. (2017), where they were
# generated using a strategy of expert-guided manual merging of automatically generated
# clusters from the FlowSOM algorithm. 
# 
# Lukas Weber, January 2018
##########################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(flowCore)
  library(FlowSOM)
  library(ConsensusClusterPlus)
  library(SummarizedExperiment)
})

DIR_RAW_DATA <- "raw_data"
if (!dir.exists(DIR_RAW_DATA)) dir.create(DIR_RAW_DATA)

####################################################################
# Run code from first few sections in Nowicka et al. (2017) workflow
####################################################################

# -----------
# Data import
# -----------
## Load metadata
url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
metadata_filename <- "PBMC8_metadata.xlsx"
download.file(file.path(url, metadata_filename), destfile = file.path(DIR_RAW_DATA, metadata_filename))
md <- read_excel(file.path(DIR_RAW_DATA, metadata_filename))

## Load .fcs files
fcs_filename <- "PBMC8_fcs_files.zip"
download.file(file.path(url, fcs_filename), destfile = file.path(DIR_RAW_DATA, fcs_filename))
unzip(file.path(DIR_RAW_DATA, fcs_filename), exdir = DIR_RAW_DATA)
fcs_raw <- read.flowSet(file.path(DIR_RAW_DATA, md$file_name), 
                        transformation = FALSE, truncate_max_range = FALSE)

## Load panel file
panel_filename <- "PBMC8_panel.xlsx"
download.file(file.path(url, panel_filename), destfile = file.path(DIR_RAW_DATA, panel_filename))
panel <- read_excel(file.path(DIR_RAW_DATA, panel_filename))

# Replace problematic characters
panel$Antigen <- gsub("-", "_", panel$Antigen)
panel_fcs <- pData(parameters(fcs_raw[[1]]))
panel_fcs$desc <- gsub("-", "_", panel_fcs$desc)

# Lineage markers
(lineage_markers <- panel$Antigen[panel$Lineage == 1])

# Functional markers
(functional_markers <- panel$Antigen[panel$Functional == 1])

# Spot checks
stopifnot(all(lineage_markers %in% panel_fcs$desc))
stopifnot(all(functional_markers %in% panel_fcs$desc))

# -------------------
# Data transformation
# -------------------
## arcsinh transformation and column subsetting
fcs <- fsApply(fcs_raw, function(x, cofactor=5){
  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  expr <- asinh(expr[, colnames(expr)] / cofactor)
  exprs(x) <- expr
  x
})

# --------------------------------------------------------------------
# Cell population identification with FlowSOM and ConsensusClusterPlus
# --------------------------------------------------------------------
## FlowSOM clustering
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = lineage_markers)

## Metaclustering into 20 clusters with ConsensusClusterPlus
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

# ------------------------
# Cluster merging (manual)
# ------------------------

## Download manual merging scheme from Nowicka et al. (2017)
cluster_merging1_filename <- "PBMC8_cluster_merging1.xlsx"
download.file(file.path(url, cluster_merging1_filename),
              destfile = file.path(DIR_RAW_DATA, cluster_merging1_filename))
cluster_merging1 <- read_excel(file.path(DIR_RAW_DATA, cluster_merging1_filename))

## Convert to factor with merged clusters in correct order
levels_merged <- c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells", "CD8 T-cells", 
                   "DC", "NK cells", "monocytes", "surface-")
cluster_merging1$new_cluster <- factor(cluster_merging1$new_cluster, 
                                       levels = levels_merged)

## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

# ---------------------------
# Create SummarizedExperiment
# ---------------------------
## Generate row (cell) info table
sample_ids <- rep(gsub("\\.fcs$", "", gsub("^PBMC8_30min_", "", basename(md$file_name))), 
                  fsApply(fcs_raw, nrow))
group_ids <- factor(gsub("^patient[0-9+]_", "", sample_ids), 
                    levels = c("Reference", "BCR-XL"))
patient_ids <- factor(gsub("_.*$", "", sample_ids))
row_data <- data.frame(sample_id = sample_ids, group_id = group_ids, 
                       patient_id = patient_ids, cluster = cell_clustering1, 
                       population = cell_clustering1m, 
                       stringsAsFactors = FALSE)

## Extract abundances
expr <- fsApply(fcs, exprs)

bcrxl <- SummarizedExperiment(
  assays = list(exprs = expr),
  rowData = row_data,
  colData = data.frame(marker = colnames(expr),
                       lineage = colnames(expr) %in% lineage_markers,
                       functional = colnames(expr) %in% functional_markers,
                       row.names = colnames(expr),
                       stringsAsFactors = FALSE)
)

save(bcrxl, file = "bcrxl.rda")
