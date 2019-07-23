##########################################################################################
# Script to prepare benchmark dataset 'Weber-BCR-XL-sim-null'
# 
# See Weber et al. (2019), Supplementary Note 1 (paper introducing 'diffcyt' framework)
# for more details
# 
# The original 'BCR-XL' dataset is sourced from Bodenmiller et al. (2012), and was
# previously used for benchmark evaluations by Bruggner et al. (2014) (Citrus paper).
# 
# Raw data downloaded from Cytobank (experiment 15713)
# - see Citrus wiki (section 'PBMC Example 1'):
# https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - direct link to Cytobank repository:
# https://community.cytobank.org/cytobank/experiments/15713/download_files
# 
# Cell population labels are reproduced from Nowicka et al. (2017), where they were
# generated using a strategy of expert-guided manual merging of automatically generated
# clusters from the FlowSOM algorithm. Code to reproduce the cell population labels is
# available in the script 'cell_population_labels_BCR_XL.R'.
# 
# The 'Weber-BCR-XL-sim' dataset in this script is generated as follows:
# - select unstimulated reference samples from the main 'BCR-XL' dataset (8 individuals)
# - randomly split each unstimulated sample into two halves
# - in one half, replace B cells with equivalent number of B cells from the corresponding
# paired sample from BCR-XL stimulated condition
# 
# Methods are then evaluated by their ability to detect the known strong differential
# expression signal of the functional marker pS6 in B cells.
# 
# Lukas Weber, Jul 2019
##########################################################################################


# modified to create 'null' simulations: no true spike-in cells


# original version of this script available at: https://github.com/lmweber/diffcyt-evaluations

# note: random number generators were changed in R version 3.6.0; we use 'RNGversion()' to
# set random number generators to R version 3.5.3 for reproducibility (see
# https://cran.r-project.org/doc/manuals/r-devel/NEWS.html)

suppressPackageStartupMessages({
  library(flowCore)
  library(SummarizedExperiment)
})

RNGversion("3.5.3")



# -------------
# Download data
# -------------

# create temporary directories
DIR_TMP <- "tmp"
dir.create(file.path(DIR_TMP))
dir.create(file.path(DIR_TMP, "fcs_files"))
dir.create(file.path(DIR_TMP, "population_IDs"))

# download from 'imlspenticton' server
URL <- "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData"
DIR <- "Bodenmiller_BCR_XL"

# download .fcs files
fcs_filename <- "Bodenmiller_BCR_XL_fcs_files.zip"
download.file(file.path(URL, DIR, fcs_filename), destfile = file.path(DIR_TMP, "fcs_files", fcs_filename))
unzip(file.path(DIR_TMP, "fcs_files", fcs_filename), exdir = file.path(DIR_TMP, "fcs_files"))

# download population IDs
pop_filename <- "Bodenmiller_BCR_XL_population_IDs.zip"
download.file(file.path(URL, DIR, pop_filename), destfile = file.path(DIR_TMP, "population_IDs", pop_filename))
unzip(file.path(DIR_TMP, "population_IDs", pop_filename), exdir = file.path(DIR_TMP, "population_IDs"))



# ---------
# Filenames
# ---------

DIR_FCS <- file.path(DIR_TMP, "fcs_files")
DIR_LABELS <- file.path(DIR_TMP, "population_IDs")

# .fcs files
files_fcs <- list.files(DIR_FCS, pattern = "\\.fcs$", full.names = TRUE)
files_BCRXL <- files_fcs[grep("patient[1-8]_BCR-XL\\.fcs$", files_fcs)]
files_ref <- files_fcs[grep("patient[1-8]_Reference\\.fcs$", files_fcs)]
files_fcs_all <- c(files_BCRXL, files_ref)

# cell population labels
files_labels <- list.files(DIR_LABELS, pattern = "\\.csv$", full.names = TRUE)
files_labels_BCRXL <- files_labels[grep("patient[1-8]_BCR-XL\\.csv$", files_labels)]
files_labels_ref <- files_labels[grep("patient[1-8]_Reference\\.csv$", files_labels)]
files_labels_all <- c(files_labels_BCRXL, files_labels_ref)



# ---------
# Load data
# ---------

data <- lapply(files_fcs_all, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

# sample IDs
sample_IDs <- gsub("\\.fcs$", "", gsub("^PBMC8_30min_", "", basename(files_fcs_all)))
sample_IDs

names(data) <- sample_IDs

# conditions
conditions <- as.factor(gsub("^patient[0-9+]_", "", sample_IDs))
conditions

# patient IDs
patient_IDs <- as.factor(gsub("_.*$", "", sample_IDs))
patient_IDs

# cell population labels
labels <- lapply(files_labels_all, read.csv)

stopifnot(all(sapply(data, nrow) == sapply(labels, nrow)))

names(labels) <- names(data)



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ---------------------
# Randomized replicates
# ---------------------

# use different random seed for each replicate

n_replicates <- 3

data_replicates <- labels_replicates <- vector("list", n_replicates)
names(data_replicates) <- names(labels_replicates) <- paste0("seed", 1:n_replicates)


for (r in 1:n_replicates) {
  
  data_replicates[[r]] <- labels_replicates[[r]] <- vector("list", 2)
  names(data_replicates[[r]]) <- names(labels_replicates[[r]]) <- c("null1", "null2")
  
  
  
  # ---------------------------------------------------------------
  # Split each reference sample into two halves: 'base' and 'spike'
  # ---------------------------------------------------------------
  
  data_ref <- data[conditions == "Reference"]
  labels_ref <- labels[conditions == "Reference"]
  
  names(data_ref) <- gsub("_Reference$", "", names(data_ref))
  names(labels_ref) <- gsub("_Reference$", "", names(labels_ref))
  
  stopifnot(all(sapply(data_ref, nrow) == sapply(labels_ref, nrow)))
  
  n_cells_ref <- sapply(data_ref, nrow)
  
  # modified random seed for each replicate
  set.seed(10000 + 100 * r)
  
  # generate random indices
  inds <- lapply(n_cells_ref, function(n) {
    i_null_1 <- sort(sample(seq_len(n), floor(n / 2)))
    i_null_2 <- setdiff(seq_len(n), i_null_1)
    list(null_1 = i_null_1, null_2 = i_null_2)
  })
  
  inds_null_1 <- lapply(inds, function(l) l[[1]])
  inds_null_2 <- lapply(inds, function(l) l[[2]])
  
  # subset data
  data_null_1 <- mapply(function(d, i) d[i, , drop = FALSE], data_ref, inds_null_1, SIMPLIFY = FALSE)
  data_null_2 <- mapply(function(d, i) d[i, , drop = FALSE], data_ref, inds_null_2, SIMPLIFY = FALSE)
  
  # subset labels
  labels_null_1 <- mapply(function(d, i) d[i, , drop = FALSE], labels_ref, inds_null_1, SIMPLIFY = FALSE)
  labels_null_2 <- mapply(function(d, i) d[i, , drop = FALSE], labels_ref, inds_null_2, SIMPLIFY = FALSE)
  
  stopifnot(all(sapply(data_null_1, nrow) == sapply(labels_null_1, nrow)))
  stopifnot(all(sapply(data_null_2, nrow) == sapply(labels_null_2, nrow)))
  
  # store data
  data_replicates[[r]][["null1"]] <- data_null_1
  data_replicates[[r]][["null2"]] <- data_null_2
  
  labels_replicates[[r]][["null1"]] <- labels_null_1
  labels_replicates[[r]][["null2"]] <- labels_null_2
  
}



# ---------------
# Create metadata
# ---------------

data_replicates_combined <- lapply(data_replicates, function(r) c(r[["null1"]], r[["null2"]]))
labels_replicates_combined <- lapply(labels_replicates, function(r) c(r[["null1"]], r[["null2"]]))

for (r in 1:n_replicates) {
  stopifnot(all(sapply(data_replicates_combined[[r]], nrow) == sapply(labels_replicates_combined[[r]], nrow)))
}

conditions_combined <- c(rep("null1", length(data_replicates[[1]][["null1"]])), rep("null2", length(data_replicates[[1]][["null2"]])))


# sample information

patient_id <- as.factor(names(data_replicates_combined[[1]]))
patient_id

group_id <- factor(conditions_combined, levels = c("null1", "null2"))
group_id

sample_id <- paste(patient_id, group_id, sep = "_")
sample_id <- factor(sample_id, levels = sample_id)
sample_id

experiment_info <- data.frame(group_id, patient_id, sample_id, stringsAsFactors = FALSE)
experiment_info


for (r in 1:n_replicates) {
  names(data_replicates_combined[[r]]) <- sample_id
  names(labels_replicates_combined[[r]]) <- sample_id
}


# marker information

# indices of all marker columns, lineage markers, and functional markers
# (10 lineage markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

stopifnot(all(sapply(seq_along(data), function(i) all(colnames(data[[i]]) == colnames(data[[1]])))))

# channel and marker names
channel_name <- as.character(colnames(data[[1]]))
marker_name <- gsub("\\(.*$", "", channel_name)

# marker classes
marker_class <- rep("none", length(marker_name))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("none", "type", "state"))

marker_info <- data.frame(channel_name, marker_name, marker_class, stringsAsFactors = FALSE)
marker_info



# -----------------------------------
# Create SummarizedExperiment objects
# -----------------------------------

# create a separate object for each replicate (since row data is different for each replicate)


d_SE_list <- vector("list", length(data_replicates_combined))
names(d_SE_list) <- names(data_replicates_combined)

for (r in 1:length(d_SE_list)) {
  
  # set up row data
  n_cells <- sapply(data_replicates_combined[[r]], dim)[1, ]
  stopifnot(length(n_cells) == nrow(experiment_info))
  stopifnot(all(names(n_cells) == experiment_info$sample_id))
  
  row_data <- as.data.frame(lapply(experiment_info, function(col) {
    as.factor(rep(col, n_cells))
  }))
  
  stopifnot(nrow(row_data) == sum(n_cells))
  
  # add population IDs
  population_id <- do.call("rbind", labels_replicates_combined[[r]])
  rownames(population_id) <- NULL
  colnames(population_id) <- "population_id"
  
  stopifnot(nrow(population_id) == sum(n_cells))
  
  row_data <- cbind(row_data, population_id)
  
  # add column indicating B cells
  row_data$B_cell <- row_data$population_id %in% c("B-cells IgM-", "B-cells IgM+")
  
  # add column indicating spike-in cells (all B cells in 'spike' samples)
  row_data$spikein <- FALSE
  
  
  # set up column data
  col_data <- marker_info
  
  
  # set up expression data
  d_exprs <- do.call("rbind", data_replicates_combined[[r]])
  
  stopifnot(nrow(d_exprs) == nrow(row_data))
  stopifnot(ncol(d_exprs) == nrow(col_data))
  stopifnot(nrow(d_exprs) == sum(n_cells))
  stopifnot(all(colnames(d_exprs) == marker_info$channel_name))
  
  
  # create SummarizedExperiment object
  d_SE_list[[r]] <- SummarizedExperiment(
    assays = list(exprs = d_exprs), 
    rowData = row_data, 
    colData = col_data, 
    metadata = list(experiment_info = experiment_info, n_cells = n_cells)
  )
}



# ----------------------
# Create flowSet objects
# ----------------------

# note: row data is stored as additional columns of data in the expression matrices;
# additional information from row data and column data (e.g. marker classes) is stored in
# 'description' slot


d_flowSet_list <- vector("list", length(data_replicates_combined))
names(d_flowSet_list) <- names(data_replicates_combined)

for (r in 1:length(d_flowSet_list)) {
  
  # extract data
  row_data <- rowData(d_SE_list[[r]])
  col_data <- colData(d_SE_list[[r]])
  d_exprs <- assay(d_SE_list[[r]])
  meta_data <- metadata(d_SE_list[[r]])
  
  # create tables to identify row data values when converted to numeric
  stopifnot(all(colnames(row_data) == c("group_id", "patient_id", "sample_id", "population_id", "B_cell", "spikein")))
  group_info <- data.frame(
    group_id = seq_len(nlevels(row_data$group_id)), 
    group_name = levels(row_data$group_id), 
    stringsAsFactors = FALSE
  )
  patient_info <- data.frame(
    patient_id = seq_len(nlevels(row_data$patient_id)), 
    patient_name = levels(row_data$patient_id), 
    stringsAsFactors = FALSE
  )
  sample_info <- data.frame(
    sample_id = seq_len(nlevels(row_data$sample_id)), 
    sample_name = levels(row_data$sample_id), 
    stringsAsFactors = FALSE
  )
  population_info <- data.frame(
    population_id = seq_len(nlevels(row_data$population_id)), 
    population_name = levels(row_data$population_id), 
    stringsAsFactors = FALSE
  )
  B_cell_info <- data.frame(
    B_cell = c(0, 1), 
    B_cell_status = c("FALSE", "TRUE"), 
    stringsAsFactors = FALSE
  )
  spikein_info <- data.frame(
    spikein = c(0, 1), 
    spikein_status = c("FALSE", "TRUE"), 
    stringsAsFactors = FALSE
  )
  
  # create extra columns of data from row data
  row_data_fs <- do.call("cbind", lapply(row_data, as.numeric))
  stopifnot(all(colnames(row_data_fs) == colnames(row_data)))
  stopifnot(nrow(row_data_fs) == nrow(d_exprs))
  
  # create marker info
  marker_info <- as.data.frame(col_data, row.names = seq_len(nrow(col_data)))
  
  # create flowFrame
  ff <- flowFrame(cbind(d_exprs, row_data_fs))
  stopifnot(all(colnames(ff) == c(colnames(d_exprs), colnames(row_data_fs))))
  # include both channel and marker names in 'pData(parameters(.))'
  stopifnot(length(c(marker_info$marker_name, colnames(row_data_fs))) == nrow(pData(parameters(ff))))
  pData(parameters(ff))$desc <- c(marker_info$marker_name, colnames(row_data_fs))
  
  # include additional information in 'description' slot
  description(ff)$GROUP_INFO <- group_info
  description(ff)$PATIENT_INFO <- patient_info
  description(ff)$SAMPLE_INFO <- sample_info
  description(ff)$POPULATION_INFO <- population_info
  description(ff)$B_CELL_INFO <- B_cell_info
  description(ff)$SPIKEIN_INFO <- spikein_info
  # experiment information and marker information
  description(ff)$EXPERIMENT_INFO <- meta_data$experiment_info
  description(ff)$MARKER_INFO <- marker_info
  # simulation replicate (seed)
  description(ff)$REPLICATE <- names(d_flowSet_list)[r]
  
  # create flowSet
  d_flowSet_list[[r]] <- flowSet(ff)
}



# ------------
# Save objects
# ------------

stopifnot(all(names(d_SE_list) == names(d_flowSet_list)))

reps <- paste0("rep", 1:3)

filenames_SE <- paste0("Weber_BCR_XL_sim_null_", reps, "_SE.rda")
filenames_flowSet <- paste0("Weber_BCR_XL_sim_null_", reps, "_flowSet.rda")

for (r in 1:3) {
  d_SE <- d_SE_list[[r]]
  d_flowSet <- d_flowSet_list[[r]]
  
  save(d_SE, file = filenames_SE[r])
  save(d_flowSet, file = filenames_flowSet[r])
}



