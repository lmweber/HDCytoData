##########################################################################################
# R script to prepare benchmark dataset Krieg_Anti_PD_1
# 
# This is a mass cytometry (CyTOF) dataset from Krieg et al. (2018), who used mass
# cytometry to characterize immune cell subsets in peripheral blood from melanoma skin
# cancer patients treated with anti-PD-1 immunotherapy. This study found that the
# frequency of CD14+CD16−HLA-DRhi monocytes in baseline samples (taken from patients prior
# to treatment) was a strong predictor of survival in response to immunotherapy treatment.
# In particular, the frequency of a small subpopulation of
# CD14+CD33+HLA-DRhiICAM-1+CD64+CD141+CD86+CD11c+CD38+PD-L1+CD11b+ monocytes in baseline
# samples was strongly associated with responder status following immunotherapy treatment.
# Note that this dataset contains a strong batch effect, due to sample acquisition on two
# different days (Krieg et al., 2018).
# 
# This dataset was subsequently used for benchmarking purposes in our paper:
# Weber et al. (2018): https://www.biorxiv.org/content/early/2018/11/22/349738
# 
# This R script loads the data, adds row and column metadata, and exports it in
# SummarizedExperiment and flowSet formats.
# 
# Source: Krieg et al. (2018): https://www.ncbi.nlm.nih.gov/pubmed/29309059
# 
# Link to raw data: http://flowrepository.org/id/FR-FCM-ZY34
# 
# Additional information: Weber et al. (2018),
# https://www.biorxiv.org/content/early/2018/11/22/349738 (see Supplementary Note 3:
# Benchmark datasets)
# 
# Data files also available from: http://flowrepository.org/id/FR-FCM-ZYL8
# 
# Lukas Weber, Jan 2019
##########################################################################################


suppressPackageStartupMessages({
  library(flowCore)
  library(SummarizedExperiment)
  library(readxl)
})



# --------------
# Download files
# --------------

# create temporary directories
DIR_FILES <- "tmp"
dir.create(file.path(DIR_FILES))

# download from 'imlspenticton' server
URL <- "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData"
DIR <- "Krieg_Anti_PD_1"

# download files
filename <- "Krieg_Anti_PD_1.zip"
download.file(file.path(URL, DIR, filename), destfile = file.path(DIR_FILES, filename))
unzip(file.path(DIR_FILES, filename), exdir = file.path(DIR_FILES))



# -----------------------------------
# Filenames and metadata spreadsheets
# -----------------------------------

# note: data from batches '23' and '29'

fn_metadata_23 <- file.path(DIR_FILES, "CK_metadata", "metadata_23_03all.xlsx")
fn_metadata_29 <- file.path(DIR_FILES, "CK_metadata", "metadata_29_03all3.xlsx")

path_23 <- file.path(DIR_FILES, "CK_2016-06-23_03all", "010_cleanfcs")
path_29 <- file.path(DIR_FILES, "CK_2016-06-29_03all3", "010_cleanfcs")

fn_panel <- file.path(DIR_FILES, "CK_panels", "panel3_v3.xlsx")

# load metadata spreadsheets for each data set ("data 23" and "data 29")
metadata_23 <- read_excel(fn_metadata_23)
metadata_29 <- read_excel(fn_metadata_29)

ix_keep <- 6:15

paths <- c(rep(path_23, length(ix_keep)), rep(path_29, length(ix_keep)))

files <- c(metadata_23$filename[ix_keep], metadata_29$filename[ix_keep])



# ------------------
# Sample information
# ------------------

# group IDs
group_id <- factor(
  gsub("^base_", "", c(metadata_23$condition[ix_keep], metadata_29$condition[ix_keep])), 
  levels = c("NR", "R")
)
group_id

# batch IDs
batch_id <- factor(
  c(rep("batch23", length(ix_keep)), rep("batch29", length(ix_keep))), 
  levels = c("batch23", "batch29")
)
batch_id

# sample IDs
sample_id <- factor(
  gsub("^base_", "", c(metadata_23$shortname[ix_keep], metadata_29$shortname[ix_keep])), 
  levels = c(paste0("NR", 1:9), paste0("R", 1:11))
)
sample_id

experiment_info <- data.frame(group_id, batch_id, sample_id, stringsAsFactors = FALSE)
experiment_info



# ---------
# Load data
# ---------

stopifnot(length(paths) == length(files))

fn <- file.path(paths, files)

d_input <- lapply(fn, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

names(d_input) <- files

# check column names match (note: exclude columns 58 and 59, which contain "beadDist" and "Time")
ix <- 1:57
check_cols <- lapply(d_input, function(d) pData(parameters(d))$name)
all(sapply(check_cols, function(ch) all(ch[ix] == check_cols[[1]][ix])))



# ------------------
# Marker information
# ------------------

# note: not including CD45 as a 'cell type' marker (since almost all cells are high in
# CD45, so does not help distinguish populations)

# load panel details spreadsheet
panel <- as.data.frame(read_excel(fn_panel))
panel

# replace NAs
panel$Antigen[is.na(panel$Antigen)] <- panel$fcs_colname[is.na(panel$Antigen)]
panel

# subset and rearrange data; so order of columns matches rows in panel

markers_ix <- match(panel$fcs_colname, pData(parameters(d_input[[1]]))$name)

stopifnot(length(markers_ix) == nrow(panel))

# create list of matrices
d_exprs_list <- lapply(d_input, function(d) {
  e <- exprs(d)[, markers_ix]
  e
})

# marker information
# (note: all surface markers in this panel)

is_marker <- as.logical(panel$transform)

# marker classes (note: CD45 is classified as 'none')
marker_class <- rep("none", nrow(panel))
marker_class[is_marker] <- "type"
marker_class <- factor(marker_class, levels = c("none", "type", "state"))

# channel names and marker names
channel_name <- as.character(panel$fcs_colname)
marker_name <- as.character(panel$Antigen)

marker_info <- data.frame(channel_name, marker_name, marker_class, stringsAsFactors = FALSE)
marker_info



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_FILES, recursive = TRUE)



# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
n_cells <- sapply(d_exprs_list, nrow)
stopifnot(all(experiment_info$sample_id == gsub(".*_([NR]+[0-9]+)\\.fcs$", "\\1", names(n_cells))))
names(n_cells) <- experiment_info$sample_id

stopifnot(nrow(experiment_info) == length(n_cells))

row_data <- data.frame(
  group_id = rep(experiment_info$group_id, n_cells), 
  batch_id = rep(experiment_info$batch_id, n_cells), 
  sample_id = rep(experiment_info$sample_id, n_cells), 
  stringsAsFactors = FALSE
)

stopifnot(nrow(row_data) == sum(n_cells))


# set up column data
col_data <- marker_info


# set up expression table
d_exprs <- do.call("rbind", d_exprs_list)

stopifnot(nrow(d_exprs) == sum(n_cells))
stopifnot(ncol(d_exprs) == nrow(marker_info))

# use marker names as column names (for SummarizedExperiment)
colnames(d_exprs) <- marker_info$marker_name


# create SummarizedExperiment object
d_SE <- SummarizedExperiment(
  assays = list(exprs = d_exprs), 
  rowData = row_data, 
  colData = col_data, 
  metadata = list(experiment_info = experiment_info, n_cells = n_cells)
)



# ---------------------
# Create flowSet object
# ---------------------

# note: row data is stored as additional columns of data in the expression matrices;
# additional information from row data and column data (e.g. marker classes) is stored in
# 'description' slot

# create list of extra columns of data for each sample
row_data_fs <- do.call("cbind", lapply(row_data, as.numeric))
stopifnot(nrow(row_data_fs) == nrow(row_data))

row_data_fs_list <- lapply(split(row_data_fs, row_data$sample_id), matrix, ncol = ncol(row_data_fs))
# rearrange samples into correct order
row_data_fs_list <- row_data_fs_list[sample_id]
stopifnot(all(names(row_data_fs_list) == sample_id))

# replace column names
row_data_fs_list <- lapply(row_data_fs_list, function(d) {
  colnames(d) <- colnames(row_data_fs)
  d
})

# add extra columns of data and create new flowSet object
d_flowFrames_list <- mapply(function(e, extra_cols) {
  stopifnot(nrow(e) == nrow(extra_cols))
  # use original column names (for flowSet)
  stopifnot(ncol(e) == length(channel_name))
  colnames(e) <- channel_name
  # combine and create flowFrame
  ff <- flowFrame(cbind(e, extra_cols))
  # include both channel and marker names in 'pData(parameters(.))'
  stopifnot(length(c(marker_info$marker_name, colnames(extra_cols))) == nrow(pData(parameters(ff))))
  pData(parameters(ff))$desc <- c(marker_info$marker_name, colnames(extra_cols))
  ff
}, d_exprs_list, row_data_fs_list)

d_flowSet <- flowSet(d_flowFrames_list)

# include additional information in 'description' slot
for (i in seq_along(d_flowSet)) {
  # filename and sample information
  stopifnot(identifier(d_flowSet[[i]]) == files[i])
  description(d_flowSet[[i]])$FILENAME <- identifier(d_flowSet[[i]])
  description(d_flowSet[[i]])$GROUP_ID <- group_id[i]
  description(d_flowSet[[i]])$BATCH_ID <- batch_id[i]
  description(d_flowSet[[i]])$SAMPLE_ID <- sample_id[i]
  # data frames of experiment information and marker information
  description(d_flowSet[[i]])$EXPERIMENT_INFO <- experiment_info
  description(d_flowSet[[i]])$MARKER_INFO <- marker_info
}



# ------------
# Save objects
# ------------

save(d_SE, file = "Krieg_Anti_PD_1_SE.rda")
save(d_flowSet, file = "Krieg_Anti_PD_1_flowSet.rda")



