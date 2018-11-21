##########################################################################################
# Dataset: Bodenmiller_BCR_XL
# 
# Script to load data from .fcs files, add row and column metadata (including reference
# population labels), and export in 'SummarizedExperiment' and 'flowSet' formats
##########################################################################################


suppressPackageStartupMessages({
  library(flowCore)
  library(SummarizedExperiment)
})



# ----------------------
# Download and load data
# ----------------------

# create temporary directories
DIR_TMP <- "tmp"
dir.create(file.path(DIR_TMP), showWarnings = FALSE)
dir.create(file.path(DIR_TMP, "fcs_files"), showWarnings = FALSE)
dir.create(file.path(DIR_TMP, "population_IDs"), showWarnings = FALSE)

# download from 'imlspenticton' server
URL <- "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData"

# load .fcs files
fcs_filename <- "Bodenmiller_BCR_XL_fcs_files.zip"
download.file(file.path(URL, fcs_filename), destfile = file.path(DIR_TMP, "fcs_files", fcs_filename))
unzip(file.path(DIR_TMP, "fcs_files", fcs_filename), exdir = file.path(DIR_TMP, "fcs_files"))

files_load_fcs <- list.files(file.path(DIR_TMP, "fcs_files"), pattern = "\\.fcs$", full.names = TRUE)

data_flowSet <- read.flowSet(files_load_fcs, transformation = FALSE, truncate_max_range = FALSE)

# load population IDs
pop_filename <- "Bodenmiller_BCR_XL_population_IDs.zip"
download.file(file.path(URL, pop_filename), destfile = file.path(DIR_TMP, "population_IDs", pop_filename))
unzip(file.path(DIR_TMP, "population_IDs", pop_filename), exdir = file.path(DIR_TMP, "population_IDs"))

files_load_pop <- list.files(file.path(DIR_TMP, "population_IDs"), pattern = "\\.csv$", full.names = TRUE)

data_population_IDs <- lapply(files_load_pop, read.csv)

# check numbers of cells match
stopifnot(all(sapply(seq_along(files_load_fcs), function(i) {
  nrow(data_population_IDs[[i]]) == nrow(data_flowSet[[i]])
})))



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ---------------
# Create metadata
# ---------------

# sample information

sample_id <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", basename(files_load_fcs)))
sample_id

group_id <- factor(gsub("^.*_", "", sample_id), levels = c("Reference", "BCR-XL"))
group_id

patient_id <- factor(gsub("_.*$", "", sample_id))
patient_id

experiment_info <- data.frame(group_id, patient_id, sample_id, stringsAsFactors = FALSE)
experiment_info


# marker information

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

channel_name <- colnames(data_flowSet)

marker_name <- gsub("\\(.*$", "", channel_name)

marker_class <- rep("none", length(marker_name))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("type", "state", "none"))

marker_info <- data.frame(channel_name, marker_name, marker_class, stringsAsFactors = FALSE)
marker_info



# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
n_cells <- sapply(as(data_flowSet, "list"), nrow)
names(n_cells) <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", names(n_cells)))

stopifnot(all(names(n_cells) == sample_id))

row_data <- as.data.frame(lapply(experiment_info, function(col) {
  as.factor(rep(col, n_cells))
}))

stopifnot(nrow(row_data) == sum(n_cells))

# add population IDs
population_id <- do.call("rbind", data_population_IDs)
colnames(population_id) <- "population_id"
stopifnot(nrow(population_id) == sum(n_cells))

row_data <- cbind(row_data, population_id)

# column data
col_data <- marker_info

# collapse expression data into a single table of values
d_exprs <- fsApply(data_flowSet, exprs)
stopifnot(all(colnames(d_exprs) == channel_name))

colnames(d_exprs) <- marker_name

stopifnot(nrow(d_exprs) == sum(n_cells), 
          ncol(d_exprs) == nrow(col_data))

# create SummarizedExperiment object
d_SE <- SummarizedExperiment(
  assays = list(exprs = d_exprs), 
  rowData = row_data, 
  colData = col_data, 
  metadata = list(experiment_info = experiment_info, n_cells = n_cells)
)



# --------------
# Create flowSet
# --------------

# note: sample information is stored as additional columns of data in the expression value
# matrices; additional marker information (channel names and marker classes) cannot be
# included, since marker information is stored in column names only

# create list of extra columns of data for each sample
row_data_fs <- do.call("cbind", lapply(row_data, as.numeric))
stopifnot(nrow(row_data_fs) == nrow(row_data))

row_data_fs_list <- lapply(split(row_data_fs, row_data$sample_id), matrix, ncol = ncol(row_data_fs))
stopifnot(all(names(row_data_fs_list) == sample_id))

# replace column names
row_data_fs_list <- lapply(row_data_fs_list, function(d) {
  colnames(d) <- colnames(row_data_fs)
  d
})

# create new flowSet object and add extra columns of data
data_flowSet_list <- as(data_flowSet, "list")

d_flowFrames_list <- mapply(function(d, extra_cols) {
  e <- exprs(d)
  stopifnot(nrow(e) == nrow(extra_cols))
  # combine and create flowFrame
  flowFrame(cbind(e, extra_cols))
}, data_flowSet_list, row_data_fs_list)

d_flowSet <- flowSet(d_flowFrames_list)



# ------------
# Save objects
# ------------

save(d_SE, file = "Bodenmiller_BCR_XL_SE.rda")
save(d_flowSet, file = "Bodenmiller_BCR_XL_flowSet.rda")



