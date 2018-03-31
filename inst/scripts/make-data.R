##########################################################################################
# Data set: Bodenmiller_BCR_XL
# 
# Script to load data from .fcs files, add meta-data and population labels, and export in
# 'SummarizedExperiment' and 'flowSet' formats
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
URL <- "http://imlspenticton.uzh.ch/robinson_lab/cytofData/"

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


# ------------------
# Meta-data: samples
# ------------------

# get meta-data from filenames

# check samples are in correct order
stopifnot(all(pData(data_flowSet)$name == basename(files_load_fcs)))

# sample information
sample_IDs <- gsub("^PBMC8_30min_", "", gsub("\\.fcs$", "", basename(files_load_fcs)))
group_IDs <- factor(gsub("^patient[0-9]+_", "", sample_IDs), levels = c("BCR-XL", "Reference"))
patient_IDs <- factor(gsub("_.*$", "", sample_IDs))

sample_info <- data.frame(group_IDs, patient_IDs, sample_IDs, stringsAsFactors = FALSE)


# ------------------
# Meta-data: markers
# ------------------

# additional meta-data to identify 'cell type' and 'cell state' markers

# indices of all marker columns, lineage markers, and functional markers
# (10 surface markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

# clean marker names
marker_names <- colnames(data_flowSet)
marker_names <- gsub("\\(.*$", "", marker_names)

# marker information
is_marker <- is_type_marker <- is_state_marker <- rep(FALSE, length(marker_names))
is_marker[cols_markers] <- TRUE
is_type_marker[cols_lineage] <- TRUE
is_state_marker[cols_func] <- TRUE

marker_info <- data.frame(marker_name = marker_names, 
                          is_marker, is_type_marker, is_state_marker, 
                          stringsAsFactors = FALSE)


# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
n_cells <- sapply(as(data_flowSet, "list"), nrow)
stopifnot(all(names(n_cells) == pData(data_flowSet)$name))

row_data <- as.data.frame(lapply(sample_info, function(col) {
  as.factor(rep(col, n_cells))
}))

colnames(row_data) <- gsub("_IDs$", "", colnames(row_data))

stopifnot(nrow(row_data) == sum(n_cells))

# add population IDs
population <- do.call("rbind", data_population_IDs)
stopifnot(nrow(population) == sum(n_cells))

row_data <- cbind(row_data, population)

# set up column data
col_data <- marker_info

# collapse data into a single table of expression values
d_exprs <- fsApply(data_flowSet, exprs)
colnames(d_exprs) <- NULL

stopifnot(nrow(d_exprs) == sum(n_cells), 
          ncol(d_exprs) == nrow(col_data))

# create SummarizedExperiment object
d_SE <- SummarizedExperiment(
  d_exprs, 
  rowData = row_data, 
  colData = col_data, 
  metadata = list(sample_info = sample_info)
)


# --------------
# Create flowSet
# --------------

# add sample information as additional columns in the expression data

# (note: additional marker information cannot be included, since marker information is
# stored in column names only)

# create list of extra columns for each sample
row_data_FS <- do.call("cbind", lapply(row_data, as.numeric))
stopifnot(nrow(row_data_FS) == nrow(row_data))

row_data_FS_list <- lapply(split(row_data_FS, row_data$sample), matrix, ncol = ncol(row_data_FS))

# replace column names
row_data_FS_list <- lapply(row_data_FS_list, function(d) {
  colnames(d) <- colnames(row_data_FS)
  d
})

# create new flowSet object and add extra columns
data_flowSet_list <- as(data_flowSet, "list")

d_FF <- mapply(function(d, extra_cols) {
  e <- exprs(d)
  stopifnot(nrow(e) == nrow(extra_cols))
  # clean column (marker) names
  colnames(e) <- gsub("\\(.*$", "", colnames(e))
  # combine and create flowFrame
  flowFrame(cbind(e, extra_cols))
}, data_flowSet_list, row_data_FS_list)

d_FS <- flowSet(d_FF)


# ------------
# Save objects
# ------------

save(d_SE, file = "Bodenmiller_BCR_XL_SE.rda")
save(d_FS, file = "Bodenmiller_BCR_XL_flowSet.rda")


