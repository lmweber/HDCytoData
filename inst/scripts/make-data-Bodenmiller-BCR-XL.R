##########################################################################################
# R script to prepare benchmark dataset Bodenmiller_BCR_XL
# 
# This is a mass cytometry (CyTOF) dataset from Bodenmiller et al. (2012), consisting of
# paired samples of peripheral blood mononuclear cells (PBMCs) from healthy individuals,
# where one sample from each pair was stimulated with B cell receptor / Fc receptor
# cross-linker (BCR-XL). There are 8 paired samples (i.e. 16 samples in total). The
# dataset contains expression levels of 24 protein markers (10 surface markers and 14
# intracellular signaling markers). Cell population labels are available from Nowicka et
# al. (2017). This dataset contains strong differential expression signals for several
# signaling markers in several cell populations. In particular, one of the strongest
# effects is differential expression of phosphorylated S6 (pS6) in B cells.
# 
# This R script loads the data, adds row and column metadata (including cell population
# labels), and exports it in SummarizedExperiment and flowSet formats.
# 
# Source:
# - Original source: Bodenmiller et al. (2012)
# - Link to paper: https://www.ncbi.nlm.nih.gov/pubmed/22902532
# - Link to raw data: https://community.cytobank.org/cytobank/experiments/15713/download_files
# - Additional information: https://github.com/nolanlab/citrus/wiki/PBMC-Example-1
# - Cell population labels from: Nowicka et al. (2017), v2: https://f1000research.com/articles/6-748/v2
# 
# Lukas Weber, Jan 2019
##########################################################################################


suppressPackageStartupMessages({
  library(flowCore)
  library(SummarizedExperiment)
})



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

files_load_fcs <- list.files(file.path(DIR_TMP, "fcs_files"), pattern = "\\.fcs$", full.names = TRUE)

# download population IDs
pop_filename <- "Bodenmiller_BCR_XL_population_IDs.zip"
download.file(file.path(URL, DIR, pop_filename), destfile = file.path(DIR_TMP, "population_IDs", pop_filename))
unzip(file.path(DIR_TMP, "population_IDs", pop_filename), exdir = file.path(DIR_TMP, "population_IDs"))

files_load_pop <- list.files(file.path(DIR_TMP, "population_IDs"), pattern = "\\.csv$", full.names = TRUE)



# ---------
# Load data
# ---------

# load .fcs files
data_flowSet <- read.flowSet(files_load_fcs, transformation = FALSE, truncate_max_range = FALSE)

# load population IDs
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
# (10 lineage markers / 14 functional markers; see Bruggner et al. 2014, Table 1)
cols_markers <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
cols_lineage <- c(3:4, 9, 11, 12, 14, 21, 29, 31, 33)
cols_func <- setdiff(cols_markers, cols_lineage)

# channel and marker names
channel_name <- as.character(pData(parameters(data_flowSet[[1]]))$name)
marker_name <- as.character(pData(parameters(data_flowSet[[1]]))$desc)
# original column names
col_names <- colnames(data_flowSet)
# check
all(channel_name == col_names)
all(marker_name == gsub("\\(.*$", "", channel_name))

# marker classes
marker_class <- rep("none", length(marker_name))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, levels = c("none", "type", "state"))

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


# set up column data
col_data <- marker_info


# set up expression data
d_exprs <- fsApply(data_flowSet, exprs)
stopifnot(all(colnames(d_exprs) == channel_name))

# use marker names as column names (for SummarizedExperiment)
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



# ---------------------
# Create flowSet object
# ---------------------

# note: row data (e.g. population IDs) is stored as additional columns of data in the
# expression matrices; additional information from row data and column data (e.g. marker
# classes, cell population names) is stored in 'description' slot

# table of cell population information
population_info <- data.frame(
  population_id = seq_len(nlevels(row_data$population_id)), 
  population_name = levels(row_data$population_id), 
  stringsAsFactors = FALSE
)

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

# add extra columns of data and create new flowSet object
data_flowSet_list <- as(data_flowSet, "list")

d_flowFrames_list <- mapply(function(d, extra_cols) {
  e <- exprs(d)
  stopifnot(nrow(e) == nrow(extra_cols))
  # use original column names (for flowSet)
  colnames(e) <- col_names
  # combine and create flowFrame
  flowFrame(cbind(e, extra_cols))
}, data_flowSet_list, row_data_fs_list)

d_flowSet <- flowSet(d_flowFrames_list)

# include additional information in 'description' slot
for (i in seq_along(d_flowSet)) {
  # filename and sample information
  description(d_flowSet[[i]])$FILENAME <- identifier(d_flowSet[[i]])
  description(d_flowSet[[i]])$GROUP_ID <- group_id[i]
  description(d_flowSet[[i]])$PATIENT_ID <- patient_id[i]
  description(d_flowSet[[i]])$SAMPLE_ID <- sample_id[i]
  # data frames of experiment information and marker information
  description(d_flowSet[[i]])$EXPERIMENT_INFO <- experiment_info
  description(d_flowSet[[i]])$MARKER_INFO <- marker_info
  # data frame of cell population information
  description(d_flowSet[[i]])$POPULATION_INFO <- population_info
}



# ------------
# Save objects
# ------------

save(d_SE, file = "Bodenmiller_BCR_XL_SE.rda")
save(d_flowSet, file = "Bodenmiller_BCR_XL_flowSet.rda")



