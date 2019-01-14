##########################################################################################
# R script to prepare benchmark dataset Levine_13dim
# 
# This is a 13-dimensional mass cytometry (CyTOF) dataset, consisting of expression levels
# of 13 surface protein markers. Cell population labels are available for 24 manually
# gated populations. Cells are healthy human bone marrow mononuclear cells (BMMCs), from 1
# patient.
#
# This R script loads the data, adds manually gated cell population labels, and exports it
# in SummarizedExperiment and flowSet formats.
#
# Source: "benchmark data set 1" in the following paper:
# Levine et al. (2015), "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like
# Cells that Correlate with Prognosis", Cell, 162, 184-197.
#
# Link to paper: https://www.ncbi.nlm.nih.gov/pubmed/26095251
# Link to raw data: https://www.cytobank.org/cytobank/experiments/46259 (download the .fcs
# files with Actions -> Export -> Download Files -> All FCS Files)
# 
# Lukas Weber, Jan 2019
##########################################################################################


# original version of this script can be found at:
# https://github.com/lmweber/cytometry-clustering-comparison


suppressPackageStartupMessages({
  library(flowCore)
  library(SummarizedExperiment)
  library(magrittr)
})



# -------------
# Download data
# -------------

# create temporary directories
DIR_TMP <- "tmp"
dir.create(file.path(DIR_TMP))
dir.create(file.path(DIR_TMP, "fcs_files"))

# download from 'imlspenticton' server
URL <- "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData"
DIR <- "Levine_13dim"

# download .fcs files
fcs_filename <- "Levine_13dim_fcs_files.zip"
download.file(file.path(URL, DIR, fcs_filename), destfile = file.path(DIR_TMP, "fcs_files", fcs_filename))
unzip(file.path(DIR_TMP, "fcs_files", fcs_filename), exdir = file.path(DIR_TMP, "fcs_files"))

files <- list.files(file.path(DIR_TMP, "fcs_files"), pattern = "\\.fcs$", full.names = TRUE)



# ---------
# Load data
# ---------

# one .fcs file per manually gated cluster
# 13 surface markers (dimensions), 24 manually gated populations

# "unassigned" cells are those where cluster labels are unavailable
files_assigned <- files[-grep("NotGated", files)]
files_unassigned <- files[grep("NotGated", files)]

# cell population names
pop_names <- 
  files_assigned %>% 
  gsub("^.*Marrow1_", "", .) %>% 
  gsub("\\.fcs$", "", .) %>% 
  gsub(" ", "_", .)

# channel and marker names
# note: '$name' and '$desc' are the wrong way around for this dataset
channel_name <- as.character(pData(parameters(read.FCS(files_assigned[1])))$desc)
marker_name <- as.character(pData(parameters(read.FCS(files_assigned[1])))$name)
# original column names
col_names <- colnames(exprs(read.FCS(files_assigned[1])))

# marker classes (cell type, cell state, or none)
marker_class <- rep("type", length(marker_name))


# load data and create vector of population IDs (for "assigned" cells)
data_assigned <- matrix(nrow = 0, ncol = length(marker_name))
pop_id_assigned <- c()

for (i in 1:length(files_assigned)) {
  # load data
  data_i <- exprs(read.FCS(files_assigned[i], transformation = FALSE, truncate_max_range = FALSE))
  data_assigned <- rbind(data_assigned, data_i)
  # population IDs
  pop_id_assigned <- c(pop_id_assigned, rep(pop_names[i], nrow(data_i)))
}

dim(data_assigned)  # 81,747 assigned cells, 13 dimensions
table(pop_id_assigned)  # 14 manually gated clusters

stopifnot(nrow(data_assigned) == length(pop_id_assigned))


# load data and create vector of population IDs (for "unassigned" cells)
data_unassigned <- matrix(nrow = 0, ncol = length(marker_name))
pop_id_unassigned <- c()

for (i in 1:length(files_unassigned)) {
  # load data
  data_i <- exprs(read.FCS(files_unassigned[i], transformation = FALSE, truncate_max_range = FALSE))
  data_unassigned <- rbind(data_unassigned, data_i)
  # population IDs
  pop_id_unassigned <- c(pop_id_unassigned, rep("unassigned", nrow(data_i)))
}

dim(data_unassigned)  # 85,297 unassigned cells

stopifnot(nrow(data_unassigned) == length(pop_id_unassigned))


# combine assigned and unassigned cells
data_all <- rbind(data_assigned, data_unassigned)
pop_id_all <- c(pop_id_assigned, pop_id_unassigned)

stopifnot(nrow(data_all) == length(pop_id_all))



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
row_data <- data.frame(
  population_id = as.factor(pop_id_all), 
  stringsAsFactors = FALSE
)

# set up column data
marker_info <- data.frame(
  channel_name = as.character(channel_name), 
  marker_name = as.character(marker_name), 
  marker_class = factor(marker_class, levels = c("none", "type", "state")), 
  stringsAsFactors = FALSE
)
col_data <- marker_info

# set up expression data
d_exprs <- data_all
# use marker names as column names (for SummarizedExperiment)
colnames(d_exprs) <- marker_name

stopifnot(nrow(d_exprs) == nrow(row_data))
stopifnot(ncol(d_exprs) == nrow(col_data))

# create SummarizedExperiment object
d_SE <- SummarizedExperiment(
  assays = list(exprs = d_exprs), 
  rowData = row_data, 
  colData = col_data
)



# ---------------------
# Create flowSet object
# ---------------------

# note: row data (e.g. population IDs) is stored as additional columns of data in the
# expression matrices; additional information from row data and column data (e.g. marker
# classes, cell population names) is stored in 'description' slot

# table of cell population information
population_info <- data.frame(
  population_id = seq_len(length(pop_names) + 1), 
  name = c(pop_names, "unassigned"), 
  stringsAsFactors = FALSE
)

# add extra columns of data and create new flowSet object
exprs_fs <- d_exprs
# use original column names (for flowSet)
colnames(exprs_fs) <- col_names
exprs_fs <- cbind(exprs_fs, population_id = as.numeric(row_data$population_id))

d_flowFrame <- flowFrame(exprs_fs)
d_flowSet <- flowSet(d_flowFrame)

# include additional information in 'description' slot
description(d_flowSet[[1]])$MARKER_INFO <- marker_info
description(d_flowSet[[1]])$POPULATION_INFO <- population_info



# ------------
# Save objects
# ------------

save(d_SE, file = "Levine_13dim_SE.rda")
save(d_flowSet, file = "Levine_13dim_flowSet.rda")



