##########################################################################################
# R script to prepare benchmark dataset Nilsson_rare
# 
# This is a 13-dimensional flow cytometry dataset containing a rare population of
# hematopoietic stem cells (HSCs), from human bone marrow cells from a single healthy
# donor.
# 
# This R script loads the data, adds manually gated cell population labels for the rare
# population, and exports it in SummarizedExperiment and flowSet formats.
# 
# Source: Figure 2 in the following paper:
# Nilsson et al. (2013), "Frequency Determination of Rare Populations by Flow Cytometry: A
# Hematopoietic Stem Cell Perspective", Cytometry Part A, 83A, 721-727.
# 
# Link to paper: http://www.ncbi.nlm.nih.gov/pubmed/23839904
# Link to data: http://flowrepository.org/id/FR-FCM-ZZ6L
# 
# Lukas Weber, Jan 2019
##########################################################################################


# original version of this script can be found at:
# https://github.com/lmweber/cytometry-clustering-comparison


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

# download from 'imlspenticton' server
URL <- "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData"
DIR <- "Nilsson_rare"

# download .fcs files
fcs_filename <- "Nilsson_rare_fcs_files.zip"
download.file(file.path(URL, DIR, fcs_filename), destfile = file.path(DIR_TMP, "fcs_files", fcs_filename))
unzip(file.path(DIR_TMP, "fcs_files", fcs_filename), exdir = file.path(DIR_TMP, "fcs_files"))

files <- list.files(file.path(DIR_TMP, "fcs_files"), pattern = "\\.fcs$", full.names = TRUE)



# -------------------
# Load data: raw data
# -------------------

file_raw <- file.path(DIR_TMP, "fcs_files", "Pronk CytometryA Figure2.fcs")
data_raw <- read.FCS(file_raw, transformation = FALSE, truncate_max_range = FALSE)

dim(data_raw)

# channel and marker names
channel_name <- as.character(pData(parameters(data_raw))$name)
marker_name <- as.character(pData(parameters(data_raw))$desc)
# clean marker names
ix_NAs <- is.na(marker_name)
marker_name[ix_NAs] <- channel_name[ix_NAs]
marker_name <- gsub(" .*$", "", marker_name)
# original column names
col_names <- colnames(exprs(data_raw))

# marker classes (cell type, cell state, or none)
marker_class <- rep("none", length(marker_name))
ix_markers <- grep("^CD", marker_name)
marker_class[ix_markers] <- "type"

# spillover matrix (not required since compensation is done automatically in Cytobank)
#description(data_raw)$SPILL



# -----------------------------------------
# Load data: single cells, non-debris, live
# -----------------------------------------

# pre-gating to exclude doublets, debris, and dead cells was done in Cytobank, following
# the gating scheme shown in the first 3 panels of Figure 2 in Nilsson et al. (2013)

# note: compensation is done automatically in Cytobank using the spillover matrix from the
# raw data file

file_single_nondeb_live <- file.path(DIR_TMP, "fcs_files", "Pronk CytometryA Figure2_Live cells.fcs")
data_single_nondeb_live <- exprs(read.FCS(file_single_nondeb_live, transformation = FALSE, truncate_max_range = FALSE))

dim(data_single_nondeb_live)  # 44,140 cells



# -------------------------------------------------------------
# Load data: rare population of hematopoietic stem cells (HSCs)
# -------------------------------------------------------------

# gating was done in Cytobank, following the gating scheme in Figure 2 in Nilsson et al. (2013)

file_HSC <- file.path(DIR_TMP, "fcs_files", "Pronk CytometryA Figure2_HSCs.fcs")
data_HSC <- exprs(read.FCS(file_HSC, transformation = FALSE, truncate_max_range = FALSE))

dim(data_HSC)  # 358 cells

# rare cell proportion
nrow(data_HSC) / nrow(data_single_nondeb_live)  # 0.8% of single, non-debris, live cells



# -----------------
# Population labels
# -----------------

# since the .fcs files do not include event numbers or identifiers, we identify cells from
# the HSC population by duplicating rows

data_dup <- rbind(data_single_nondeb_live, data_HSC)
dim(data_dup)

ix_dup <- duplicated(data_dup, fromLast = TRUE)

sum(ix_dup)  # 358 cells
length(ix_dup)

# create vector of labels
labels <- ix_dup[1:nrow(data_single_nondeb_live)]

sum(labels)  # 358 HSCs
length(labels)  # 44,140 single, non-debris, live cells
table(labels)

stopifnot(nrow(data_single_nondeb_live) == length(labels))



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
population_id <- factor(as.numeric(labels), labels = c("other", "HSCs"))

row_data <- data.frame(
  population_id = population_id, 
  stringsAsFactors = FALSE
)

# set up column data
col_data <- data.frame(
  channel_name = as.character(channel_name), 
  marker_name = as.character(marker_name), 
  marker_class = factor(marker_class, levels = c("none", "type", "state")), 
  stringsAsFactors = FALSE
)

# set up expression data
d_exprs <- data_single_nondeb_live
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



# --------------
# Create flowSet
# --------------

# note: population IDs are stored as an additional column of data in the expression
# matrices; additional marker information (marker names and marker classes) cannot be
# included here, since marker information is stored in column names only

# create table of cell population names
df_population_names <- data.frame(
  population_id = c(0, 1), 
  name = c("other", "HSCs"), 
  stringsAsFactors = FALSE
)

# add extra columns of data and create new flowSet object
exprs_fs <- d_exprs
# use original column names (for flowSet)
colnames(exprs_fs) <- col_names
exprs_fs <- cbind(exprs_fs, population_id = as.numeric(labels))

d_flowFrame <- flowFrame(exprs_fs)

# include both channel and marker names in 'pData(parameters(.))'
stopifnot(length(marker_name) + 1 == nrow(pData(parameters(d_flowFrame))))
pData(parameters(d_flowFrame))$desc <- c(marker_name, "population_id")

d_flowSet <- flowSet(d_flowFrame)

# include table of population names in 'description' slot
description(d_flowSet[[1]])$POPULATION_NAMES <- df_population_names



# ------------
# Save objects
# ------------

save(d_SE, file = "Nilsson_rare_SE.rda")
save(d_flowSet, file = "Nilsson_rare_flowSet.rda")



