##########################################################################################
# R script to prepare benchmark dataset Mosmann_rare
# 
# This is a 14-dimensional flow cytometry dataset containing a rare population of
# activated (cytokine-producing) memory CD4 T cells, from human peripheral blood
# mononuclear cells (PBMCs) exposed to influenza antigens, from a single healthy donor.
# 
# This R script loads the data, adds manually gated cell population labels for the rare
# population, and exports it in SummarizedExperiment and flowSet formats.
# 
# Source: Figure 4 in the following paper:
# Mosmann et al. (2014), "SWIFT — Scalable Clustering for Automated Identification of Rare
# Cell Populations in Large, High-Dimensional Flow Cytometry Datasets, Part 2: Biological
# Evaluation", Cytometry Part A, 85A, 422-433.
#
# Link to paper: https://www.ncbi.nlm.nih.gov/pubmed/24532172
# Link to data: http://flowrepository.org/id/FR-FCM-ZZ8J
# (filename: "JMW034-J16OFVQX_G2 0o1 3_D07.fcs"; see Supplementary Information file 3 for
# full list of filenames)
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
DIR <- "Mosmann_rare"

# download .fcs files
fcs_filename <- "Mosmann_rare_fcs_files.zip"
download.file(file.path(URL, DIR, fcs_filename), destfile = file.path(DIR_TMP, "fcs_files", fcs_filename))
unzip(file.path(DIR_TMP, "fcs_files", fcs_filename), exdir = file.path(DIR_TMP, "fcs_files"))

files <- list.files(file.path(DIR_TMP, "fcs_files"), pattern = "\\.fcs$", full.names = TRUE)



# -------------------
# Load data: raw data
# -------------------

file_raw <- file.path(DIR_TMP, "fcs_files", "JMW034-J16OFVQX_G2 0o1 3_D07.fcs")
data_raw <- read.FCS(file_raw, transformation = FALSE, truncate_max_range = FALSE)

dim(data_raw)

# channel and marker names
channel_name <- as.character(pData(parameters(data_raw))$name)
marker_name <- as.character(pData(parameters(data_raw))$desc)
# clean marker names
ix_NAs <- is.na(marker_name)
marker_name[ix_NAs] <- channel_name[ix_NAs]
marker_name <- gsub(" .*$", "", marker_name)
marker_name[10] <- gsub("/", "_", marker_name[10])
marker_name[c(11, 19)] <- c("GZB-SA", "CCL4")  # from Supporting Information file 1
# original column names
col_names <- colnames(exprs(data_raw))

# marker classes (cell type, cell state, or none)
marker_class <- rep("none", length(marker_name))
marker_class[c(7, 12, 13, 14, 15, 18, 20)] <- "type"  # from paper and Supporting Information file 1
marker_class[(c(8, 9, 11, 16, 17, 19, 21))] <- "state"  # from paper and Supporting Information file 1

# spillover matrix (not required since compensation is done automatically in Cytobank)
#description(data_raw)$SPILL



# -----------------------------
# Load data: single, live cells
# -----------------------------

# pre-gating to exclude debris and doublets was done in Cytobank, following the gating
# scheme shown in the first two panels of Figure 4A in Mosmann et al. (2014)

# note: compensation is done automatically in Cytobank using the spillover matrix from the
# raw data file

file_single_live <- file.path(DIR_TMP, "fcs_files", "JMW034-J16OFVQX_G2 0o1 3_D07_Live_cells.fcs")
data_single_live <- exprs(read.FCS(file_single_live, transformation = FALSE, truncate_max_range = FALSE))

dim(data_single_live)  # 396,460 cells



# ---------------------------------------------------
# Load data: rare population of activated CD4 T cells
# ---------------------------------------------------

# gating was done in Cytobank, following the gating scheme in Figure 4A in Mosmann et al. (2014)

file_activ <- file.path(DIR_TMP, "fcs_files", "JMW034-J16OFVQX_G2 0o1 3_D07_IFNg_vs_TNFa.fcs")
data_activ <- exprs(read.FCS(file_activ, transformation = FALSE, truncate_max_range = FALSE))

dim(data_activ)  # 109 cells

# rare cell proportion
nrow(data_activ) / nrow(data_single_live)  # <0.03% of single, live cells



# -----------------
# Population labels
# -----------------

# since the .fcs files do not include event numbers or identifiers, we identify cells from
# the activated CD4 T cell population by duplicating rows

data_dup <- rbind(data_single_live, data_activ)
dim(data_dup)

ix_dup <- duplicated(data_dup, fromLast = TRUE)

sum(ix_dup)  # 109 cells
length(ix_dup)

# create vector of labels
labels <- ix_dup[1:nrow(data_single_live)]

sum(labels)  # 109 rare cells
length(labels)  # 396,460 single, live cells
table(labels)

stopifnot(nrow(data_single_live) == length(labels))



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
population_id <- factor(as.numeric(labels), labels = c("other", "activated"))

row_data <- data.frame(
  population_id = population_id, 
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
d_exprs <- data_single_live
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
  population_id = c(0, 1), 
  name = c("other", "activated"), 
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

# include additional information in 'description' slot
description(d_flowSet[[1]])$MARKER_INFO <- marker_info
description(d_flowSet[[1]])$POPULATION_INFO <- population_info



# ------------
# Save objects
# ------------

save(d_SE, file = "Mosmann_rare_SE.rda")
save(d_flowSet, file = "Mosmann_rare_flowSet.rda")



