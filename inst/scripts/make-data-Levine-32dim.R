##########################################################################################
# R script to prepare benchmark data set Levine_32dim
# 
# This is a 32-dimensional mass cytometry (CyTOF) data set, consisting of expression
# levels of 32 surface protein markers. Cell population labels are available for 14
# manually gated populations. Cells are healthy human bone marrow mononuclear cells
# (BMMCs), from 2 patients.
#
# This R script loads the data, adds manually gated cell population labels, and exports it
# in SummarizedExperiment and flowSet formats.
#
# Source: "benchmark data set 2" in the following paper:
# Levine et al. (2015), "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like
# Cells that Correlate with Prognosis", Cell, 162, 184-197.
#
# Link to paper: https://www.ncbi.nlm.nih.gov/pubmed/26095251
# Link to raw data: https://www.cytobank.org/cytobank/experiments/46102 (download the .zip
# file shown under "Exported Files")
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
DIR <- "Levine_32dim"

# download .fcs files
fcs_filename <- "Levine_32dim_fcs_files.zip"
download.file(file.path(URL, DIR, fcs_filename), destfile = file.path(DIR_TMP, "fcs_files", fcs_filename))
unzip(file.path(DIR_TMP, "fcs_files", fcs_filename), exdir = file.path(DIR_TMP, "fcs_files"))

files <- list.files(file.path(DIR_TMP, "fcs_files"), pattern = "\\.fcs$", full.names = TRUE)



# ---------
# Load data
# ---------

# one .fcs file per manually gated cluster, per patient (H1 and H2)
# 32 surface markers (dimensions), 14 manually gated populations, 2 patients (H1 and H2)

# "unassigned" cells are those where cluster labels are unavailable
files_assigned <- files[-grep("NotDebrisSinglets", files)]
files_unassigned <- files[grep("NotDebrisSinglets", files)]

# cell population names
pop_names <- 
  grep("H1\\.fcs$", files_assigned) %>% 
  files_assigned[.] %>% 
  gsub("(_H1\\.fcs$)", "", .) %>% 
  gsub("^.*_", "", .) %>% 
  gsub(" ", "_", .)

# column names (protein markers and others)
col_names <- 
  read.FCS(files_assigned[1], transformation = FALSE, truncate_max_range = FALSE) %>% 
  exprs %>% 
  colnames %>% 
  unname

# clean column names
channel_name <- col_names
marker_name <- gsub("\\(.*", "", col_names)

# marker classes (cell type, cell state, or none)
marker_class <- rep("none", length(col_names))
marker_class[5:36] <- "type"

# vector of labels for patients H1 and H2
patient_names_assigned <- rep(NA, length(files_assigned))
patient_names_assigned[grep("_H1\\.fcs$", files_assigned)] <- "H1"
patient_names_assigned[grep("_H2\\.fcs$", files_assigned)] <- "H2"

patient_names_unassigned <- c("H1", "H2")


# load data and create vectors of population IDs and sample IDs (for "assigned" cells)
data_assigned <- matrix(nrow = 0, ncol = length(col_names))
pop_id_assigned <- c()
patient_id_assigned <- c()

for (i in 1:length(files_assigned)) {
  # load data
  data_i <- exprs(read.FCS(files_assigned[i], transformation = FALSE, truncate_max_range = FALSE))
  data_assigned <- rbind(data_assigned, data_i)
  # population IDs
  pop_id_assigned <- c(pop_id_assigned, rep(pop_names[((i - 1) %% length(pop_names)) + 1], nrow(data_i)))
  # patient IDs
  patient_id_assigned <- c(patient_id_assigned, rep(patient_names_assigned[i], nrow(data_i)))
}

dim(data_assigned)  # 104,184 assigned cells, 32 dimensions
table(pop_id_assigned)  # 14 manually gated clusters
table(patient_id_assigned) # 2 patient (72,463 and 31,721 assigned cells each)

stopifnot(nrow(data_assigned) == length(pop_id_assigned))
stopifnot(nrow(data_assigned) == length(patient_id_assigned))


# load data and create vectors of population IDs and sample IDs (for "unassigned" cells)
data_unassigned <- matrix(nrow = 0, ncol = length(col_names))
pop_id_unassigned <- c()
patient_id_unassigned <- c()

for (i in 1:length(files_unassigned)) {
  # load data
  data_i <- exprs(read.FCS(files_unassigned[i], transformation = FALSE, truncate_max_range = FALSE))
  data_unassigned <- rbind(data_unassigned, data_i)
  # population IDs
  pop_id_unassigned <- c(pop_id_unassigned, rep("unassigned", nrow(data_i)))
  # patient IDs
  patient_id_unassigned <- c(patient_id_unassigned, rep(patient_names_unassigned[i], nrow(data_i)))
}

dim(data_unassigned)  # 161,443 unassigned cells
table(patient_id_unassigned) # 2 patients (118,888 and 42,555 unassigned cells each)

stopifnot(nrow(data_unassigned) == length(pop_id_unassigned))
stopifnot(nrow(data_unassigned) == length(patient_id_unassigned))


# combine assigned and unassigned cells
data_all <- rbind(data_assigned, data_unassigned)
pop_id_all <- c(pop_id_assigned, pop_id_unassigned)
patient_id_all <- c(patient_id_assigned, patient_id_unassigned)

stopifnot(nrow(data_all) == length(pop_id_all))
stopifnot(nrow(data_all) == length(patient_id_all))



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
row_data <- data.frame(
  patient_id = as.factor(patient_id_all), 
  population_id = as.factor(pop_id_all), 
  stringsAsFactors = FALSE
)

# set up column data
col_data <- data.frame(
  channel_name = as.character(channel_name), 
  marker_name = as.character(marker_name), 
  marker_class = as.factor(marker_class), 
  stringsAsFactors = FALSE
)

# set up expression data
d_exprs <- data_all
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
  population_id = 1:(length(pop_names) + 1), 
  population_name = c(pop_names, "unassigned"), 
  stringsAsFactors = FALSE
)

# create list of extra columns of data for each sample
cols_keep <- 2
row_data_fs <- do.call("cbind", lapply(row_data[, cols_keep, drop = FALSE], as.numeric))
stopifnot(nrow(row_data_fs) == nrow(row_data))

row_data_fs_list <- lapply(split(row_data_fs, row_data$patient_id), matrix, ncol = ncol(row_data_fs))

patient_id_names <- names(table(row_data$patient_id))
stopifnot(all(names(row_data_fs_list) == patient_id_names))

# replace column names
row_data_fs_list <- lapply(row_data_fs_list, function(d) {
  colnames(d) <- colnames(row_data_fs)
  d
})

# create new flowSet object and add extra columns of data
exprs_fs_list <- lapply(split(d_exprs, row_data$patient_id), matrix, ncol = ncol(d_exprs))
exprs_fs_list <- lapply(exprs_fs_list, function(e) {
  colnames(e) <- channel_name
  e
})

d_flowFrames_list <- mapply(function(e, extra_cols) {
  # combine and create flowFrame
  stopifnot(nrow(e) == nrow(extra_cols))
  flowFrame(cbind(e, extra_cols))
}, exprs_fs_list, row_data_fs_list)

d_flowSet <- flowSet(d_flowFrames_list)

# include patient IDs and table of population names in 'description' slot
for (i in seq_along(d_flowSet)) {
  description(d_flowSet[[i]])$PATIENT_ID <- identifier(d_flowSet[[i]])
  description(d_flowSet[[i]])$POPULATION_NAMES <- df_population_names
}



# ------------
# Save objects
# ------------

save(d_SE, file = "Levine_32dim_SE.rda")
save(d_flowSet, file = "Levine_32dim_flowSet.rda")



