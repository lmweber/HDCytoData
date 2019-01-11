##########################################################################################
# R script to prepare benchmark datasets Samusik_01 and Samusik_all
# 
# This is a 39-dimensional mass cytometry (CyTOF) dataset, consisting of expression levels
# of 39 surface protein markers. Cell population labels are available for 24 manually
# gated populations. The dataset contains cells from 10 replicate bone marrow samples from
# C57BL/6J mice (samples from 10 different mice): Samusik_01 contains data from sample 01
# only, and Samusik_all contains data from all samples.
# 
# This R script loads the data, adds manually gated cell population labels, and exports it
# in SummarizedExperiment and flowSet formats.
#
# Source: Samusik et al. (2016), "Automated mapping of phenotype space with single-cell
# data", Nature Methods, 13(6), 493-496.
# 
# Link to paper: https://www.ncbi.nlm.nih.gov/pubmed/27183440
# Link to data (.zip file): "https://web.stanford.edu/~samusik/Panorama BM 1-10.zip"
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

# one FCS file per sample (samples 01 to 10)
# 39 surface markers (dimensions), 24 manually gated populations

# create temporary directories
DIR_TMP <- "tmp"
dir.create(file.path(DIR_TMP))
dir.create(file.path(DIR_TMP, "fcs_files"))
dir.create(file.path(DIR_TMP, "population_IDs"))

# download from 'imlspenticton' server
URL <- "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData"
DIR <- "Samusik"

# download .fcs files
fcs_filename <- "Samusik_fcs_files.zip"
download.file(file.path(URL, DIR, fcs_filename), destfile = file.path(DIR_TMP, "fcs_files", fcs_filename))
unzip(file.path(DIR_TMP, "fcs_files", fcs_filename), exdir = file.path(DIR_TMP, "fcs_files"))

files_load_fcs <- list.files(file.path(DIR_TMP, "fcs_files"), pattern = "\\.fcs$", full.names = TRUE)

# download population IDs
pop_filename <- "Samusik_population_IDs.zip"
download.file(file.path(URL, DIR, pop_filename), destfile = file.path(DIR_TMP, "population_IDs", pop_filename))
unzip(file.path(DIR_TMP, "population_IDs", pop_filename), exdir = file.path(DIR_TMP, "population_IDs"))

file_gating <- file.path(DIR_TMP, "population_IDs", "population_assignments.txt")



# ----------------------
# Load population labels
# ----------------------

data_gating <- read.table(file_gating, header = FALSE, stringsAsFactors = FALSE, sep = "\t")

# remove line with error (line is partially cut off) and reset row names
data_gating[164401, ]
data_gating <- data_gating[-164401, ]
rownames(data_gating) <- NULL

dim(data_gating)  # 514,386 cells

# extract sample numbers
sample <- as.factor(sapply(strsplit(data_gating[, 1], split = " "), function(s) s[1]))
str(sample)
length(sample)

# extract event (cell) numbers
event <- sapply(strsplit(data_gating[, 1], split = " "), function(s) s[3])
event <- as.numeric(event)
str(event)

# add 1 to event numbers (event numbers are provided as index-0, but R-based row numbers
# in the .fcs files are index-1)
event <- event + 1

all(event == floor(event))  # check: all integers
sum(event == 1)  # check: multiple events no. 1 (but not necessarily one for every sample due to unassigned cells)

# split event numbers into one data frame per sample
event <- split(event, sample)
str(event)

# extract population labels
population <- data_gating[, 2]
str(population)

# split population labels into one data frame per sample
population <- split(population, sample)
str(population)



# ---------------------------------
# Load data and create row metadata
# ---------------------------------

# note: 'unassigned' refers to cells not assigned to any population by manual gating

# channel and marker names
channel_name <- as.character(pData(parameters(read.FCS(files_load_fcs[1])))$name)
marker_name <- as.character(pData(parameters(read.FCS(files_load_fcs[1])))$desc)
# original column names
col_names <- colnames(exprs(read.FCS(files_load_fcs[1])))

# marker classes (cell type, cell state, or none)
marker_class <- rep("none", length(marker_name))
marker_class[c(9:47)] <- "type"

data <- list()
row_data <- list()

for (i in 1:length(files_load_fcs)) {
  data_i <- exprs(read.FCS(files_load_fcs[i], transformation = FALSE, truncate_max_range = FALSE))
  
  names_i <- basename(files_load_fcs[i])
  
  # sample IDs
  sample_i <- data.frame(sample_id = rep(gsub("^.*_([0-9]+)_.*$", "\\1", files_load_fcs[[i]]), nrow(data_i)))
  # population IDs ('unassigned' = not assigned to any population by manual gating)
  pop_i <- data.frame(population_id = rep("unassigned", nrow(data_i)), stringsAsFactors = FALSE)
  pop_i[event[[i]], "population_id"] <- population[[i]]
  # convert to factor
  pop_i[, "population_id"] <- as.factor(pop_i[, "population_id"])
  
  data[[i]] <- data_i
  names(data)[i] <- names_i
  
  row_data[[i]] <- cbind(sample_i, pop_i)
  names(row_data)[i] <- names_i
}

# number of cells per sample (check against Supp. Table. S2 Excel spreadsheet)
n_cells <- sapply(data, nrow)
n_cells
# number of assigned cells per sample
n_assigned <- sapply(event, length)
n_assigned
# proportion assigned
prop_assigned <- n_assigned / n_cells
prop_assigned
# total number of cells
sum(n_cells)

# check column names
for (i in seq_along(data)) {
  stopifnot(all(colnames(data[[i]]) == channel_name))
}



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ----------------------------------
# Create SummarizedExperiment object
# ----------------------------------

# set up row data
row_data <- do.call("rbind", row_data)
rownames(row_data) <- NULL

stopifnot(nrow(row_data) == sum(n_cells))

# set up column data
col_data <- data.frame(
  channel_name = as.character(channel_name), 
  marker_name = as.character(marker_name), 
  marker_class = factor(marker_class, levels = c("none", "type", "state")), 
  stringsAsFactors = FALSE
)

# set up expression data
d_exprs <- do.call("rbind", data)
# use marker names as column names (for SummarizedExperiment)
colnames(d_exprs) <- marker_name

stopifnot(nrow(d_exprs) == nrow(row_data))
stopifnot(ncol(d_exprs) == nrow(col_data))

# create SummarizedExperiment object for 'Samusik_all' (full dataset)
d_SE_all <- SummarizedExperiment(
  assays = list(exprs = d_exprs), 
  rowData = row_data, 
  colData = col_data
)

# create SummarizedExperiment object for 'Samusik_01' (sample 01)
d_SE_01 <- d_SE_all[rowData(d_SE_all)$sample_id == "01", ]



# --------------
# Create flowSet
# --------------

# note: population IDs are stored as an additional column of data in the expression
# matrices; additional marker information (marker names and marker classes) cannot be
# included here, since marker information is stored in column names only

# create table of cell population names
df_population_names <- data.frame(
  population_id = 1:nlevels(row_data$population_id), 
  population_name = levels(row_data$population_id), 
  stringsAsFactors = FALSE
)

# create list of extra columns of data for each sample
row_data_fs <- row_data[, "population_id", drop = FALSE]
row_data_fs$population_id <- as.numeric(row_data_fs$population_id)

row_data_fs_list <- lapply(split(row_data_fs, row_data$sample_id), as.matrix)

sample_id_names <- names(table(row_data$sample_id))
stopifnot(all(names(row_data_fs_list) == sample_id_names))

stopifnot(all(sapply(data, nrow) == sapply(row_data_fs_list, nrow)))

# add extra columns of data and create new flowSet object
exprs_fs_list <- data
exprs_fs_list <- lapply(exprs_fs_list, function(e) {
  # use original column names (for flowSet)
  colnames(e) <- col_names
  e
})
d_flowFrames_list <- mapply(function(e, extra_cols) {
  # combine and create flowFrame
  stopifnot(nrow(e) == nrow(extra_cols))
  ff <- flowFrame(cbind(e, extra_cols))
  # include both channel and marker names in 'pData(parameters(.))'
  stopifnot(length(c(marker_name, colnames(extra_cols))) == nrow(pData(parameters(ff))))
  pData(parameters(ff))$desc <- c(marker_name, colnames(extra_cols))
  ff
}, exprs_fs_list, row_data_fs_list)


# create flowSet object for 'Samusik_all' (full dataset)
d_flowSet_all <- flowSet(d_flowFrames_list)

# include filenames, sample IDs, and table of population names in 'description' slot
for (i in seq_along(d_flowSet_all)) {
  description(d_flowSet_all[[i]])$FILENAME <- identifier(d_flowSet_all[[i]])
  description(d_flowSet_all[[i]])$SAMPLE_ID <- gsub("^.*_([0-9]+)_.*$", "\\1", identifier(d_flowSet_all[[i]]))
  stopifnot(description(d_flowSet_all[[i]])$SAMPLE_ID == sample_id_names[i])
  description(d_flowSet_all[[i]])$POPULATION_NAMES <- df_population_names
}

# create flowSet object for 'Samusik_01' (sample 01)
d_flowSet_01 <- flowSet(d_flowFrames_list[1])

# include filenames, sample IDs, and table of population names in 'description' slot
for (i in seq_along(d_flowSet_01)) {
  description(d_flowSet_01[[i]])$FILENAME <- identifier(d_flowSet_01[[i]])
  description(d_flowSet_01[[i]])$SAMPLE_ID <- gsub("^.*_([0-9]+)_.*$", "\\1", identifier(d_flowSet_01[[i]]))
  stopifnot(description(d_flowSet_all[[i]])$SAMPLE_ID == sample_id_names[i])
  description(d_flowSet_01[[i]])$POPULATION_NAMES <- df_population_names
}



# ------------
# Save objects
# ------------

save(d_SE_all, file = "Samusik_all_SE.rda")
save(d_flowSet_all, file = "Samusik_all_flowSet.rda")

save(d_SE_01, file = "Samusik_01_SE.rda")
save(d_flowSet_01, file = "Samusik_01_flowSet.rda")



