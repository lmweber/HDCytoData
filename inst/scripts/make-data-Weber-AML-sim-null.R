##########################################################################################
# R script to prepare benchmark dataset 'Weber-AML-sim-null'
# 
# See Weber et al. (2019), Supplementary Note 1 (paper introducing 'diffcyt' framework)
# for more details
# 
# The 'Weber-AML-sim' dataset is constructed by computationally 'spiking in' small
# percentages of AML (acute myeloid leukemia) blast cells into samples of healthy BMMCs
# (bone marrow mononuclear cells). This simulates the phenotype of minimal residual
# disease (MRD) in AML patients. Raw data is sourced from Levine et al. (2015) (PhenoGraph
# paper). The data generation strategy is modified from Arvaniti et al. (2017) (CellCnn
# paper), who generated a similar benchmark dataset for their evaluations.
# 
# Raw data downloaded from Cytobank:
# - all cells (also contains gating scheme for CD34+CD45mid cells, i.e. blasts):
# https://community.cytobank.org/cytobank/experiments/46098/illustrations/121588
# - blasts (repository cloned from the one for 'all cells' above, using the gating scheme
# for CD34+CD45mid cells; this allows .fcs files for the subset to be exported):
# https://community.cytobank.org/cytobank/experiments/63534/illustrations/125318
# 
# Notes:
# - Gating plots for blasts are also shown in Levine et al. (2015), Supplemental Data S3B.
# - Individuals SJ1, SJ2, and SJ3 each contain two replicates in separate .fcs files. The
# original Cytobank repository combines the two replicates for each individual (see
# 'Individuals' dimension setup); so the combined cells from both .fcs files should be
# used for downstream analysis. However, for SJ1, the total percentage of blasts does not
# match to the published numbers (Levine et al. 2015, Supplemental Data S3B); so we have
# not used these samples.
# - Arvaniti et al. (2017) (CellCnn paper) classified patients SJ10, SJ12, SJ13 as CN
# (cytogenetically normal), and SJ1, SJ2, SJ3, SJ4, SJ5 as CBF (core-binding factor
# translocation); we re-use these classifications here.
# - Sample names and filenames in the raw data are shuffled (e.g. file H3 actually refers
# to sample H1). The matching scheme can be seen in the 'Individuals' setup in Cytobank,
# or in the downloaded .tsv files 'experiment_46098_annotations.tsv' and
# 'experiment_63534_annotations.tsv'.
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

# note: 'experiment_46098' contains all cells; 'experiment_63534' contains subset of CD34+CD45mid blast cells

# create temporary directories
DIR_TMP <- "tmp"
dir.create(file.path(DIR_TMP))
dir.create(file.path(DIR_TMP, "experiment_46098"))
dir.create(file.path(DIR_TMP, "experiment_63534"))

# download from 'imlspenticton' server
URL <- "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData"
DIR <- "Levine_AML"

# download files
fn_experiment_46098 <- "Levine_AML_experiment_46098_files.zip"
fn_experiment_63534 <- "Levine_AML_experiment_63534_files.zip"

download.file(file.path(URL, DIR, fn_experiment_46098), destfile = file.path(DIR_TMP, "experiment_46098", fn_experiment_46098))
download.file(file.path(URL, DIR, fn_experiment_63534), destfile = file.path(DIR_TMP, "experiment_63534", fn_experiment_63534))

unzip(file.path(DIR_TMP, "experiment_46098", fn_experiment_46098), exdir = file.path(DIR_TMP, "experiment_46098"))
unzip(file.path(DIR_TMP, "experiment_63534", fn_experiment_63534), exdir = file.path(DIR_TMP, "experiment_63534"))



# ---------
# Filenames
# ---------

DIR_RAW_DATA_ALL <- file.path(DIR_TMP, "experiment_46098")
DIR_RAW_DATA_BLASTS <- file.path(DIR_TMP, "experiment_63534")

files_all <- list.files(DIR_RAW_DATA_ALL, pattern = "\\.fcs$", full.names = TRUE)
files_healthy <- files_all[grep("H[0-9]+", files_all)]

files_blasts <- list.files(DIR_RAW_DATA_BLASTS, pattern = "\\.fcs$", full.names = TRUE)

# metadata spreadsheets to match shuffled sample names
file_match_samples_all <- file.path(DIR_RAW_DATA_ALL, "experiment_46098_annotations.tsv")
file_match_samples_blasts <- file.path(DIR_RAW_DATA_BLASTS, "experiment_63534_annotations.tsv")



# -----------------------------------------------
# Load data for healthy samples H1-H5 (all cells)
# -----------------------------------------------

data_healthy <- lapply(files_healthy, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

# show shuffled sample names
tbl_match_healthy <- read.delim(file_match_samples_all)
tbl_match_healthy_sub <- tbl_match_healthy[grep("H[0-9]+", tbl_match_healthy[, "FCS.Filename"]), ]
tbl_match_healthy_sub[, c("FCS.Filename", "Individuals")]

# match correct sample names; store as names of list items
names(data_healthy) <- tbl_match_healthy_sub[, "Individuals"]

length(data_healthy)
sapply(data_healthy, dim)



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# ---------------------
# Randomized replicates
# ---------------------

# use different random seed for each replicate

n_replicates <- 3

data_replicates <- vector("list", n_replicates)
names(data_replicates) <- paste0("seed", 1:n_replicates)


for (r in 1:n_replicates) {
  
  data_replicates[[r]] <- vector("list", 2)
  names(data_replicates[[r]]) <- c("null1", "null2")
  
  
  
  # ---------------------
  # Split healthy samples
  # ---------------------
  
  # Null simulation: split each healthy sample (H1-H5) into 2 equal parts.
  
  data_healthy_null_1 <- data_healthy_null_2 <- vector("list", length(data_healthy))
  names(data_healthy_null_1) <- names(data_healthy_null_2) <- names(data_healthy)
  
  # modified random seed for each replicate
  seed <- 10000 + 100 * r
  
  for (i in 1:length(data_healthy)) {
    data_i <- data_healthy[[i]]
    
    # modified random seed for each replicate
    set.seed(seed + i)
    
    # null simulation: 2 equal parts
    n <- round(nrow(data_i) / 2)
    
    ix_null_1 <- sample(1:nrow(data_i), n)
    ix_null_2 <- setdiff(1:nrow(data_i), ix_null_1)
    
    data_healthy_null_1[[i]] <- data_i[ix_null_1, ]
    data_healthy_null_2[[i]] <- data_i[ix_null_2, ]
  }
  
  sapply(data_healthy_null_1, dim)
  sapply(data_healthy_null_2, dim)
  
  # store data
  data_replicates[[r]][[1]] <- data_healthy_null_1
  data_replicates[[r]][[2]] <- data_healthy_null_2
  
}



# ---------------
# Create metadata
# ---------------

# sample information

conditions <- c("null1", "null2")

sample_id <- paste0(rep(conditions, each = length(data_healthy)), "_", names(data_healthy))
# convert to factor with levels in expected order
sample_id <- factor(sample_id, levels = sample_id)

group_id <- factor(gsub("_.*$", "", sample_id), levels = conditions)
group_id

patient_id <- factor(gsub("^.*_", "", sample_id))
patient_id

experiment_info <- data.frame(group_id, patient_id, sample_id, stringsAsFactors = FALSE)
experiment_info


# marker information

# indices of all marker columns, lineage markers, and functional markers
# (16 surface markers / 15 functional markers; see Levine et al. 2015, Supplemental 
# Information, p. 4)
cols_markers <- 11:41
cols_lineage <- c(35, 29, 14, 30, 12, 26, 17, 33, 41, 32, 22, 40, 27, 37, 23, 39)
cols_func <- setdiff(cols_markers, cols_lineage)

stopifnot(all(sapply(seq_along(data_healthy), function(i) all(colnames(data_healthy[[i]]) == colnames(data_healthy[[1]])))))

# channel and marker names
channel_name <- colnames(data_healthy[[1]])
names(channel_name) <- NULL
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

# create a single object containing separate assays for each replicate


# check numbers of cells are identical across replicates
for (r in 1:n_replicates) {
  stopifnot(identical(sapply(data_replicates[[r]][[1]], dim), sapply(data_replicates[[1]][[1]], dim)))
  stopifnot(identical(sapply(data_replicates[[r]][[2]], dim), sapply(data_replicates[[1]][[2]], dim)))
}

patients_nm <- names(data_healthy)


# set up row data
n_cells_null1 <- sapply(data_replicates[[1]][[1]], dim)[1, ]
n_cells_null2 <- sapply(data_replicates[[1]][[2]], dim)[1, ]

n_cells_z <- c(n_cells_null1, n_cells_null2)
stopifnot(length(n_cells_z) == nrow(experiment_info))

stopifnot(all(names(n_cells_z) == gsub("^.*_", "", experiment_info$sample_id)))
names(n_cells_z) <- experiment_info$sample_id

row_data <- as.data.frame(lapply(experiment_info, function(col) {
  as.factor(rep(col, n_cells_z))
}))

stopifnot(nrow(row_data) == sum(n_cells_z))

# add spike-in column
row_data$spikein <- FALSE


# set up column data
col_data <- marker_info


# set up expression data
# note: one assay per replicate
d_exprs <- vector("list", length(data_replicates))
names(d_exprs) <- names(data_replicates)

for (r in 1:length(d_exprs)) {
  data_z <- matrix(, nrow = 0, ncol = ncol(data_replicates[[1]][[1]][[1]]))
  colnames(data_z) <- colnames(data_replicates[[1]][[1]][[1]])
  
  for (i in patients_nm) {
    stopifnot(colnames(data_z) == colnames(data_replicates[[r]][["null1"]][[i]]))
    data_z <- rbind(data_z, data_replicates[[r]][["null1"]][[i]])
  }
  for (i in patients_nm) {
    stopifnot(colnames(data_z) == colnames(data_replicates[[r]][["null2"]][[i]]))
    data_z <- rbind(data_z, data_replicates[[r]][["null2"]][[i]])
  }
  
  stopifnot(nrow(data_z) == nrow(row_data))
  stopifnot(nrow(data_z) == sum(n_cells_z))
  stopifnot(all(colnames(data_z) == marker_info$channel_name))
  stopifnot(ncol(data_z) == nrow(col_data))
  
  d_exprs[[r]] <- data_z
}


# create SummarizedExperiment objects
d_SE <- SummarizedExperiment(
  assays = d_exprs, 
  rowData = row_data, 
  colData = col_data, 
  metadata = list(experiment_info = experiment_info, n_cells = n_cells_z)
)



# ----------------------
# Create flowSet objects
# ----------------------

# note: row data is stored as additional columns of data in the expression matrices;
# additional information from row data and column data (e.g. marker classes) is stored in
# 'description' slot

# create a single flowSet object containing separate flowFrames for each replicate


# extract data
row_data <- rowData(d_SE)
col_data <- colData(d_SE)
d_exprs <- assays(d_SE)
meta_data <- metadata(d_SE)

# create tables to identify row data values when converted to numeric
stopifnot(all(colnames(row_data) == c("group_id", "patient_id", "sample_id", "spikein")))
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
spikein_info <- data.frame(
  spikein = c(0, 1), 
  spikein_status = c("FALSE", "TRUE"), 
  stringsAsFactors = FALSE
)

# create extra columns of data from row data
row_data_fs <- do.call("cbind", lapply(row_data, as.numeric))
stopifnot(all(colnames(row_data_fs) == colnames(row_data)))
stopifnot(nrow(row_data_fs) == nrow(d_exprs[[1]]))

# create marker info
marker_info <- as.data.frame(col_data, row.names = seq_len(nrow(col_data)))

# create flowFrames
ffs <- vector("list", length(d_exprs))
names(ffs) <- names(d_exprs)

for (r in 1:length(d_exprs)) {
  ff <- flowFrame(cbind(d_exprs[[r]], row_data_fs))
  stopifnot(all(colnames(ff) == c(colnames(d_exprs[[r]]), colnames(row_data_fs))))
  # include both channel and marker names in 'pData(parameters(.))'
  stopifnot(length(c(marker_info$marker_name, colnames(row_data_fs))) == nrow(pData(parameters(ff))))
  pData(parameters(ff))$desc <- c(marker_info$marker_name, colnames(row_data_fs))
  
  # include additional information in 'description' slot
  description(ff)$GROUP_INFO <- group_info
  description(ff)$PATIENT_INFO <- patient_info
  description(ff)$SAMPLE_INFO <- sample_info
  description(ff)$SPIKEIN_INFO <- spikein_info
  # experiment information and marker information
  description(ff)$EXPERIMENT_INFO <- meta_data$experiment_info
  description(ff)$MARKER_INFO <- marker_info
  # simulation replicate (seed)
  description(ff)$REPLICATE <- names(d_exprs)[r]
  
  # store
  ffs[[r]] <- ff
}

# create flowSet
d_flowSet <- flowSet(ffs)



# ------------
# Save objects
# ------------

filename_SE <- "Weber_AML_sim_null_SE.rda"
filename_flowSet <- "Weber_AML_sim_null_flowSet.rda"

save(d_SE, file = filename_SE)
save(d_flowSet, file = filename_flowSet)



