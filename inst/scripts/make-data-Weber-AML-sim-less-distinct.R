##########################################################################################
# R script to prepare benchmark dataset 'Weber-AML-sim-less-distinct'
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


# modified to create 'less distinct' spike-in cells: reduce differences in expression
# profiles (medians and standard deviations of arcsinh-transformed values) between
# diseased and healthy blast cells


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



# -------------------------------------------------
# Load data for healthy samples H1-H5 (blast cells)
# -------------------------------------------------

# note sample names and filenames are shuffled
tbl_match_blasts <- read.delim(file_match_samples_blasts)
tbl_match_blasts[grep("H[0-9]+", tbl_match_blasts[, "FCS.Filename"]), c("FCS.Filename", "Individuals")]

files_healthy_blasts <- files_blasts[1:5]

data_healthy_blasts <- lapply(files_healthy_blasts, function(f) exprs(read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)))

names(data_healthy_blasts) <- names(data_healthy)

# check numbers of cells
sapply(data_healthy_blasts, dim)



# ----------------------------------------------------------------------------
# Load data for AML patients: blast cells (CN: patient SJ10; CBF: patient SJ4)
# ----------------------------------------------------------------------------

# note sample names and filenames are shuffled
tbl_match_blasts <- read.delim(file_match_samples_blasts)
tbl_match_blasts[-grep("H[0-9]+", tbl_match_blasts[, "FCS.Filename"]), c("FCS.Filename", "Individuals")]

file_SJ10 <- files_blasts[6]
file_SJ10  # note: sample 'SJ10' has filename 'SJ11'

file_SJ4 <- files_blasts[20]
file_SJ4  # note: sample 'SJ4' has filename 'SJ5'

# load data for SJ10 (CN)
data_SJ10 <- exprs(read.FCS(file_SJ10, transformation = FALSE, truncate_max_range = FALSE))

# load data for SJ4 (CBF)
data_SJ4 <- exprs(read.FCS(file_SJ4, transformation = FALSE, truncate_max_range = FALSE))


# check column names match for all samples (healthy and blasts)
for (i in 1:5) {
  print(all.equal(colnames(data_healthy[[i]]), colnames(data_SJ4)))
}
all.equal(colnames(data_SJ10), colnames(data_SJ4))



# ----------------------
# Delete temporary files
# ----------------------

unlink(DIR_TMP, recursive = TRUE)



# --------------------------------------------
# Replicates: reduced levels of 'distinctness'
# --------------------------------------------

# replicates: create 'less distinct' data sets by reducing difference in median and
# standard deviation of arcsinh-transformed expression by various proportions (e.g. 50%),
# along each dimension (protein marker), for the blast population of interest


# cofactor for arcsinh transform
cofactor <- 5


distinctness <- c(0.5, 0.75)  # 50%, 75%
names(distinctness) <- c("less_50pc", "less_75pc")


data_distinctness <- vector("list", length(distinctness))
names(data_distinctness) <- names(distinctness)


for (di in 1:length(distinctness)) {
  
  data_distinctness[[di]] <- vector("list", 3)
  names(data_distinctness[[di]]) <- c("healthy", "CN", "CBF")
  
  
  
  # ---------------------
  # Split healthy samples
  # ---------------------
  
  # Split each healthy sample (H1-H5) into 3 equal parts. One part will be used as the
  # healthy sample, and the other parts will each have spike-in cells added (for conditions
  # CN and CBF).
  
  data_healthy_base <- data_healthy_CN <- data_healthy_CBF <- 
    vector("list", length(data_healthy))
  names(data_healthy_base) <- names(data_healthy_CN) <- names(data_healthy_CBF) <- 
    names(data_healthy)
  
  # note: use same random seed as in main results (i.e. select same cells for comparability)
  seed <- 100
  
  for (i in 1:length(data_healthy)) {
    data_i <- data_healthy[[i]]
    
    # note: use same random seed as in main results (i.e. select same cells for comparability)
    set.seed(seed + i)
    
    n <- round(nrow(data_i) / 3)
    
    ix_base <- sample(1:nrow(data_i), n)
    ix_CN <- sample((1:nrow(data_i))[-ix_base], n)
    ix_CBF <- setdiff(1:nrow(data_i), c(ix_base, ix_CN))
    
    data_healthy_base[[i]] <- data_i[ix_base, ]
    data_healthy_CN[[i]] <- data_i[ix_CN, ]
    data_healthy_CBF[[i]] <- data_i[ix_CBF, ]
  }
  
  sapply(data_healthy_base, dim)
  sapply(data_healthy_CN, dim)
  sapply(data_healthy_CBF, dim)
  
  
  
  # -------------------------
  # Create spike-in data sets
  # -------------------------
  
  # AML blast cells are subsampled at various thresholds (5%, 1%, 0.1%) of the total number
  # of healthy cells for each sample, and combined with the healthy cells to create the
  # spike-in data sets.
  
  data_spike_CN <- data_spike_CBF <- 
    vector("list", length(data_healthy))
  names(data_spike_CN) <- names(data_spike_CBF) <- 
    names(data_healthy)
  
  thresholds <- c(0.05, 0.01, 0.001)  # 5%, 1%, 0.1%
  thresholds_nm <- paste0(thresholds * 100, "pc")
  
  
  # condition CN (patient SJ10)
  
  data_blasts_AML <- data_SJ10
  cnd <- "CN"
  
  # note: use same random seed as in main results (i.e. select same cells for comparability)
  seed <- 200
  
  for (i in 1:length(data_healthy_CN)) {
    data_i <- data_healthy_CN[[i]]
    nm_i <- names(data_healthy_CN)[i]
    
    data_spike_CN[[i]] <- vector("list", length(thresholds))
    names(data_spike_CN[[i]]) <- thresholds_nm
    
    for (z in 1:length(thresholds)) {
      th <- thresholds[z]
      
      # note: use same random seed as in main results (i.e. select same cells for comparability)
      set.seed(seed + 10 * th + i)
      
      n_spikein <- ceiling(th * nrow(data_i))
      is_spikein <- c(rep(0, nrow(data_i)), rep(1, n_spikein))
      
      cat("n =", n_spikein, "\n")
      
      # subsample blasts
      spikein_i <- data_blasts_AML[sample(1:nrow(data_blasts_AML), n_spikein), , drop = FALSE]
      
      
      # calculate difference in median and standard deviation between AML blasts and
      # healthy blasts along each dimension, then reduce by certain proportion
      
      # note: use arcsinh-transformed values; then convert back to non-transformed
      
      stopifnot(all(colnames(data_blasts_AML) == colnames(data_healthy_blasts[[i]])))
      stopifnot(all(colnames(spikein_i) == colnames(data_healthy_blasts[[i]])))
      
      medians_H <- apply(asinh(data_healthy_blasts[[i]] / cofactor), 2, median)
      medians_AML <- apply(asinh(data_blasts_AML / cofactor), 2, median)
      
      sds_H <- apply(asinh(data_healthy_blasts[[i]] / cofactor), 2, sd)
      sds_AML <- apply(asinh(data_blasts_AML / cofactor), 2, sd)
      
      # use transpose to allow vectorized subtraction
      spikein_i <- t((t(asinh(spikein_i / cofactor)) - (distinctness[di] * (medians_AML - medians_H))) * (1 - (distinctness[di] * (sds_AML - sds_H)) / sds_AML))
      # convert back to non-arcsinh-transformed
      spikein_i <- sinh(spikein_i) * cofactor
      
      
      data_out_i <- rbind(data_i, spikein_i)
      data_out_i <- cbind(data_out_i, spikein = is_spikein)
      
      data_spike_CN[[i]][[z]] <- data_out_i
      
      names(data_spike_CN[[i]][z]) <- paste0(cnd, "_", nm_i, "_", th * 100, "pc")
    }
  }
  
  
  
  # condition CBF (patient SJ4)
  
  data_blasts_AML <- data_SJ4
  cnd <- "CBF"
  
  # note: use same random seed as in main results (i.e. select same cells for comparability)
  seed <- 300
  
  for (i in 1:length(data_healthy_CBF)) {
    data_i <- data_healthy_CBF[[i]]
    nm_i <- names(data_healthy_CBF)[i]
    
    data_spike_CBF[[i]] <- vector("list", length(thresholds))
    names(data_spike_CBF[[i]]) <- thresholds_nm
    
    for (z in 1:length(thresholds)) {
      th <- thresholds[z]
      
      # note: use same random seed as in main results (i.e. select same cells for comparability)
      set.seed(seed + 10 * th + i)
      
      n_spikein <- ceiling(th * nrow(data_i))
      is_spikein <- c(rep(0, nrow(data_i)), rep(1, n_spikein))
      
      cat("n =", n_spikein, "\n")
      
      # subsample blasts
      spikein_i <- data_blasts_AML[sample(1:nrow(data_blasts_AML), n_spikein), , drop = FALSE]
      
      
      # calculate difference in median and standard deviation between AML blasts and
      # healthy blasts along each dimension, then reduce by certain proportion
      
      # note: use arcsinh-transformed values; then convert back to non-transformed
      
      stopifnot(all(colnames(data_blasts_AML) == colnames(data_healthy_blasts[[i]])))
      stopifnot(all(colnames(spikein_i) == colnames(data_healthy_blasts[[i]])))
      
      medians_H <- apply(asinh(data_healthy_blasts[[i]] / cofactor), 2, median)
      medians_AML <- apply(asinh(data_blasts_AML / cofactor), 2, median)
      
      sds_H <- apply(asinh(data_healthy_blasts[[i]] / cofactor), 2, sd)
      sds_AML <- apply(asinh(data_blasts_AML / cofactor), 2, sd)
      
      # use transpose to allow vectorized subtraction
      spikein_i <- t((t(asinh(spikein_i / cofactor)) - (distinctness[di] * (medians_AML - medians_H))) * (1 - (distinctness[di] * (sds_AML - sds_H)) / sds_AML))
      # convert back to non-arcsinh-transformed
      spikein_i <- sinh(spikein_i) * cofactor
      
      
      data_out_i <- rbind(data_i, spikein_i)
      data_out_i <- cbind(data_out_i, spikein = is_spikein)
      
      data_spike_CBF[[i]][[z]] <- data_out_i
      
      names(data_spike_CBF[[i]][z]) <- paste0(cnd, "_", nm_i, "_", th * 100, "pc")
    }
  }
  
  
  # store data
  data_distinctness[[di]][["healthy"]] <- data_healthy_base
  data_distinctness[[di]][["CN"]] <- data_spike_CN
  data_distinctness[[di]][["CBF"]] <- data_spike_CBF
  
}



# ---------------
# Create metadata
# ---------------

# sample information

conditions <- c("healthy", "CN", "CBF")

sample_id <- paste0(rep(conditions, each = length(data_healthy)), "_", names(data_healthy))
# convert to factor with levels in expected order
sample_id <- factor(sample_id, levels = sample_id)

group_id <- factor(gsub("_.*$", "", sample_id), levels = c("healthy", "CN", "CBF"))
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

for (i in 1:length(distinctness)) {
  d_all <- c(data_distinctness[[i]][["healthy"]], data_distinctness[[i]][["CN"]], data_distinctness[[i]][["CBF"]])
  stopifnot(all(sapply(seq_along(d_all), function(i) all(colnames(d_all[[i]]) == colnames(d_all[[1]])))))
}

# channel and marker names
channel_name <- colnames(data_distinctness[[1]][["healthy"]][[1]])
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

# create a separate object for each threshold (simulation), with each object containing
# separate assays for each replicate (distinctness)

d_SE_list <- vector("list", length(thresholds))
names(d_SE_list) <- thresholds_nm

patients_nm <- names(data_healthy)


for (z in 1:length(d_SE_list)) {
  
  # check numbers of cells are identical across replicates
  for (di in 1:length(distinctness)) {
    stopifnot(identical(sapply(data_distinctness[[di]][["healthy"]], dim), sapply(data_distinctness[[1]][["healthy"]], dim)))
    stopifnot(identical(sapply(data_distinctness[[di]][["CN"]][[z]], dim), sapply(data_distinctness[[1]][["CN"]][[z]], dim)))
    stopifnot(identical(sapply(data_distinctness[[di]][["CBF"]][[z]], dim), sapply(data_distinctness[[1]][["CBF"]][[z]], dim)))
  }
  
  
  # set up row data
  n_cells_healthy <- sapply(data_distinctness[[di]][["healthy"]], dim)[1, ]
  n_cells_CN_z <- sapply(patients_nm, function(i) dim(data_distinctness[[di]][["CN"]][[i]][[z]]))[1, ]
  n_cells_CBF_z <- sapply(patients_nm, function(i) dim(data_distinctness[[di]][["CBF"]][[i]][[z]]))[1, ]
  
  n_cells_z <- c(n_cells_healthy, n_cells_CN_z, n_cells_CBF_z)
  stopifnot(length(n_cells_z) == nrow(experiment_info))
  
  stopifnot(all(names(n_cells_z) == gsub("^.*_", "", experiment_info$sample_id)))
  names(n_cells_z) <- experiment_info$sample_id
  
  row_data <- as.data.frame(lapply(experiment_info, function(col) {
    as.factor(rep(col, n_cells_z))
  }))
  
  stopifnot(nrow(row_data) == sum(n_cells_z))
  
  
  # set up column data
  col_data <- marker_info
  
  
  # set up expression data
  # note: one assay per replicate (distinctness)
  d_exprs <- vector("list", length(data_distinctness))
  names(d_exprs) <- names(data_distinctness)
  
  for (di in 1:length(d_exprs)) {
    data_healthy_base_z <- sapply(data_distinctness[[di]][["healthy"]], function(d) {
      cbind(d, spikein = 0)
    })
    
    data_z <- matrix(, nrow = 0, ncol = ncol(data_healthy_base_z[[1]]))
    colnames(data_z) <- colnames(data_healthy_base_z[[1]])
    
    for (i in patients_nm) {
      stopifnot(colnames(data_z) == colnames(data_healthy_base_z[[i]]))
      data_z <- rbind(data_z, data_healthy_base_z[[i]])
    }
    for (i in patients_nm) {
      stopifnot(colnames(data_z) == colnames(data_distinctness[[di]][["CN"]][[i]][[z]]))
      data_z <- rbind(data_z, data_distinctness[[di]][["CN"]][[i]][[z]])
    }
    for (i in patients_nm) {
      stopifnot(colnames(data_z) == colnames(data_distinctness[[di]][["CBF"]][[i]][[z]]))
      data_z <- rbind(data_z, data_distinctness[[di]][["CBF"]][[i]][[z]])
    }
    
    stopifnot(nrow(data_z) == nrow(row_data))
    stopifnot(nrow(data_z) == sum(n_cells_z))
    stopifnot(all(colnames(data_z)[-ncol(data_z)] == marker_info$channel_name))
    
    d_exprs[[di]] <- data_z
  }
  
  
  # move spike-in column to row data
  for (di in 1:length(d_exprs)) {
    stopifnot(identical(d_exprs[[di]][, "spikein", drop = FALSE], d_exprs[[1]][, "spikein", drop = FALSE]))
  }
  
  row_data <- cbind(row_data, spikein = as.logical(d_exprs[[1]][, "spikein"]))
  
  for (di in 1:length(d_exprs)) {
    d_exprs[[di]] <- d_exprs[[di]][, -ncol(d_exprs[[di]])]
    
    stopifnot(all(colnames(d_exprs[[di]]) == marker_info$channel_name))
    stopifnot(ncol(d_exprs[[di]]) == nrow(col_data))
  }
  
  
  # create SummarizedExperiment objects
  d_SE_list[[z]] <- SummarizedExperiment(
    assays = d_exprs, 
    rowData = row_data, 
    colData = col_data, 
    metadata = list(experiment_info = experiment_info, n_cells = n_cells_z)
  )
}



# ----------------------
# Create flowSet objects
# ----------------------

# note: row data is stored as additional columns of data in the expression matrices;
# additional information from row data and column data (e.g. marker classes) is stored in
# 'description' slot

# create a separate flowSet object for each threshold (simulation), with each flowSet
# containing separate flowFrames for each replicate (distinctness)

d_flowSet_list <- vector("list", length(d_SE_list))
names(d_flowSet_list) <- names(d_SE_list)

for (i in seq_along(d_SE_list)) {
  
  # extract data
  row_data <- rowData(d_SE_list[[i]])
  col_data <- colData(d_SE_list[[i]])
  d_exprs <- assays(d_SE_list[[i]])
  meta_data <- metadata(d_SE_list[[i]])
  
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
  stopifnot(nrow(row_data_fs) == nrow(d_exprs))
  
  # create marker info
  marker_info <- as.data.frame(col_data, row.names = seq_len(nrow(col_data)))
  
  # create flowFrames
  ffs <- vector("list", length(d_exprs))
  names(ffs) <- names(d_exprs)
  
  for (di in 1:length(d_exprs)) {
    ff <- flowFrame(cbind(d_exprs[[di]], row_data_fs))
    stopifnot(all(colnames(ff) == c(colnames(d_exprs[[di]]), colnames(row_data_fs))))
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
    description(ff)$REPLICATE <- names(d_exprs)[di]
    
    # store data
    ffs[[di]] <- ff
  }
  
  # create flowSet
  d_flowSet_list[[i]] <- flowSet(ffs)
}



# ------------
# Save objects
# ------------

stopifnot(all(names(d_SE_list) == names(d_flowSet_list)))

filenames_SE <- paste0("Weber_AML_sim_less_distinct_", names(d_SE_list), "_SE.rda")
filenames_flowSet <- paste0("Weber_AML_sim_less_distinct_", names(d_flowSet_list), "_flowSet.rda")

for (i in 1:length(d_SE_list)) {
  d_SE <- d_SE_list[[i]]
  d_flowSet <- d_flowSet_list[[i]]
  
  save(d_SE, file = filenames_SE[i])
  save(d_flowSet, file = filenames_flowSet[i])
}



