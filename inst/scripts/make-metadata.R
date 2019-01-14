# ------------------------------
# Base metadata for all datasets
# ------------------------------

df_base <- data.frame(
  BiocVersion = "3.8", 
  Genome = NA, 
  SourceType = "FCS", 
  SourceUrl = "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/", 
  SourceVersion = NA, 
  Coordinate_1_based = TRUE, 
  DataProvider = "", 
  Maintainer = "Lukas M. Weber <lukmweber@gmail.com>", 
  DispatchClass = "Rda", 
  stringsAsFactors = FALSE
)


# --------------------
# Levine_32dim dataset
# --------------------

# SummarizedExperiment
df_Levine_32dim_SE <- cbind(
  df_base, 
  Title = "Levine_32dim_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Levine et al. (2015), ", 
    "containing 32 dimensions (surface protein markers). ", 
    "Manually gated cell population labels are available for 14 populations. ", 
    "Cells are human bone marrow cells from 2 healthy donors. ", 
    "This dataset can be used to benchmark clustering algorithms."), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Levine_32dim/Levine_32dim_SE.rda"
)

# flowSet
df_Levine_32dim_flowSet <- df_Levine_32dim_SE
df_Levine_32dim_flowSet$Title <- "Levine_32dim_flowSet"
df_Levine_32dim_flowSet$RDataClass <- "flowSet"
df_Levine_32dim_flowSet$RDataPath <- "HDCytoData/Levine_32dim/Levine_32dim_flowSet.rda"



# --------------------
# Levine_13dim dataset
# --------------------

# SummarizedExperiment
df_Levine_13dim_SE <- cbind(
  df_base, 
  Title = "Levine_13dim_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Levine et al. (2015), ", 
    "containing 13 dimensions (surface protein markers). ", 
    "Manually gated cell population labels are available for 24 populations. ", 
    "Cells are human bone marrow cells from a single healthy donor. ", 
    "This dataset can be used to benchmark clustering algorithms."), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Levine_13dim/Levine_13dim_SE.rda"
)

# flowSet
df_Levine_13dim_flowSet <- df_Levine_13dim_SE
df_Levine_13dim_flowSet$Title <- "Levine_13dim_flowSet"
df_Levine_13dim_flowSet$RDataClass <- "flowSet"
df_Levine_13dim_flowSet$RDataPath <- "HDCytoData/Levine_13dim/Levine_13dim_flowSet.rda"



# ------------------
# Samusik_01 dataset
# ------------------

# SummarizedExperiment
df_Samusik_01_SE <- cbind(
  df_base, 
  Title = "Samusik_01_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Samusik et al. (2016), ", 
    "containing 39 dimensions (surface protein markers). ", 
    "Manually gated cell population labels are available for 24 populations. ", 
    "The full dataset (Samusik_all) contains cells from 10 replicate bone marrow samples from ", 
    "C57BL/6J mice (i.e. samples from 10 different mice); this dataset (Samusik_01) contains ", 
    "the data from sample 01 only. ", 
    "This dataset can be used to benchmark clustering algorithms."), 
  Species = "Mus musculus", 
  TaxonomyId = "10090", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Samusik_01/Samusik_01_SE.rda"
)

# flowSet
df_Samusik_01_flowSet <- df_Samusik_01_SE
df_Samusik_01_flowSet$Title <- "Samusik_01_flowSet"
df_Samusik_01_flowSet$RDataClass <- "flowSet"
df_Samusik_01_flowSet$RDataPath <- "HDCytoData/Samusik_01/Samusik_01_flowSet.rda"



# -------------------
# Samusik_all dataset
# -------------------

# SummarizedExperiment
df_Samusik_all_SE <- cbind(
  df_base, 
  Title = "Samusik_all_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Samusik et al. (2016), ", 
    "containing 39 dimensions (surface protein markers). ", 
    "Manually gated cell population labels are available for 24 populations. ", 
    "This dataset contains cells from 10 replicate bone marrow samples from ", 
    "C57BL/6J mice (i.e. samples from 10 different mice). ", 
    "This dataset can be used to benchmark clustering algorithms."), 
  Species = "Mus musculus", 
  TaxonomyId = "10090", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Samusik_all/Samusik_all_SE.rda"
)

# flowSet
df_Samusik_all_flowSet <- df_Samusik_all_SE
df_Samusik_all_flowSet$Title <- "Samusik_all_flowSet"
df_Samusik_all_flowSet$RDataClass <- "flowSet"
df_Samusik_all_flowSet$RDataPath <- "HDCytoData/Samusik_all/Samusik_all_flowSet.rda"



# --------------------
# Nilsson_rare dataset
# --------------------

# SummarizedExperiment
df_Nilsson_rare_SE <- cbind(
  df_base, 
  Title = "Nilsson_rare_SE", 
  Description = paste0(
    "Flow cytometry dataset from Nilsson et al. (2013), ", 
    "containing 13 dimensions (surface protein markers). ", 
    "Manually gated cell population labels are available for one rare population of ", 
    "hematopoietic stem cells (HSCs). ", 
    "Cells are human bone marrow cells from a single healthy donor. ", 
    "This dataset can be used to benchmark clustering algorithms for rare cell populations."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Nilsson_rare/Nilsson_rare_SE.rda"
)

# flowSet
df_Nilsson_rare_flowSet <- df_Nilsson_rare_SE
df_Nilsson_rare_flowSet$Title <- "Nilsson_rare_flowSet"
df_Nilsson_rare_flowSet$RDataClass <- "flowSet"
df_Nilsson_rare_flowSet$RDataPath <- "HDCytoData/Nilsson_rare/Nilsson_rare_flowSet.rda"



# --------------------
# Mosmann_rare dataset
# --------------------

# SummarizedExperiment
df_Mosmann_rare_SE <- cbind(
  df_base, 
  Title = "Mosmann_rare_SE", 
  Description = paste0(
    "Flow cytometry dataset from Mosmann et al. (2014), ", 
    "containing 14 dimensions (7 surface protein markers and 7 signaling markers). ", 
    "Manually gated cell population labels are available for one rare population of ", 
    "activated (cytokine-producing) memory CD4 T cells. ", 
    "Cells are human peripheral blood cells exposed to influenza antigens, from a single healthy donor. ", 
    "This dataset can be used to benchmark clustering algorithms for rare cell populations."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Mosmann_rare/Mosmann_rare_SE.rda"
)

# flowSet
df_Mosmann_rare_flowSet <- df_Mosmann_rare_SE
df_Mosmann_rare_flowSet$Title <- "Mosmann_rare_flowSet"
df_Mosmann_rare_flowSet$RDataClass <- "flowSet"
df_Mosmann_rare_flowSet$RDataPath <- "HDCytoData/Mosmann_rare/Mosmann_rare_flowSet.rda"



# -----------------------
# Krieg_Anti_PD_1 dataset
# -----------------------

# SummarizedExperiment
df_Krieg_Anti_PD_1_SE <- cbind(
  df_base, 
  Title = "Krieg_Anti_PD_1_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Krieg et al. (2018), consisting of 20 baseline samples ", 
    "(prior to treatment) of peripheral blood from melanoma skin cancer patients subsequently treated ", 
    "with anti-PD-1 immunotherapy. The samples are split across 2 conditions (non-responders and ", 
    "responders) and 2 batches. This dataset can be used to benchmark algorithms for differential ", 
    "analysis, in particular detection of differentially abundant rare cell populations."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Krieg_Anti_PD_1/Krieg_Anti_PD_1_SE.rda"
)

# flowSet
df_Krieg_Anti_PD_1_flowSet <- df_Krieg_Anti_PD_1_SE
df_Krieg_Anti_PD_1_flowSet$Title <- "Krieg_Anti_PD_1_flowSet"
df_Krieg_Anti_PD_1_flowSet$RDataClass <- "flowSet"
df_Krieg_Anti_PD_1_flowSet$RDataPath <- "HDCytoData/Krieg_Anti_PD_1/Krieg_Anti_PD_1_flowSet.rda"



# --------------------------
# Bodenmiller_BCR_XL dataset
# --------------------------

# SummarizedExperiment
df_Bodenmiller_BCR_XL_SE <- cbind(
  df_base, 
  Title = "Bodenmiller_BCR_XL_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Bodenmiller et al. (2012), consisting of 8 paired samples ", 
    "(16 samples) of stimulated (BCR-XL) and unstimulated peripheral blood cells from healthy ", 
    "individuals. This dataset can be used to benchmark algorithms for differential analysis, ", 
    "in particular detection of differential states within cell populations."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Bodenmiller_BCR_XL/Bodenmiller_BCR_XL_SE.rda"
)

# flowSet
df_Bodenmiller_BCR_XL_flowSet <- df_Bodenmiller_BCR_XL_SE
df_Bodenmiller_BCR_XL_flowSet$Title <- "Bodenmiller_BCR_XL_flowSet"
df_Bodenmiller_BCR_XL_flowSet$RDataClass <- "flowSet"
df_Bodenmiller_BCR_XL_flowSet$RDataPath <- "HDCytoData/Bodenmiller_BCR_XL/Bodenmiller_BCR_XL_flowSet.rda"


# ---------------------------------
# Combine for all datasets and save
# ---------------------------------

# Combine all datasets
df_all <- rbind(
  df_Levine_32dim_SE, 
  df_Levine_32dim_flowSet, 
  df_Levine_13dim_SE, 
  df_Levine_13dim_flowSet, 
  df_Samusik_01_SE, 
  df_Samusik_01_flowSet, 
  df_Samusik_all_SE, 
  df_Samusik_all_flowSet, 
  df_Nilsson_rare_SE, 
  df_Nilsson_rare_flowSet, 
  df_Mosmann_rare_SE, 
  df_Mosmann_rare_flowSet, 
  df_Krieg_Anti_PD_1_SE, 
  df_Krieg_Anti_PD_1_flowSet, 
  df_Bodenmiller_BCR_XL_SE, 
  df_Bodenmiller_BCR_XL_flowSet
)

# Save as .csv file
write.csv(df_all, file = "../extdata/metadata.csv", row.names = FALSE)


