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


# --------------------------
# Bodenmiller_BCR_XL dataset
# --------------------------

# SummarizedExperiment
df_Bodenmiller_BCR_XL_SE <- cbind(
  df_base, 
  Title = "Bodenmiller_BCR_XL_SE", 
  Description = paste0(
    "Mass cytometry data from Bodenmiller et al. (2012). Paired samples of healthy PBMCs, ", 
    "where one sample from each pair was stimulated with B cell receptor / Fc receptor ", 
    "cross-linker (BCR-XL). This dataset can be used to benchmark algorithms for differential ", 
    "analysis."), 
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


# --------------------
# Levine_32dim dataset
# --------------------

# SummarizedExperiment
df_Levine_32dim_SE <- cbind(
  df_base, 
  Title = "Levine_32dim_SE", 
  Description = paste0(
    "32-dimensional mass cytometry (CyTOF) dataset from Levine et al. (2015), containing 14 manually ", 
    "gated immune cell populations from healthy human bone marrow mononuclear cells (BMMCs). This ", 
    "dataset can be used to benchmark clustering algorithms."), 
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
    "13-dimensional mass cytometry (CyTOF) dataset from Levine et al. (2015), containing 24 manually ", 
    "gated immune cell populations from healthy human bone marrow mononuclear cells (BMMCs). This ", 
    "dataset can be used to benchmark clustering algorithms."), 
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



# ---------------------------------
# Combine for all datasets and save
# ---------------------------------

# Combine all datasets
df_all <- rbind(
  df_Bodenmiller_BCR_XL_SE, 
  df_Bodenmiller_BCR_XL_flowSet, 
  df_Levine_32dim_SE, 
  df_Levine_32dim_flowSet, 
  df_Levine_13dim_SE, 
  df_Levine_13dim_flowSet
)

# Save as .csv file
write.csv(df_all, file = "../extdata/metadata.csv", row.names = FALSE)


