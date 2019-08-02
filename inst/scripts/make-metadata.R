# ------------------------------
# Base metadata for all datasets
# ------------------------------

df_base <- data.frame(
  BiocVersion = "3.10", 
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


# ------------
# Levine_32dim
# ------------

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
  RDataPath = "HDCytoData/Levine_32dim/Levine_32dim_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Levine_32dim_flowSet <- df_Levine_32dim_SE
df_Levine_32dim_flowSet$Title <- "Levine_32dim_flowSet"
df_Levine_32dim_flowSet$RDataClass <- "flowSet"
df_Levine_32dim_flowSet$RDataPath <- "HDCytoData/Levine_32dim/Levine_32dim_flowSet.rda"



# ------------
# Levine_13dim
# ------------

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
  RDataPath = "HDCytoData/Levine_13dim/Levine_13dim_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Levine_13dim_flowSet <- df_Levine_13dim_SE
df_Levine_13dim_flowSet$Title <- "Levine_13dim_flowSet"
df_Levine_13dim_flowSet$RDataClass <- "flowSet"
df_Levine_13dim_flowSet$RDataPath <- "HDCytoData/Levine_13dim/Levine_13dim_flowSet.rda"



# ----------
# Samusik_01
# ----------

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
  RDataPath = "HDCytoData/Samusik_01/Samusik_01_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Samusik_01_flowSet <- df_Samusik_01_SE
df_Samusik_01_flowSet$Title <- "Samusik_01_flowSet"
df_Samusik_01_flowSet$RDataClass <- "flowSet"
df_Samusik_01_flowSet$RDataPath <- "HDCytoData/Samusik_01/Samusik_01_flowSet.rda"



# -----------
# Samusik_all
# -----------

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
  RDataPath = "HDCytoData/Samusik_all/Samusik_all_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Samusik_all_flowSet <- df_Samusik_all_SE
df_Samusik_all_flowSet$Title <- "Samusik_all_flowSet"
df_Samusik_all_flowSet$RDataClass <- "flowSet"
df_Samusik_all_flowSet$RDataPath <- "HDCytoData/Samusik_all/Samusik_all_flowSet.rda"



# ------------
# Nilsson_rare
# ------------

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
  RDataPath = "HDCytoData/Nilsson_rare/Nilsson_rare_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Nilsson_rare_flowSet <- df_Nilsson_rare_SE
df_Nilsson_rare_flowSet$Title <- "Nilsson_rare_flowSet"
df_Nilsson_rare_flowSet$RDataClass <- "flowSet"
df_Nilsson_rare_flowSet$RDataPath <- "HDCytoData/Nilsson_rare/Nilsson_rare_flowSet.rda"



# ------------
# Mosmann_rare
# ------------

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
  RDataPath = "HDCytoData/Mosmann_rare/Mosmann_rare_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Mosmann_rare_flowSet <- df_Mosmann_rare_SE
df_Mosmann_rare_flowSet$Title <- "Mosmann_rare_flowSet"
df_Mosmann_rare_flowSet$RDataClass <- "flowSet"
df_Mosmann_rare_flowSet$RDataPath <- "HDCytoData/Mosmann_rare/Mosmann_rare_flowSet.rda"



# ---------------
# Krieg_Anti_PD_1
# ---------------

# SummarizedExperiment
df_Krieg_Anti_PD_1_SE <- cbind(
  df_base, 
  Title = "Krieg_Anti_PD_1_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Krieg et al. (2018), consisting of 20 baseline samples ", 
    "(prior to treatment) of peripheral blood from melanoma skin cancer patients subsequently treated ", 
    "with anti-PD-1 immunotherapy. The samples are split across 2 conditions (non-responders and ", 
    "responders) and 2 batches. This dataset can be used to benchmark differential analysis algorithms ", 
    "used to test for differentially abundant rare cell populations."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Krieg_Anti_PD_1/Krieg_Anti_PD_1_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Krieg_Anti_PD_1_flowSet <- df_Krieg_Anti_PD_1_SE
df_Krieg_Anti_PD_1_flowSet$Title <- "Krieg_Anti_PD_1_flowSet"
df_Krieg_Anti_PD_1_flowSet$RDataClass <- "flowSet"
df_Krieg_Anti_PD_1_flowSet$RDataPath <- "HDCytoData/Krieg_Anti_PD_1/Krieg_Anti_PD_1_flowSet.rda"



# ------------------
# Bodenmiller_BCR_XL
# ------------------

# SummarizedExperiment
df_Bodenmiller_BCR_XL_SE <- cbind(
  df_base, 
  Title = "Bodenmiller_BCR_XL_SE", 
  Description = paste0(
    "Mass cytometry (CyTOF) dataset from Bodenmiller et al. (2012), consisting of 8 paired samples ", 
    "(16 samples) of stimulated (BCR-XL) and unstimulated peripheral blood cells from healthy ", 
    "individuals. This dataset can be used to benchmark differential analysis algorithms ", 
    "used to test for differential states within cell populations."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Bodenmiller_BCR_XL/Bodenmiller_BCR_XL_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Bodenmiller_BCR_XL_flowSet <- df_Bodenmiller_BCR_XL_SE
df_Bodenmiller_BCR_XL_flowSet$Title <- "Bodenmiller_BCR_XL_flowSet"
df_Bodenmiller_BCR_XL_flowSet$RDataClass <- "flowSet"
df_Bodenmiller_BCR_XL_flowSet$RDataPath <- "HDCytoData/Bodenmiller_BCR_XL/Bodenmiller_BCR_XL_flowSet.rda"



# ------------------
# Weber_AML_sim_main
# ------------------

# main simulations

# SummarizedExperiments
df_Weber_AML_sim_main_5pc_SE <- cbind(
  df_base, 
  Title = "Weber_AML_sim_main_5pc_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), ", 
    "constructed by computationally 'spiking in' small percentages of AML (acute myeloid leukemia) ", 
    "blast cells into samples of healthy BMMCs (bone marrow mononuclear cells), simulating ", 
    "the phenotype of minimal residual disease (MRD) in AML patients. ", 
    "Main simulations. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differentially abundant rare cell populations. ", 
    "Raw data sourced from Levine et al. (2015), and data generation strategy modified from ", 
    "Arvaniti et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_AML_sim/Weber_AML_sim_main_5pc_SE.rda", 
  stringsAsFactors = FALSE
)

# additional thresholds
df_Weber_AML_sim_main_1pc_SE <- df_Weber_AML_sim_main_5pc_SE
df_Weber_AML_sim_main_1pc_SE$Title <- gsub("5pc", "1pc", df_Weber_AML_sim_main_1pc_SE$Title)
df_Weber_AML_sim_main_1pc_SE$RDataPath <- gsub("5pc", "1pc", df_Weber_AML_sim_main_1pc_SE$RDataPath)

df_Weber_AML_sim_main_0.1pc_SE <- df_Weber_AML_sim_main_5pc_SE
df_Weber_AML_sim_main_0.1pc_SE$Title <- gsub("5pc", "0.1pc", df_Weber_AML_sim_main_0.1pc_SE$Title)
df_Weber_AML_sim_main_0.1pc_SE$RDataPath <- gsub("5pc", "0.1pc", df_Weber_AML_sim_main_0.1pc_SE$RDataPath)

# additional object for all blasts
df_Weber_AML_sim_main_blasts_all_SE <- df_Weber_AML_sim_main_5pc_SE
df_Weber_AML_sim_main_blasts_all_SE$Title <- gsub("5pc", "blasts_all", df_Weber_AML_sim_main_blasts_all_SE$Title)
df_Weber_AML_sim_main_blasts_all_SE$RDataPath <- gsub("5pc", "blasts_all", df_Weber_AML_sim_main_blasts_all_SE$RDataPath)

# flowSets
df_Weber_AML_sim_main_5pc_flowSet <- df_Weber_AML_sim_main_5pc_SE
df_Weber_AML_sim_main_5pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_main_5pc_flowSet$Title)
df_Weber_AML_sim_main_5pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_main_5pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_main_5pc_flowSet$RDataPath)

df_Weber_AML_sim_main_1pc_flowSet <- df_Weber_AML_sim_main_1pc_SE
df_Weber_AML_sim_main_1pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_main_1pc_flowSet$Title)
df_Weber_AML_sim_main_1pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_main_1pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_main_1pc_flowSet$RDataPath)

df_Weber_AML_sim_main_0.1pc_flowSet <- df_Weber_AML_sim_main_0.1pc_SE
df_Weber_AML_sim_main_0.1pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_main_0.1pc_flowSet$Title)
df_Weber_AML_sim_main_0.1pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_main_0.1pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_main_0.1pc_flowSet$RDataPath)

df_Weber_AML_sim_main_blasts_all_flowSet <- df_Weber_AML_sim_main_blasts_all_SE
df_Weber_AML_sim_main_blasts_all_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_main_blasts_all_flowSet$Title)
df_Weber_AML_sim_main_blasts_all_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_main_blasts_all_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_main_blasts_all_flowSet$RDataPath)



# ------------------
# Weber_AML_sim_null
# ------------------

# additional simulations: null simulations

# SummarizedExperiments
df_Weber_AML_sim_null_SE <- cbind(
  df_base, 
  Title = "Weber_AML_sim_null_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), ", 
    "constructed by computationally 'spiking in' small percentages of AML (acute myeloid leukemia) ", 
    "blast cells into samples of healthy BMMCs (bone marrow mononuclear cells), simulating ", 
    "the phenotype of minimal residual disease (MRD) in AML patients. ", 
    "Additional simulations: null simulations. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differentially abundant rare cell populations. ", 
    "Raw data sourced from Levine et al. (2015), and data generation strategy modified from ", 
    "Arvaniti et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_AML_sim/Weber_AML_sim_null_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSets
df_Weber_AML_sim_null_flowSet <- df_Weber_AML_sim_null_SE
df_Weber_AML_sim_null_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_null_flowSet$Title)
df_Weber_AML_sim_null_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_null_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_null_flowSet$RDataPath)



# --------------------------
# Weber_AML_sim_random_seeds
# --------------------------

# additional simulations: modified random seeds

# SummarizedExperiments
df_Weber_AML_sim_random_seeds_5pc_SE <- cbind(
  df_base, 
  Title = "Weber_AML_sim_random_seeds_5pc_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), ", 
    "constructed by computationally 'spiking in' small percentages of AML (acute myeloid leukemia) ", 
    "blast cells into samples of healthy BMMCs (bone marrow mononuclear cells), simulating ", 
    "the phenotype of minimal residual disease (MRD) in AML patients. ", 
    "Additional simulations: modified random seeds. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differentially abundant rare cell populations. ", 
    "Raw data sourced from Levine et al. (2015), and data generation strategy modified from ", 
    "Arvaniti et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_AML_sim/Weber_AML_sim_random_seeds_5pc_SE.rda", 
  stringsAsFactors = FALSE
)

# additional thresholds
df_Weber_AML_sim_random_seeds_1pc_SE <- df_Weber_AML_sim_random_seeds_5pc_SE
df_Weber_AML_sim_random_seeds_1pc_SE$Title <- gsub("5pc", "1pc", df_Weber_AML_sim_random_seeds_1pc_SE$Title)
df_Weber_AML_sim_random_seeds_1pc_SE$RDataPath <- gsub("5pc", "1pc", df_Weber_AML_sim_random_seeds_1pc_SE$RDataPath)

df_Weber_AML_sim_random_seeds_0.1pc_SE <- df_Weber_AML_sim_random_seeds_5pc_SE
df_Weber_AML_sim_random_seeds_0.1pc_SE$Title <- gsub("5pc", "0.1pc", df_Weber_AML_sim_random_seeds_0.1pc_SE$Title)
df_Weber_AML_sim_random_seeds_0.1pc_SE$RDataPath <- gsub("5pc", "0.1pc", df_Weber_AML_sim_random_seeds_0.1pc_SE$RDataPath)

# flowSets
df_Weber_AML_sim_random_seeds_5pc_flowSet <- df_Weber_AML_sim_random_seeds_5pc_SE
df_Weber_AML_sim_random_seeds_5pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_random_seeds_5pc_flowSet$Title)
df_Weber_AML_sim_random_seeds_5pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_random_seeds_5pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_random_seeds_5pc_flowSet$RDataPath)

df_Weber_AML_sim_random_seeds_1pc_flowSet <- df_Weber_AML_sim_random_seeds_1pc_SE
df_Weber_AML_sim_random_seeds_1pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_random_seeds_1pc_flowSet$Title)
df_Weber_AML_sim_random_seeds_1pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_random_seeds_1pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_random_seeds_1pc_flowSet$RDataPath)

df_Weber_AML_sim_random_seeds_0.1pc_flowSet <- df_Weber_AML_sim_random_seeds_0.1pc_SE
df_Weber_AML_sim_random_seeds_0.1pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_random_seeds_0.1pc_flowSet$Title)
df_Weber_AML_sim_random_seeds_0.1pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_random_seeds_0.1pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_random_seeds_0.1pc_flowSet$RDataPath)



# ---------------------------
# Weber_AML_sim_less_distinct
# ---------------------------

# additional simulations: 'less distinct' spike-in cells

# SummarizedExperiments
df_Weber_AML_sim_less_distinct_5pc_SE <- cbind(
  df_base, 
  Title = "Weber_AML_sim_less_distinct_5pc_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), ", 
    "constructed by computationally 'spiking in' small percentages of AML (acute myeloid leukemia) ", 
    "blast cells into samples of healthy BMMCs (bone marrow mononuclear cells), simulating ", 
    "the phenotype of minimal residual disease (MRD) in AML patients. ", 
    "Additional simulations: 'less distinct' spike-in cells. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differentially abundant rare cell populations. ", 
    "Raw data sourced from Levine et al. (2015), and data generation strategy modified from ", 
    "Arvaniti et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_AML_sim/Weber_AML_sim_less_distinct_5pc_SE.rda", 
  stringsAsFactors = FALSE
)

# additional thresholds
df_Weber_AML_sim_less_distinct_1pc_SE <- df_Weber_AML_sim_less_distinct_5pc_SE
df_Weber_AML_sim_less_distinct_1pc_SE$Title <- gsub("5pc", "1pc", df_Weber_AML_sim_less_distinct_1pc_SE$Title)
df_Weber_AML_sim_less_distinct_1pc_SE$RDataPath <- gsub("5pc", "1pc", df_Weber_AML_sim_less_distinct_1pc_SE$RDataPath)

df_Weber_AML_sim_less_distinct_0.1pc_SE <- df_Weber_AML_sim_less_distinct_5pc_SE
df_Weber_AML_sim_less_distinct_0.1pc_SE$Title <- gsub("5pc", "0.1pc", df_Weber_AML_sim_less_distinct_0.1pc_SE$Title)
df_Weber_AML_sim_less_distinct_0.1pc_SE$RDataPath <- gsub("5pc", "0.1pc", df_Weber_AML_sim_less_distinct_0.1pc_SE$RDataPath)

# flowSets
df_Weber_AML_sim_less_distinct_5pc_flowSet <- df_Weber_AML_sim_less_distinct_5pc_SE
df_Weber_AML_sim_less_distinct_5pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_less_distinct_5pc_flowSet$Title)
df_Weber_AML_sim_less_distinct_5pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_less_distinct_5pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_less_distinct_5pc_flowSet$RDataPath)

df_Weber_AML_sim_less_distinct_1pc_flowSet <- df_Weber_AML_sim_less_distinct_1pc_SE
df_Weber_AML_sim_less_distinct_1pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_less_distinct_1pc_flowSet$Title)
df_Weber_AML_sim_less_distinct_1pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_less_distinct_1pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_less_distinct_1pc_flowSet$RDataPath)

df_Weber_AML_sim_less_distinct_0.1pc_flowSet <- df_Weber_AML_sim_less_distinct_0.1pc_SE
df_Weber_AML_sim_less_distinct_0.1pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_AML_sim_less_distinct_0.1pc_flowSet$Title)
df_Weber_AML_sim_less_distinct_0.1pc_flowSet$RDataClass <- "flowSet"
df_Weber_AML_sim_less_distinct_0.1pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_AML_sim_less_distinct_0.1pc_flowSet$RDataPath)



# ---------------------
# Weber_BCR_XL_sim_main
# ---------------------

# main simulation

# SummarizedExperiment
df_Weber_BCR_XL_sim_main_SE <- cbind(
  df_base, 
  Title = "Weber_BCR_XL_sim_main_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), constructed by ", 
    "randomly splitting unstimulated (reference) samples of PBMCs (peripheral blood mononuclear cells) ", 
    "into two halves, and replacing B cells in one half with stimulated (BCR-XL) B cells from ", 
    "corresponding paired samples. ", 
    "Main simulation. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differential states within cell populations. ", 
    "Raw data sourced from Bodenmiller et al. (2012); cell population labels reproduced from ", 
    "Nowicka et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_BCR_XL_sim/Weber_BCR_XL_sim_main_SE.rda", 
  stringsAsFactors = FALSE
)

# flowSet
df_Weber_BCR_XL_sim_main_flowSet <- df_Weber_BCR_XL_sim_main_SE
df_Weber_BCR_XL_sim_main_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_main_flowSet$Title)
df_Weber_BCR_XL_sim_main_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_main_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_main_flowSet$RDataPath)



# ---------------------
# Weber_BCR_XL_sim_null
# ---------------------

# additional simulations: null simulations

# SummarizedExperiments
df_Weber_BCR_XL_sim_null_rep1_SE <- cbind(
  df_base, 
  Title = "Weber_BCR_XL_sim_null_rep1_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), constructed by ", 
    "randomly splitting unstimulated (reference) samples of PBMCs (peripheral blood mononuclear cells) ", 
    "into two halves, and replacing B cells in one half with stimulated (BCR-XL) B cells from ", 
    "corresponding paired samples. ", 
    "Additional simulations: null simulations. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differential states within cell populations. ", 
    "Raw data sourced from Bodenmiller et al. (2012); cell population labels reproduced from ", 
    "Nowicka et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_BCR_XL_sim/Weber_BCR_XL_sim_null_rep1_SE.rda", 
  stringsAsFactors = FALSE
)

# additional replicates
df_Weber_BCR_XL_sim_null_rep2_SE <- df_Weber_BCR_XL_sim_null_rep1_SE
df_Weber_BCR_XL_sim_null_rep2_SE$Title <- gsub("rep1", "rep2", df_Weber_BCR_XL_sim_null_rep2_SE$Title)
df_Weber_BCR_XL_sim_null_rep2_SE$RDataPath <- gsub("rep1", "rep2", df_Weber_BCR_XL_sim_null_rep2_SE$RDataPath)

df_Weber_BCR_XL_sim_null_rep3_SE <- df_Weber_BCR_XL_sim_null_rep1_SE
df_Weber_BCR_XL_sim_null_rep3_SE$Title <- gsub("rep1", "rep3", df_Weber_BCR_XL_sim_null_rep3_SE$Title)
df_Weber_BCR_XL_sim_null_rep3_SE$RDataPath <- gsub("rep1", "rep3", df_Weber_BCR_XL_sim_null_rep3_SE$RDataPath)

# flowSets
df_Weber_BCR_XL_sim_null_rep1_flowSet <- df_Weber_BCR_XL_sim_null_rep1_SE
df_Weber_BCR_XL_sim_null_rep1_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_null_rep1_flowSet$Title)
df_Weber_BCR_XL_sim_null_rep1_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_null_rep1_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_null_rep1_flowSet$RDataPath)

df_Weber_BCR_XL_sim_null_rep2_flowSet <- df_Weber_BCR_XL_sim_null_rep2_SE
df_Weber_BCR_XL_sim_null_rep2_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_null_rep2_flowSet$Title)
df_Weber_BCR_XL_sim_null_rep2_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_null_rep2_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_null_rep2_flowSet$RDataPath)

df_Weber_BCR_XL_sim_null_rep3_flowSet <- df_Weber_BCR_XL_sim_null_rep3_SE
df_Weber_BCR_XL_sim_null_rep3_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_null_rep3_flowSet$Title)
df_Weber_BCR_XL_sim_null_rep3_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_null_rep3_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_null_rep3_flowSet$RDataPath)



# -----------------------------
# Weber_BCR_XL_sim_random_seeds
# -----------------------------

# additional simulations: modified random seeds

# SummarizedExperiments
df_Weber_BCR_XL_sim_random_seeds_rep1_SE <- cbind(
  df_base, 
  Title = "Weber_BCR_XL_sim_random_seeds_rep1_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), constructed by ", 
    "randomly splitting unstimulated (reference) samples of PBMCs (peripheral blood mononuclear cells) ", 
    "into two halves, and replacing B cells in one half with stimulated (BCR-XL) B cells from ", 
    "corresponding paired samples. ", 
    "Additional simulations: modified random seeds. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differential states within cell populations. ", 
    "Raw data sourced from Bodenmiller et al. (2012); cell population labels reproduced from ", 
    "Nowicka et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_BCR_XL_sim/Weber_BCR_XL_sim_random_seeds_rep1_SE.rda", 
  stringsAsFactors = FALSE
)

# additional replicates
df_Weber_BCR_XL_sim_random_seeds_rep2_SE <- df_Weber_BCR_XL_sim_random_seeds_rep1_SE
df_Weber_BCR_XL_sim_random_seeds_rep2_SE$Title <- gsub("rep1", "rep2", df_Weber_BCR_XL_sim_random_seeds_rep2_SE$Title)
df_Weber_BCR_XL_sim_random_seeds_rep2_SE$RDataPath <- gsub("rep1", "rep2", df_Weber_BCR_XL_sim_random_seeds_rep2_SE$RDataPath)

df_Weber_BCR_XL_sim_random_seeds_rep3_SE <- df_Weber_BCR_XL_sim_random_seeds_rep1_SE
df_Weber_BCR_XL_sim_random_seeds_rep3_SE$Title <- gsub("rep1", "rep3", df_Weber_BCR_XL_sim_random_seeds_rep3_SE$Title)
df_Weber_BCR_XL_sim_random_seeds_rep3_SE$RDataPath <- gsub("rep1", "rep3", df_Weber_BCR_XL_sim_random_seeds_rep3_SE$RDataPath)

# flowSets
df_Weber_BCR_XL_sim_random_seeds_rep1_flowSet <- df_Weber_BCR_XL_sim_random_seeds_rep1_SE
df_Weber_BCR_XL_sim_random_seeds_rep1_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_random_seeds_rep1_flowSet$Title)
df_Weber_BCR_XL_sim_random_seeds_rep1_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_random_seeds_rep1_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_random_seeds_rep1_flowSet$RDataPath)

df_Weber_BCR_XL_sim_random_seeds_rep2_flowSet <- df_Weber_BCR_XL_sim_random_seeds_rep2_SE
df_Weber_BCR_XL_sim_random_seeds_rep2_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_random_seeds_rep2_flowSet$Title)
df_Weber_BCR_XL_sim_random_seeds_rep2_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_random_seeds_rep2_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_random_seeds_rep2_flowSet$RDataPath)

df_Weber_BCR_XL_sim_random_seeds_rep3_flowSet <- df_Weber_BCR_XL_sim_random_seeds_rep3_SE
df_Weber_BCR_XL_sim_random_seeds_rep3_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_random_seeds_rep3_flowSet$Title)
df_Weber_BCR_XL_sim_random_seeds_rep3_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_random_seeds_rep3_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_random_seeds_rep3_flowSet$RDataPath)



# ------------------------------
# Weber_BCR_XL_sim_less_distinct
# ------------------------------

# additional simulations: 'less distinct' spike-in cells

# SummarizedExperiments
df_Weber_BCR_XL_sim_less_distinct_less_50pc_SE <- cbind(
  df_base, 
  Title = "Weber_BCR_XL_sim_less_distinct_less_50pc_SE", 
  Description = paste0(
    "Semi-simulated mass cytometry (CyTOF) dataset from Weber et al. (2019), constructed by ", 
    "randomly splitting unstimulated (reference) samples of PBMCs (peripheral blood mononuclear cells) ", 
    "into two halves, and replacing B cells in one half with stimulated (BCR-XL) B cells from ", 
    "corresponding paired samples. ", 
    "Additional simulations: 'less distinct' spike-in cells. ", 
    "This dataset can be used to benchmark differential analysis algorithms used to test for ", 
    "differential states within cell populations. ", 
    "Raw data sourced from Bodenmiller et al. (2012); cell population labels reproduced from ", 
    "Nowicka et al. (2017). ", 
    "See Weber et al. (2019) (paper introducing 'diffcyt' framework), Supplementary Note 1, ", 
    "for more details."
  ), 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataClass = "SummarizedExperiment", 
  RDataPath = "HDCytoData/Weber_BCR_XL_sim/Weber_BCR_XL_sim_less_distinct_less_50pc_SE.rda", 
  stringsAsFactors = FALSE
)

# additional replicates
df_Weber_BCR_XL_sim_less_distinct_less_75pc_SE <- df_Weber_BCR_XL_sim_less_distinct_less_50pc_SE
df_Weber_BCR_XL_sim_less_distinct_less_75pc_SE$Title <- gsub("less_50pc", "less_75pc", df_Weber_BCR_XL_sim_less_distinct_less_75pc_SE$Title)
df_Weber_BCR_XL_sim_less_distinct_less_75pc_SE$RDataPath <- gsub("less_50pc", "less_75pc", df_Weber_BCR_XL_sim_less_distinct_less_75pc_SE$RDataPath)

# flowSets
df_Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet <- df_Weber_BCR_XL_sim_less_distinct_less_50pc_SE
df_Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet$Title)
df_Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet$RDataPath)

df_Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet <- df_Weber_BCR_XL_sim_less_distinct_less_75pc_SE
df_Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet$Title <- gsub("_SE$", "_flowSet", df_Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet$Title)
df_Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet$RDataClass <- "flowSet"
df_Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet$RDataPath <- gsub("_SE.rda", "_flowSet.rda", df_Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet$RDataPath)



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
  df_Bodenmiller_BCR_XL_flowSet, 
  df_Weber_AML_sim_main_5pc_SE, 
  df_Weber_AML_sim_main_5pc_flowSet, 
  df_Weber_AML_sim_main_1pc_SE, 
  df_Weber_AML_sim_main_1pc_flowSet, 
  df_Weber_AML_sim_main_0.1pc_SE, 
  df_Weber_AML_sim_main_0.1pc_flowSet, 
  df_Weber_AML_sim_main_blasts_all_SE, 
  df_Weber_AML_sim_main_blasts_all_flowSet, 
  df_Weber_AML_sim_null_SE, 
  df_Weber_AML_sim_null_flowSet, 
  df_Weber_AML_sim_random_seeds_5pc_SE, 
  df_Weber_AML_sim_random_seeds_5pc_flowSet, 
  df_Weber_AML_sim_random_seeds_1pc_SE, 
  df_Weber_AML_sim_random_seeds_1pc_flowSet, 
  df_Weber_AML_sim_random_seeds_0.1pc_SE, 
  df_Weber_AML_sim_random_seeds_0.1pc_flowSet, 
  df_Weber_AML_sim_less_distinct_5pc_SE, 
  df_Weber_AML_sim_less_distinct_5pc_flowSet, 
  df_Weber_AML_sim_less_distinct_1pc_SE, 
  df_Weber_AML_sim_less_distinct_1pc_flowSet, 
  df_Weber_AML_sim_less_distinct_0.1pc_SE, 
  df_Weber_AML_sim_less_distinct_0.1pc_flowSet, 
  df_Weber_BCR_XL_sim_main_SE, 
  df_Weber_BCR_XL_sim_main_flowSet, 
  df_Weber_BCR_XL_sim_null_rep1_SE, 
  df_Weber_BCR_XL_sim_null_rep1_flowSet, 
  df_Weber_BCR_XL_sim_null_rep2_SE, 
  df_Weber_BCR_XL_sim_null_rep2_flowSet, 
  df_Weber_BCR_XL_sim_null_rep3_SE, 
  df_Weber_BCR_XL_sim_null_rep3_flowSet, 
  df_Weber_BCR_XL_sim_random_seeds_rep1_SE, 
  df_Weber_BCR_XL_sim_random_seeds_rep1_flowSet, 
  df_Weber_BCR_XL_sim_random_seeds_rep2_SE, 
  df_Weber_BCR_XL_sim_random_seeds_rep2_flowSet, 
  df_Weber_BCR_XL_sim_random_seeds_rep3_SE, 
  df_Weber_BCR_XL_sim_random_seeds_rep3_flowSet, 
  df_Weber_BCR_XL_sim_less_distinct_less_50pc_SE, 
  df_Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet, 
  df_Weber_BCR_XL_sim_less_distinct_less_75pc_SE, 
  df_Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet
)

# Save as .csv file
write.csv(df_all, file = "../extdata/metadata.csv", row.names = FALSE)


