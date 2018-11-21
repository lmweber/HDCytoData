# Metadata for 'SummarizedExperiment' format
df_Bodenmiller_BCR_XL_SE <- data.frame(
  Title = "Bodenmiller_BCR_XL_SE", 
  Description = paste0(
    "Mass cytometry data from Bodenmiller et al. (2012). Paired samples of healthy PBMCs, ", 
    "where one sample from each pair was stimulated with B cell receptor / Fc receptor ", 
    "cross-linker (BCR-XL)."), 
  BiocVersion = "3.8", 
  Genome = NA, 
  SourceType = "FCS", 
  SourceUrl = "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/", 
  SourceVersion = NA, 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  Coordinate_1_based = TRUE, 
  DataProvider = "", 
  Maintainer = "Lukas M. Weber <lukmweber@gmail.com>", 
  RDataClass = "SummarizedExperiment", 
  DispatchClass = "Rda", 
  RDataPath = "HDCytoData/Bodenmiller_BCR_XL_SE.rda", 
  stringsAsFactors = FALSE
)

# Metadata for 'flowSet' format
df_Bodenmiller_BCR_XL_flowSet <- df_Bodenmiller_BCR_XL_SE
df_Bodenmiller_BCR_XL_flowSet$Title <- "Bodenmiller_BCR_XL_flowSet"
df_Bodenmiller_BCR_XL_flowSet$RDataClass <- "flowSet"
df_Bodenmiller_BCR_XL_flowSet$RDataPath <- "HDCytoData/Bodenmiller_BCR_XL_flowSet.rda"

# Combined metadata for both formats
df_all <- rbind(
  df_Bodenmiller_BCR_XL_SE, 
  df_Bodenmiller_BCR_XL_flowSet
)

# Save .csv file
write.csv(df_all, file = "../extdata/metadata.csv", row.names = FALSE)

