Bodenmiller_BCR_XL_SE_desc <- paste0(
  "Mass cytometry data from Bodenmiller et al. (2012). Paired samples of healthy PBMCs, where ", 
  "one sample from each pair was stimulated with B cell receptor / Fc receptor cross-linker ", 
  "(BCR-XL)."
)

df_Bodenmiller_BCR_XL_SE <- data.frame(
  Title = "Bodenmiller_BCR_XL_SE", 
  Description = Bodenmiller_BCR_XL_SE_desc, 
  BiocVersion = "3.7", 
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
  RDataPath = "HDCytoData", 
  stringsAsFactors = FALSE
)

df_Bodenmiller_BCR_XL_flowSet <- df_Bodenmiller_BCR_XL_SE
df_Bodenmiller_BCR_XL_flowSet$Title <- "Bodenmiller_BCR_XL_flowSet"
df_Bodenmiller_BCR_XL_flowSet$RDataClass <- "flowSet"

df_all <- rbind(
  df_Bodenmiller_BCR_XL_SE, 
  df_Bodenmiller_BCR_XL_flowSet
)

write.csv(df_all, file = "../extdata/metadata.csv", row.names = FALSE)

