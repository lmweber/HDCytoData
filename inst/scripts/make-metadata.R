df <- data.frame(Title = "bodenmiller2012bcrxl",
                 Description = "Abundance measurements from Bodenmiller et al (2012)",
                 BiocVersion = "3.7",
                 Genome = NA,
                 SourceType = "FCS",
                 SourceUrl = "",
                 SourceVersion = "",
                 Species = "Homo sapiens",
                 TaxonomyId = "9606",
                 Coordinate_1_based = NA,
                 DataProvider = "",
                 Maintainer = "",
                 RDataClass = "SummarizedExperiment",
                 DispatchClass = "Rda",
                 RDataPath = "",
                 ResourceName = "",
                 stringsAsFactors = FALSE)


write.csv(df, file = "../extdata/metadata.csv", row.names = FALSE)

