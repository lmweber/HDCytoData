.onLoad <- function(libname, pkgname) {
  fl <- system.file("extdata", "metadata.csv", package = "HDCytoData")
  titles <- read.csv(fl, stringsAsFactors = FALSE)$Title
  # note: temporarily disabled for Bioconductor version 3.8 (due to re-numbering of previous ExperimentHub IDs)
  #createHubAccessors(pkgname, titles)
}

# temporary functions to load Bodenmiller_BCR_XL datasets in Bioconductor version 3.8
# (required due to re-numbering of previous ExperimentHub IDs 1119 and 1120)
# (removed from Bioconductor version 3.9 onwards)
Bodenmiller_BCR_XL_SE <- function(metadata = FALSE) {
  load(url("http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/ExperimentHub/Bodenmiller_BCR_XL_SE.rda"))
  assign("d_SE", d_SE, envir = .GlobalEnv)
  d_SE
}

Bodenmiller_BCR_XL_flowSet <- function(metadata = FALSE) {
  load(url("http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/ExperimentHub/Bodenmiller_BCR_XL_flowSet.rda"))
  assign("d_flowSet", d_flowSet, envir = .GlobalEnv)
  d_flowSet
}

Levine_32dim_SE <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Levine_32dim_flowSet <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Levine_13dim_SE <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Levine_13dim_flowSet <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Samusik_01_SE <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Samusik_01_flowSet <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Samusik_all_SE <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Samusik_all_flowSet <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Nilsson_rare_SE <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Nilsson_rare_flowSet <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Mosmann_rare_SE <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Mosmann_rare_flowSet <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Krieg_Anti_PD_1_SE <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}
Krieg_Anti_PD_1_flowSet <- function(metadata = FALSE) {
  message("Please use Bioconductor version 3.9 or later to access this dataset.")
}

