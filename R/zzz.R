.onLoad <- function(libname, pkgname) {
  fl <- system.file("extdata", "metadata.csv", package = "HDCytoData")
  titles <- read.csv(fl, stringsAsFactors = FALSE)$Title
  createHubAccessors(pkgname, titles)
}