.onLoad <- function(libname, pkgname) {
  fl <- system.file("extdata", "metadata.csv", package="cytofData")
  titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
  createHubAccessors(pkgname, titles)
}