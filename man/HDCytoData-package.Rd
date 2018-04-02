\name{HDCytoData}
\docType{data}
\alias{HDCytoData}
\alias{HDCytoData-package}
\title{
Data package of high-dimensional cytometry data sets
}
\description{
Data package containing high-dimensional cytometry data sets saved in 
\code{\link[SummarizedExperiment]{SummarizedExperiment}} and \code{\link[flowCore]{flowSet}} 
Bioconductor object formats, hosted on Bioconductor ExperimentHub.
}
\details{
This package contains a set of publicly available high-dimensional flow cytometry and mass 
cytometry (CyTOF) data sets, which have been formatted into the 
\code{\link[SummarizedExperiment]{SummarizedExperiment}} and 
\code{\link[flowCore]{flowSet}} Bioconductor object formats.

The objects contain the cell-level expression values, as well as row and column meta-data, 
including sample IDs, group IDs, true cell population or cluster labels (where available), 
and protein marker names.

These data sets have been used for benchmarking purposes in our previous work and publications, 
e.g. to evaluate the performance of clustering algorithms. They are provided here in the 
\code{SummarizedExperiment} and \code{flowSet} formats to make them easier to access.

Currently, the package contains the following data sets:
\itemize{
\item \code{\link{Bodenmiller_BCR_XL}}
}

For additional details on each data set, including references and raw data sources, see the 
help files for the data sets.

For a short example workflow demonstrating how to load the data objects and use them in an 
analysis workflow, see the package vignette.

The steps to prepare each data object from the raw data files are included in the \code{make-data} 
scripts in the directory \code{inst/scripts}.
}
\keyword{datasets}
