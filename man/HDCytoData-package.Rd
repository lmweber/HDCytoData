\name{HDCytoData}
\docType{data}
\alias{HDCytoData}
\alias{HDCytoData-package}
\title{
Data package of high-dimensional cytometry datasets
}
\description{
Data package containing a collection of high-dimensional cytometry datasets saved in 
\code{\link{SummarizedExperiment}} and \code{\link{flowSet}} Bioconductor object formats, 
hosted on Bioconductor ExperimentHub.
}
\details{
This package contains a set of publicly available high-dimensional flow cytometry and 
mass cytometry (CyTOF) datasets, which have been formatted into \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}} Bioconductor object formats.

The objects contain the cell-level expression values, as well as row and column metadata. 
The row metadata includes sample IDs, group IDs, and true cell population labels or cluster labels 
(where available). The column metadata includes channel names, protein marker names, and 
protein marker classes (cell type or cell state).

These datasets have been used in our previous work and publications for benchmarking purposes, 
e.g. to benchmark clustering algorithms or methods for differential analysis. They are provided 
here in the \code{SummarizedExperiment} and \code{flowSet} formats to make them easier to access.

The package contains the following datasets:
\itemize{
\item \code{\link{Levine_32dim}}
\item \code{\link{Levine_13dim}}
\item \code{\link{Samusik_01}}
\item \code{\link{Samusik_all}}
\item \code{\link{Bodenmiller_BCR_XL}}
}

For additional details on each dataset, including references and raw data sources, see the 
help files for each dataset.

For a short tutorial showing how to load the data objects, see the package vignette.

The steps to prepare each data object from the raw data files are included in the \code{make-data} 
scripts in the directory \code{inst/scripts}.
}
\keyword{datasets}
