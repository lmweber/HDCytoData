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

\emph{Overview}

This package contains a set of publicly available high-dimensional flow cytometry and 
mass cytometry (CyTOF) datasets, which have been formatted into \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}} Bioconductor object formats.

The objects contain the cell-level expression values, as well as row and column metadata. 
The row metadata includes sample IDs, group IDs, and true cell population labels or cluster labels 
(where available). The column metadata includes channel names, protein marker names, and 
protein marker classes (cell type, cell state, as well as non protein marker columns).

These datasets have been used in our previous work and publications for benchmarking purposes, 
e.g. to benchmark clustering algorithms or methods for differential analysis. They are provided 
here in the \code{SummarizedExperiment} and \code{flowSet} formats to make them easier to access.

The package contains the following datasets, which can be grouped into datasets useful for 
benchmarking either (i) clustering algorithms or (ii) methods for differential analysis.

Clustering:

\itemize{
\item \code{\link{Levine_32dim}}
\item \code{\link{Levine_13dim}}
\item \code{\link{Samusik_01}}
\item \code{\link{Samusik_all}}
\item \code{\link{Nilsson_rare}}
\item \code{\link{Mosmann_rare}}
}

Differential analysis:

\itemize{
\item \code{\link{Krieg_Anti_PD_1}}
\item \code{\link{Bodenmiller_BCR_XL}}
}


\emph{Programmatic access to list of datasets}

An updated list of all available datasets can also be obtained programmatically using the 
\code{ExperimentHub} accessor functions, as follows. This retrieves a table of metadata from 
the \code{ExperimentHub} database, which includes information such as the ExperimentHub ID, 
title, and description for each dataset.

\code{ehub <- ExperimentHub()  # create ExperimentHub instance} \cr
\code{ehub <- query(ehub, "HDCytoData")  # find HDCytoData datasets} \cr
\code{md <- as.data.frame(mcols(ehub))  # retrieve metadata table}


\emph{Additional details}

For additional details on each dataset, including references and raw data sources, see the 
help files for each dataset.

For a short tutorial showing how to load the data objects, see the "HDCytoData package" vignette.

Note that flow and mass cytometry datasets should be transformed prior to performing any 
downstream analyses, such as clustering. Standard transforms include the \code{\link{asinh}} with 
\code{cofactor} parameter equal to 5 (for mass cytometry data) or 150 (for flow cytometry data).

The steps to prepare each data object from the raw data files are included in the \code{make-data} 
scripts in the directory \code{inst/scripts}.

File sizes are listed in the help files for the datasets. The \code{removeCache} 
function from the \code{ExperimentHub} package can be used to clear the local download cache.
}
\keyword{datasets}
