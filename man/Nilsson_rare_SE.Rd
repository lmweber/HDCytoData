\name{Nilsson_rare}
\docType{data}
\alias{Nilsson_rare}
\alias{Nilsson_rare_SE}
\alias{Nilsson_rare_flowSet}
\title{
'Nilsson_rare' dataset
}
\description{
Flow cytometry dataset from Nilsson et al. (2013), containing 13 dimensions (surface protein 
markers). Manually gated cell population labels are available for one rare population of 
hematopoietic stem cells (HSCs). Cells are human bone marrow cells from a single healthy donor. 
This dataset can be used to benchmark clustering algorithms for rare cell populations.
}
\details{
This is a 13-dimensional flow cytometry dataset, consisting of expression levels of 
13 surface protein markers. Cell population labels are available for one rare population of 
hematopoietic stem cells (HSCs). Cells are human bone marrow cells from a single healthy donor.

This dataset can be used to benchmark clustering algorithms for rare cell populations.

The dataset contains cells from a single patient; 
a total of 44,140 cells (including 358 manually gated cells from the rare population of interest); 
and a total of 13 surface marker proteins.

The dataset is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column metadata can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains the manually gated 
cell population IDs. The column data contains channel names, protein marker names, and a 
factor \code{marker_class} to identify the class of each protein marker ('cell type' or 'cell state'; 
for this dataset, all proteins are cell type markers). The expression values for each cell can be 
accessed with \code{\link{assay}}. The expression values are formatted as a single table.

For the \code{\link{flowSet}}, the expression values are stored in a separate table for each 
sample. Each sample is represented by one \code{\link{flowFrame}} object within the overall 
\code{flowSet} (note that for this dataset, there is only one sample). 
The expression values can be accessed with the \code{\link{exprs}} function from the 
\code{\link{flowCore}} package. Row metadata is stored as additional columns of data within 
the \code{flowFrame} for each sample; note that factor values are converted to numeric values, 
since the expression tables must be numeric matrices. Channel names are stored in the column names 
of the expression tables. Additional row and column metadata is stored in the \code{description} 
slots, which can be accessed with the \code{\link{description}} accessor function for the 
individual \code{flowFrames}; this includes additional sample information (where available), 
marker information, and cell population information.

Prior to performing any downstream analyses, the expression values should be transformed. 
A standard transformation used for flow cytometry data is the \code{\link{asinh}} with 
\code{cofactor = 150}.

File sizes: 2.4 MB (\code{SummarizedExperiment} and \code{flowSet}).

Original source: Figure 2 in Nilsson et al. (2013): https://www.ncbi.nlm.nih.gov/pubmed/23839904

Original link to raw data: http://flowrepository.org/id/FR-FCM-ZZ6L

This dataset was previously used to benchmark clustering algorithms for high-dimensional 
cytometry in our article, Weber and Robinson (2016): https://www.ncbi.nlm.nih.gov/pubmed/27992111

Data files are also available from FlowRepository (FR-FCM-ZZPH): http://flowrepository.org/id/FR-FCM-ZZPH
}
\usage{
Nilsson_rare_SE(metadata = FALSE)
Nilsson_rare_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Nilsson_rare_SE()
Nilsson_rare_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Nilsson et al. (2013), "Frequency determination of rare populations by flow cytometry: A 
hematopoietic stem cell perspective", Cytometry Part A, 83A, 721-727: http://www.ncbi.nlm.nih.gov/pubmed/23839904

Weber and Robinson (2016), "Comparison of clustering methods for high-dimensional single-cell 
flow and mass cytometry data", Cytometry Part A, 89A, 1084-1096: https://www.ncbi.nlm.nih.gov/pubmed/27992111
}
\keyword{datasets}
