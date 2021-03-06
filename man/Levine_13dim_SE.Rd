\name{Levine_13dim}
\docType{data}
\alias{Levine_13dim}
\alias{Levine_13dim_SE}
\alias{Levine_13dim_flowSet}
\title{
'Levine_13dim' dataset
}
\description{
Mass cytometry (CyTOF) dataset from Levine et al. (2015), containing 13 dimensions (surface 
protein markers). Manually gated cell population labels are available for 24 populations. 
Cells are human bone marrow cells from a single healthy donor. This dataset can be used to benchmark 
clustering algorithms.
}
\details{
This is a 13-dimensional mass cytometry (CyTOF) data set, consisting of expression levels of 
13 surface marker proteins. Cell population labels are available for 24 manually gated populations. 
Cells are human bone marrow cells from a single healthy donor. Manually gated cell population labels 
were provided by the original authors.

This dataset can be used to benchmark clustering algorithms.

The dataset contains cells from a single patient; 
a total of 167,044 cells (81,747 manually gated and 85,297 unclassified); 
24 manually gated cell population IDs (as well as 'unassigned'); 
and a total of 13 surface marker proteins.

The dataset is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column metadata can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains the manually gated 
cell population IDs. The column data contains channel names, protein marker names, and a 
factor \code{marker_class} to identify the class of each protein marker ('cell type', 'cell state', 
as well as 'none' for any non protein marker columns that are not needed for downstream analyses; 
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
A standard transformation used for mass cytometry data is the \code{\link{asinh}} with 
\code{cofactor = 5}.

File sizes: 10.0 MB (\code{SummarizedExperiment} and \code{flowSet}).

Original source: "benchmark data set 1" in Levine et al. (2015): https://www.ncbi.nlm.nih.gov/pubmed/26095251

Original link to raw data: https://www.cytobank.org/cytobank/experiments/46259 (download the 
FCS files with Actions -> Export -> Download Files -> All FCS Files)

This dataset was previously used to benchmark clustering algorithms for high-dimensional 
cytometry in our article, Weber and Robinson (2016): https://www.ncbi.nlm.nih.gov/pubmed/27992111

Data files are also available from FlowRepository (FR-FCM-ZZPH): http://flowrepository.org/id/FR-FCM-ZZPH
}
\usage{
Levine_13dim_SE(metadata = FALSE)
Levine_13dim_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Levine_13dim_SE()
Levine_13dim_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Levine et al. (2015), "Data-driven phenotypic dissection of AML reveals progenitor-like cells 
that correlate with prognosis", Cell, 162, 184-197: https://www.ncbi.nlm.nih.gov/pubmed/26095251

Weber and Robinson (2016), "Comparison of clustering methods for high-dimensional single-cell 
flow and mass cytometry data", Cytometry Part A, 89A, 1084-1096: https://www.ncbi.nlm.nih.gov/pubmed/27992111
}
\keyword{datasets}
