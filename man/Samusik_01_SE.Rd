\name{Samusik_01}
\docType{data}
\alias{Samusik_01}
\alias{Samusik_01_SE}
\alias{Samusik_01_flowSet}
\title{
'Samusik_01' dataset
}
\description{
Mass cytometry (CyTOF) dataset from Samusik et al. (2016), containing 39 dimensions (surface 
protein markers). Manually gated cell population labels are available for 24 populations. 
The full dataset ('Samusik_all') contains cells from 10 replicate bone marrow samples from 
C57BL/6J mice (i.e. samples from 10 different mice); this dataset ('Samusik_01') contains 
the data from sample 01 only. This dataset can be used to benchmark clustering algorithms.
}
\details{
This is a 39-dimensional mass cytometry (CyTOF) data set, consisting of expression levels of 
39 surface marker proteins. Cell population labels are available for 24 manually gated populations. 
Manually gated cell population labels were provided by the original authors. 
The full dataset ('Samusik_all') contains cells from 10 replicate bone marrow samples from 
C57BL/6J mice (i.e. samples from 10 different mice); this dataset ('Samusik_01') contains the data 
from sample 01 only.

This dataset can be used to benchmark clustering algorithms.

The 'Samusik_01' dataset contains cells from 1 mouse (sample '01'); 
a total of 86,864 cells (53,173 manually gated and 33,691 unclassified); 
24 manually gated cell population IDs (as well as 'unassigned'); 
and a total of 39 surface marker proteins.

The dataset is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column metadata can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains sample IDs and manually gated 
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
A standard transformation used for mass cytometry data is the \code{\link{asinh}} with 
\code{cofactor = 5}.

File sizes: 19.9 MB (\code{SummarizedExperiment}), 20.0 MB (\code{flowSet}).

Original source: Samusik et al. (2016): https://www.ncbi.nlm.nih.gov/pubmed/27183440

Original link to raw data (.zip file): "https://web.stanford.edu/~samusik/Panorama BM 1-10.zip"

This dataset was previously used to benchmark clustering algorithms for high-dimensional 
cytometry in our article, Weber and Robinson (2016): https://www.ncbi.nlm.nih.gov/pubmed/27992111

Data files are also available from FlowRepository (FR-FCM-ZZPH): http://flowrepository.org/id/FR-FCM-ZZPH
}
\usage{
Samusik_01_SE(metadata = FALSE)
Samusik_01_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Samusik_01_SE()
Samusik_01_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Samusik et al. (2016), "Automated mapping of phenotype space with single-cell data", Nature Methods, 
13(6), 493-496: https://www.ncbi.nlm.nih.gov/pubmed/27183440

Weber and Robinson (2016), "Comparison of clustering methods for high-dimensional single-cell 
flow and mass cytometry data", Cytometry Part A, 89A, 1084-1096: https://www.ncbi.nlm.nih.gov/pubmed/27992111
}
\keyword{datasets}
