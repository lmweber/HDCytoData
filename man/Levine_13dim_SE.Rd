\name{Levine_13dim}
\docType{data}
\alias{Levine_13dim}
\alias{Levine_13dim_SE}
\alias{Levine_13dim_flowSet}
\title{
'Levine_13dim' dataset
}
\description{
13-dimensional mass cytometry (CyTOF) dataset from Levine et al. (2015), containing 24 manually 
gated immune cell populations from healthy human bone marrow mononuclear cells (BMMCs). This 
dataset can be used to benchmark clustering algorithms.
}
\details{
This is a 13-dimensional mass cytometry (CyTOF) data set, consisting of expression levels of 
13 surface marker proteins. Cell population labels are available for 24 manually gated populations. 
Cells are healthy human bone marrow mononuclear cells (BMMCs), from 1 patient. Manually gated 
cell populations were provided by the original authors.

This dataset can be used to benchmark clustering algorithms.

The dataset contains cells from a single patient; a total of 167,044 cells (81,747 manually gated 
and 85,297 unclassified); and a total of 13 surface marker proteins.

The dataset is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column metadata can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains the manually gated 
cell population labels. The column data contains channel names, protein marker names, and a 
factor \code{marker_class} to identify the class of each protein marker ('cell type' or 'cell state'). 
The expression values for each cell can be accessed with \code{\link{assay}}. The expression values 
are formatted as a single table.

For the \code{\link{flowSet}}, the row metadata is stored as an additional column of data 
within the \code{\link{flowFrame}} object. The row metadata includes manually 
gated cell population labels only. Note that the factor values are converted to numeric values, 
since the tables of expression values must be numeric matrices. The column metadata consists of 
channel names only, which are stored in the column names of the \code{flowFrame} object for 
each sample (since column names cannot contain multiple entries, the \code{marker_class} 
information to identify cell type and cell state markers cannot be included). The expression 
values for each cell are stored in the remaining columns of the \code{flowFrame} object 
(within the overall \code{flowSet} object). The table of expression values can be accessed with 
the \code{\link{exprs}} function from the \code{\link{flowCore}} package. A table to match 
cell population indices to cell population names can be accessed with the \code{\link{description}} 
accessor function for the \code{flowFrame}.

Prior to performing any analysis, the expression values should be transformed. A standard 
transformation used for mass cytometry data is the \code{arcsinh} with \code{cofactor = 5}.

Original source: "benchmark data set 1" in Levine et al. (2015): https://www.ncbi.nlm.nih.gov/pubmed/26095251

Original link to raw data: https://www.cytobank.org/cytobank/experiments/46259 (download the 
FCS files with Actions -> Export -> Download Files -> All FCS Files)

This dataset was previously used to benchmark clustering algorithms for high-dimensional 
cytometry in our article: Weber and Robinson (2016).
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
Levine et al. (2015), "Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells 
that Correlate with Prognosis", Cell, 162, 184-197.

Weber and Robinson (2016), "Comparison of Clustering Methods for High-Dimensional Single-Cell 
Flow and Mass Cytometry Data", Cytometry Part A, 89A, 1084-1096.
}
\keyword{datasets}
