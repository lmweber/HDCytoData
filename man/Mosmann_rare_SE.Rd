\name{Mosmann_rare}
\docType{data}
\alias{Mosmann_rare}
\alias{Mosmann_rare_SE}
\alias{Mosmann_rare_flowSet}
\title{
'Mosmann_rare' dataset
}
\description{
Flow cytometry dataset from Mosmann et al. (2014), containing 14 dimensions (7 surface protein 
markers and 7 signaling markers). Manually gated cell population labels are available for one 
rare population of activated (cytokine-producing) memory CD4 T cells. Cells are human peripheral 
blood cells exposed to influenza antigens, from a single healthy donor. This dataset can be 
used to benchmark clustering algorithms for rare cell populations.
}
\details{
This is a 14-dimensional flow cytometry dataset, consisting of expression levels of 7 surface 
protein markers and 7 signaling markers. Cell population labels are available for one rare 
population of activated (cytokine-producing) memory CD4 T cells. Cells are human peripheral 
blood cells exposed to influenza antigens, from a single healthy donor.

This dataset can be used to benchmark clustering algorithms for rare cell populations.

The dataset contains cells from a single patient; 
a total of 396,460 cells (including 109 manually gated cells from the rare population of interest); 
and a total of 14 protein markers (7 surface protein markers and 7 signaling markers).

The dataset is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column metadata can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains the manually gated 
cell population IDs. The column data contains channel names, protein marker names, and a 
factor \code{marker_class} to identify the class of each protein marker ('cell type' or 'cell state'). 
The expression values for each cell can be accessed with \code{\link{assay}}. The expression 
values are formatted as a single table.

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

Original source: Figure 4 in Mosmann et al. (2014): https://www.ncbi.nlm.nih.gov/pubmed/24532172

Original link to raw data: http://flowrepository.org/id/FR-FCM-ZZ8J (filename: 
"JMW034-J16OFVQX_G2 0o1 3_D07.fcs"; see Supplementary Information file 3 for full list of filenames)

This dataset was previously used to benchmark clustering algorithms for high-dimensional 
cytometry in our article, Weber and Robinson (2016): https://www.ncbi.nlm.nih.gov/pubmed/27992111

Data files are also available from FlowRepository (FR-FCM-ZZPH): http://flowrepository.org/id/FR-FCM-ZZPH
}
\usage{
Mosmann_rare_SE()
Mosmann_rare_flowSet()
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Mosmann_rare_SE()
Mosmann_rare_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Mosmann et al. (2014), "SWIFT - Scalable clustering for automated identification of rare 
cell populations in large, high-dimensional flow cytometry datasets, Part 2: Biological 
evaluation", Cytometry Part A, 85A, 422-433: https://www.ncbi.nlm.nih.gov/pubmed/24532172

Weber and Robinson (2016), "Comparison of clustering methods for high-dimensional single-cell 
flow and mass cytometry data", Cytometry Part A, 89A, 1084-1096: https://www.ncbi.nlm.nih.gov/pubmed/27992111
}
\keyword{datasets}
