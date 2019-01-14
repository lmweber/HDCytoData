\name{Bodenmiller_BCR_XL}
\docType{data}
\alias{Bodenmiller_BCR_XL}
\alias{Bodenmiller_BCR_XL_SE}
\alias{Bodenmiller_BCR_XL_flowSet}
\title{
'Bodenmiller_BCR_XL' dataset
}
\description{
Mass cytometry (CyTOF) dataset from Bodenmiller et al. (2012), consisting of paired samples 
of peripheral blood cells from healthy individuals, where one sample from each pair was 
stimulated with B cell receptor / Fc receptor cross-linker (BCR-XL). This dataset can be 
used to benchmark algorithms for differential analysis, in particular discovery of differential 
states within cell populations. The dataset contains strong differential expression signals 
for several signaling markers in several cell populations; one of the strongest effects is 
differential expression of phosphorylated S6 (pS6) in B cells.
}
\details{
This is a mass cytometry (CyTOF) dataset from Bodenmiller et al. (2012), consisting of paired 
samples of peripheral blood cells from healthy individuals, where one sample from each pair was 
stimulated with B cell receptor / Fc receptor cross-linker (BCR-XL).

This dataset can be used to benchmark algorithms for differential analysis, in particular 
discovery of differential states within cell populations. The dataset contains strong 
differential expression signals for several signaling markers in several cell populations; 
one of the strongest effects is differential expression of phosphorylated S6 (pS6) in B cells.

There are 8 paired samples (i.e. 16 samples in total), and a total of 172,791 cells. The 
dataset contains expression levels of 24 protein markers (10 surface markers or 'cell type' 
markers used to define cell populations or clusters, and 14 intracellular signaling markers or 
'cell state' markers).

Cell population or cluster labels are available from Nowicka et al. (2017), where 
these were generated using a strategy of expert-guided manual merging of automatically generated 
clusters from the FlowSOM clustering algorithm (Van Gassen et al., 2015).

The dataset is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column metadata can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains group IDs, patient IDs, sample IDs, 
and cell population IDs. The column data contains channel names, protein marker names, and a 
factor \code{marker_class} to identify the class of each protein marker ('cell type' or 'cell state'). 
The expression values for each cell can be accessed with \code{\link{assay}}. The expression 
values are formatted as a single table.

For the \code{\link{flowSet}}, the expression values are stored in a separate table for each 
sample. Each sample is represented by one \code{\link{flowFrame}} object within the overall 
\code{flowSet}. The expression values can be accessed with the \code{\link{exprs}} function from the 
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

Original source: Bodenmiller et al. (2012): https://www.ncbi.nlm.nih.gov/pubmed/22902532

Original link to raw data (Cytobank, experiment 15713): https://community.cytobank.org/cytobank/experiments/15713/download_files

Additional information (\code{Citrus} wiki page): https://github.com/nolanlab/citrus/wiki/PBMC-Example-1

Cell population labels from: Nowicka et al. (2017), v2: https://f1000research.com/articles/6-748/v2

This dataset has previously been used to benchmark algorithms for differential analysis by 
ourselves and other authors, including 
Bruggner et al. (2014) (https://www.ncbi.nlm.nih.gov/pubmed/24979804/), 
Nowicka et al. (2017) (https://f1000research.com/articles/6-748/v2), and 
Weber et al. (2018) (https://www.biorxiv.org/content/early/2018/11/22/349738).

Data files are also available from FlowRepository (FR-FCM-ZYL8): http://flowrepository.org/id/FR-FCM-ZYL8
}
\usage{
Bodenmiller_BCR_XL_SE(metadata = FALSE)
Bodenmiller_BCR_XL_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Bodenmiller_BCR_XL_SE()
Bodenmiller_BCR_XL_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Bodenmiller et al. (2012). "Multiplexed mass cytometry profiling of cellular states perturbed 
by small-molecule regulators." Nature Biotechnology, 30(9), 858-867: https://www.ncbi.nlm.nih.gov/pubmed/22902532

Bruggner et al. (2014), "Automated identification of stratifying signatures in cellular subpopulations." 
PNAS, 111(26), E2770-E2777: https://www.ncbi.nlm.nih.gov/pubmed/24979804/

Nowicka et al. (2017). "CyTOF workflow: differential discovery in high-throughput high-dimensional 
cytometry datasets." F1000Research, v2: https://f1000research.com/articles/6-748/v2

Van Gassen et al. (2015). "FlowSOM: Using self-organizing maps for visualization and interpretation 
of cytometry data." Cytometry Part A, 87A, 636-645: https://www.ncbi.nlm.nih.gov/pubmed/25573116

Weber et al. (2018). "diffcyt: Differential discovery in high-dimensional cytometry via high-resolution 
clustering." bioRxiv: https://www.biorxiv.org/content/early/2018/11/22/349738
}
\keyword{datasets}
