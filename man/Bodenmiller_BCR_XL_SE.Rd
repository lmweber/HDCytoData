\name{Bodenmiller_BCR_XL}
\docType{data}
\alias{Bodenmiller_BCR_XL}
\alias{Bodenmiller_BCR_XL_SE}
\alias{Bodenmiller_BCR_XL_flowSet}
\title{
'Bodenmiller_BCR_XL' data set
}
\description{
Mass cytometry (CyTOF) data set from Bodenmiller et al. (2012), consisting of paired samples 
of healthy immune cells in unstimulated and BCR-XL stimulated condition.
}
\details{
This is a mass cytometry (CyTOF) data set originally from Bodenmiller et al. (2012). The data 
set consists of paired samples of healthy peripheral blood mononuclear cells (PBMCs), where 
one sample from each pair was stimulated with B cell receptor / Fc receptor cross-linker (BCR-XL). 
This creates a strong differential expression signal for several signaling markers in several 
cell populations, especially B cells. The strongest observed differential signal is for the 
signaling marker phosphorylated S6 (pS6) in B cells (see Nowicka et al., 2017, Figure 29).

The data set contains 16 samples (8 paired samples); a total of 172,791 cells; and a total of 
24 protein markers. The markers consist of 10 'cell type' markers (which can be used to define 
cell populations or clusters), and 14 'cell state' or signaling markers.

Reference cell population or cluster labels are available from Nowicka et al. (2017), where 
these were generated using a strategy of expert-guided manual merging of automatically generated 
clusters from the FlowSOM clustering algorithm (Van Gassen et al., 2015).

The data set is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column meta-data can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains group IDs, patient IDs, 
sample IDs, and reference population or cluster IDs. The column data contains channel 
names, protein marker names, and a factor \code{marker_class} to identify the class of each 
protein marker ('cell type' or 'cell state'). In this data set, the cell type markers are 
used to define cell populations, and the cell state markers define signaling states. The 
expression values for each cell can be accessed with \code{\link{assay}}. The expression values 
are formatted as a single table (i.e. the 16 samples are stacked into a single matrix).

For the \code{\link{flowSet}}, the row meta-data is stored as additional columns of data 
within the \code{\link{flowFrame}} object for each sample. The row meta-data includes group IDs, 
patient IDs, sample IDs, and reference population or cluster IDs. Note that the factor values 
are converted to numeric values, since the tables of expression values must be numeric matrices. 
The column meta-data consists of protein marker names only, which are stored in the column names 
of the \code{flowFrame} object for each sample (since column names cannot contain multiple 
entries, the \code{marker_class} information to identify cell type and cell state markers 
cannot be included). The expression values for each cell are stored in the remaining columns 
of the \code{flowFrame} object for each sample. The tables of expression values can be accessed 
with the \code{\link{exprs}} function from the \code{\link{flowCore}} package. The expression 
values are stored in one table per sample (i.e. one \code{flowFrame} object per sample, within 
the overall \code{flowSet} object).

Prior to performing any analysis, the expression values should be transformed. A standard 
transformation used for mass cytometry data is the \code{arcsinh} with \code{cofactor = 5}. 
See the vignette for an example.

Raw data files for the original data set are available from Cytobank, experiment 15713 
(https://community.cytobank.org/cytobank/experiments/15713/download_files; additional information 
is also available from the \code{Citrus} wiki page at: https://github.com/nolanlab/citrus/wiki/PBMC-Example-1).
}
\usage{
Bodenmiller_BCR_XL_SE(metadata = FALSE)
Bodenmiller_BCR_XL_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall data set) should be returned only, or if the whole data set should be loaded. Default = FALSE, 
which loads the whole data set.}
}
\examples{
Bodenmiller_BCR_XL_SE()
Bodenmiller_BCR_XL_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Bodenmiller, B., Zunder, E. R., Finck, R., Chen, T. J., Savig, E. S., Bruggner, R. V., 
Simonds, E. F., Bendall, S. C., Sachs, K., Krutzik, P. O., and Nolan, G. P. (2012). \emph{Multiplexed 
mass cytometry profiling of cellular states perturbed by small-molecule regulators.} Nature 
Biotechnology, 30(9):858-867.

Nowicka, M., Krieg, C., Weber, L. M., Hartmann, F. J., Guglietta, S., Becher, B., 
Levesque, M. P., and Robinson, M. D. (2017). \emph{CyTOF workflow: differential discovery in 
high-throughput high-dimensional cytometry datasets.} F1000Research.

Van Gassen, S., Callebaut, B., Van Helden, M. J., Lambrecht, B. N., Demeester, P., Dhaene, T., 
and Saeys, Y. (2015). \emph{FlowSOM: Using Self-Organizing Maps for Visualization and Interpretation 
of Cytometry Data.} Cytometry Part A, 87A:636-645.
}
\keyword{datasets}
