\name{Krieg_Anti_PD_1}
\docType{data}
\alias{Krieg_Anti_PD_1}
\alias{Krieg_Anti_PD_1_SE}
\alias{Krieg_Anti_PD_1_flowSet}
\title{
'Krieg_Anti_PD_1' dataset
}
\description{
Mass cytometry (CyTOF) dataset from Krieg et al. (2018), consisting of 20 baseline samples 
(prior to treatment) of peripheral blood from melanoma skin cancer patients subsequently treated 
with anti-PD-1 immunotherapy. The samples are split across 2 conditions (non-responders and 
responders) and 2 batches. This dataset can be used to benchmark algorithms for differential 
analysis, in particular detection of differentially abundant rare cell populations.
}
\details{
This is a mass cytometry (CyTOF) dataset from Krieg et al. (2018), who used mass cytometry to 
characterize immune cell subsets in peripheral blood from melanoma skin cancer patients treated 
with anti-PD-1 immunotherapy. This study found that the frequency of CD14+CD16−HLA-DRhi monocytes 
in baseline samples (taken from patients prior to treatment) was a strong predictor of survival 
in response to immunotherapy treatment. In particular, the frequency of a small subpopulation of 
CD14+CD33+HLA-DRhiICAM-1+CD64+CD141+CD86+CD11c+CD38+PD-L1+CD11b+ monocytes in baseline 
samples was strongly associated with responder status following immunotherapy treatment. 
Note that this dataset contains a strong batch effect, due to sample acquisition on two 
different days (Krieg et al., 2018).

This dataset can be used to benchmark algorithms for differential analysis, in particular 
detection of differentially abundant rare cell populations (i.e. the small subpopulation of 
CD14+CD33+HLA-DRhiICAM-1+CD64+CD141+CD86+CD11c+CD38+PD-L1+CD11b+ monocytes).

The dataset contains 20 baseline samples (i.e. samples taken prior to treatment), from patients 
subsequently classified into 2 groups (9 non-responders and 11 responders). Samples are also 
split across 2 batches ('batch23' and 'batch29'), due to sample acquisition on two different days. 
The total number of cells is 85,715.

There are 24 'cell type' markers used to characterize cell subpopulations. (One additional 
cell type marker -- CD45 -- is also available, but should be excluded from most analyses since 
almost all cells show very high expression of CD45; so it does not help distinguish 
subpopulations, and may dominate other signals. Therefore, CD45 has been classified as 'none' 
in the \code{marker_info} table.)

The dataset is provided in two Bioconductor object formats: \code{\link{SummarizedExperiment}} 
and \code{\link{flowSet}}. In each case, cells are stored in rows, and protein markers in 
columns (this is the usual format used for flow and mass cytometry data).

For the \code{link{SummarizedExperiment}}, row and column metadata can be accessed with the 
\code{\link{rowData}} and \code{\link{colData}} accessor functions from the 
\code{SummarizedExperiment} package. The row data contains group IDs, batch IDs, and 
sample IDs. The column data contains channel names, protein marker names, and a 
factor \code{marker_class} to identify the class of each protein marker ('cell type' or 'cell state'; 
for this dataset, all proteins are cell type markers). The expression values for each cell can be 
accessed with \code{\link{assay}}. The expression values are formatted as a single table.

For the \code{\link{flowSet}}, the expression values are stored in a separate table for each 
sample. Each sample is represented by one \code{\link{flowFrame}} object within the overall 
\code{flowSet}. The expression values can be accessed with the \code{\link{exprs}} function from the 
\code{\link{flowCore}} package. Row metadata is stored as additional columns of data within 
the \code{flowFrame} for each sample; note that factor values are converted to numeric values, 
since the expression tables must be numeric matrices. Channel names are stored in the column names 
of the expression tables. Additional row and column metadata is stored in the \code{description} 
slots, which can be accessed with the \code{\link{description}} accessor function for the 
individual \code{flowFrames}; this includes filenames, additional sample information, and 
additional marker information.

Prior to performing any downstream analyses, the expression values should be transformed. 
A standard transformation used for mass cytometry data is the \code{\link{asinh}} with 
\code{cofactor = 5}.

Original source: Krieg et al. (2018): https://www.ncbi.nlm.nih.gov/pubmed/29309059

Original link to raw data (FlowRepository, FR-FCM-ZY34): http://flowrepository.org/id/FR-FCM-ZY34

This dataset was previously used to benchmark algorithms for differential analysis in our article, 
Weber et al. (2018): https://www.biorxiv.org/content/early/2018/11/22/349738. (For additional 
details on the dataset, see Supplementary Note 3: Benchmark datasets.)

Data files are also available from FlowRepository (FR-FCM-ZYL8): http://flowrepository.org/id/FR-FCM-ZYL8
}
\usage{
Krieg_Anti_PD_1_SE(metadata = FALSE)
Krieg_Anti_PD_1_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Krieg_Anti_PD_1_SE()
Krieg_Anti_PD_1_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Krieg et al. (2018), "High-dimensional single-cell analysis predicts response to anti-PD-1 immunotherapy." 
Nature Medicine, 24, 144-153: https://www.ncbi.nlm.nih.gov/pubmed/29309059

Weber et al. (2018). "diffcyt: Differential discovery in high-dimensional cytometry via high-resolution 
clustering." bioRxiv: https://www.biorxiv.org/content/early/2018/11/22/349738
}
\keyword{datasets}
