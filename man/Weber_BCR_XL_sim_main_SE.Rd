\name{Weber_BCR_XL_sim}
\docType{data}
\alias{Weber_BCR_XL_sim}
\alias{Weber_BCR_XL_sim_main}
\alias{Weber_BCR_XL_sim_main_SE}
\alias{Weber_BCR_XL_sim_main_flowSet}
\alias{Weber_BCR_XL_sim_null}
\alias{Weber_BCR_XL_sim_null_rep1}
\alias{Weber_BCR_XL_sim_null_rep1_SE}
\alias{Weber_BCR_XL_sim_null_rep1_flowSet}
\alias{Weber_BCR_XL_sim_null_rep2}
\alias{Weber_BCR_XL_sim_null_rep2_SE}
\alias{Weber_BCR_XL_sim_null_rep2_flowSet}
\alias{Weber_BCR_XL_sim_null_rep3}
\alias{Weber_BCR_XL_sim_null_rep3_SE}
\alias{Weber_BCR_XL_sim_null_rep3_flowSet}
\alias{Weber_BCR_XL_sim_random_seeds}
\alias{Weber_BCR_XL_sim_random_seeds_rep1}
\alias{Weber_BCR_XL_sim_random_seeds_rep1_SE}
\alias{Weber_BCR_XL_sim_random_seeds_rep1_flowSet}
\alias{Weber_BCR_XL_sim_random_seeds_rep2}
\alias{Weber_BCR_XL_sim_random_seeds_rep2_SE}
\alias{Weber_BCR_XL_sim_random_seeds_rep2_flowSet}
\alias{Weber_BCR_XL_sim_random_seeds_rep3}
\alias{Weber_BCR_XL_sim_random_seeds_rep3_SE}
\alias{Weber_BCR_XL_sim_random_seeds_rep3_flowSet}
\alias{Weber_BCR_XL_sim_less_distinct}
\alias{Weber_BCR_XL_sim_less_distinct_less_50pc}
\alias{Weber_BCR_XL_sim_less_distinct_less_50pc_SE}
\alias{Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet}
\alias{Weber_BCR_XL_sim_less_distinct_less_75pc}
\alias{Weber_BCR_XL_sim_less_distinct_less_75pc_SE}
\alias{Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet}
\title{
'Weber_BCR_XL_sim' semi-simulated datasets
}
\description{
Semi-simulated mass cytometry (CyTOF) datasets from Weber et al. (2019), constructed by 
randomly splitting unstimulated (reference) samples of PBMCs (peripheral blood mononuclear cells) 
into two halves, and replacing B cells in one half with stimulated (BCR-XL) B cells from 
corresponding paired samples. These datasets can be used to benchmark differential analysis 
algorithms used to test for differential states within cell populations. Raw data sourced from 
Bodenmiller et al. (2012); cell population labels reproduced from Nowicka et al. (2017). 
See Weber et al. (2019) Supplementary Note 1, for more details.
}
\details{
This is a set of semi-simulated mass cytometry (CyTOF) datasets, generated for benchmarking 
purposes in our paper introducing the 'diffcyt' framework (Weber et al., 2019).

The datasets are constructed by randomly splitting unstimulated (reference) samples of PBMCs 
(peripheral blood mononuclear cells) into two halves, and replacing B cells in one half with 
stimulated (BCR-XL) B cells from corresponding paired samples. Strong differential expression 
signals exist for several signaling state markers in B cells between the stimulated (BCR-XL) 
and unstimulated (reference) conditions; in particular phosphorylated S6 (pS6).

These datasets can be used to benchmark differential analysis algorithms used to test for 
differential states within cell populations.

The raw data consists of 8 paired samples (i.e. 16 samples in total), and a total of 172,791 
cells. The dataset contains expression levels of 24 protein markers (10 surface markers 
used to define cell populations, and 14 intracellular signaling markers). Cell population 
labels are reproduced from Nowicka et al. (2017). For more details, see Weber et al. (2019), 
Supplementary Note 1 (in particular Supplementary Tables 3 and 4).

Multiple simulations are available, as described in our paper (Weber et al., 2019). These are 
stored in the objects listed below.

In each case, the objects are available in both \code{\link{SummarizedExperiment}} and 
\code{\link{flowSet}} formats, with cells stored in rows, and protein markers in columns (i.e. 
the usual format for cytometry data). After loading the datasets, they can be inspected using the 
standard accessor functions for either \code{SummarizedExperiments} or \code{flowSets} (e.g. for 
\code{SummarizedExperiments}: \code{\link{rowData}}, \code{\link{colData}}, \code{\link{assays}}, 
and \code{\link{metadata}}).

For the \code{SummarizedExperiments}: \code{assays} contain tables of expression values 
(with multiple objects for datasets with multiple replicates; note that the replicates cannot be 
combined as multiple \code{assays} within a single object because each replicate has different row data); 
\code{rowData} contains group IDs, patient IDs, sample IDs, cell population IDs, and columns identifying 
B cells and spike-in cells; 
\code{colData} contains channel names, marker names, and marker classes; and 
\code{metadata} contains experiment information and number of cells.

For the \code{flowSets}: individual \code{flowFrames} within the \code{flowSet} contain tables 
of expression values (with multiple \code{flowSet} objects for datasets with multiple replicates); 
row data is stored as additional columns of numeric values within the expression tables; 
column data is stored in the \code{pData(parameters())} slot of the individual \code{flowFrames}; and 
additional information (e.g. experiment information, marker information, replicate information, 
and lookup tables to identify row data values) is stored in the \code{description()} slot of 
the \code{flowFrames}.

\strong{Main simulations}
\itemize{
\item Weber_BCR_XL_sim_main_SE
\item Weber_BCR_XL_sim_main_flowSet
}

\strong{Additional simulations: null simulations}
\emph{
Separate files for each replicate.
}
\itemize{
\item Weber_BCR_XL_sim_null_rep1_SE
\item Weber_BCR_XL_sim_null_rep1_flowSet
\item Weber_BCR_XL_sim_null_rep2_SE
\item Weber_BCR_XL_sim_null_rep2_flowSet
\item Weber_BCR_XL_sim_null_rep3_SE
\item Weber_BCR_XL_sim_null_rep3_flowSet
}

\strong{Additional simulations: modified random seeds}
\emph{
Separate files for each replicate.
}
\itemize{
\item Weber_BCR_XL_sim_random_seeds_rep1_SE
\item Weber_BCR_XL_sim_random_seeds_rep1_flowSet
\item Weber_BCR_XL_sim_random_seeds_rep2_SE
\item Weber_BCR_XL_sim_random_seeds_rep2_flowSet
\item Weber_BCR_XL_sim_random_seeds_rep3_SE
\item Weber_BCR_XL_sim_random_seeds_rep3_flowSet
}

\strong{Additional simulations: 'less distinct' spike-in cells}
\emph{
Separate files for each replicate.
}
\itemize{
\item Weber_BCR_XL_sim_less_distinct_less_50pc_SE
\item Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet
\item Weber_BCR_XL_sim_less_distinct_less_75pc_SE
\item Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet
}

Note that prior to performing any downstream analyses, the expression values should be 
transformed. A standard transformation used for mass cytometry data is the \code{\link{asinh}} 
with \code{cofactor = 5}.

The raw data is sourced from Bodenmiller et al. (2012), and cell population labels are reproduced 
from Nowicka et al. (2017). See Weber et al. (2019), Supplementary Note 1, for more details.

Original link to raw data (Cytobank, experiment 15713): https://community.cytobank.org/cytobank/experiments/15713/download_files

Additional information (\code{Citrus} wiki page): https://github.com/nolanlab/citrus/wiki/PBMC-Example-1

Data files are also available from FlowRepository (FR-FCM-ZYL8): http://flowrepository.org/id/FR-FCM-ZYL8
}
\usage{
Weber_BCR_XL_sim_main_SE(metadata = FALSE)
Weber_BCR_XL_sim_main_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_null_rep1_SE(metadata = FALSE)
Weber_BCR_XL_sim_null_rep1_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_null_rep2_SE(metadata = FALSE)
Weber_BCR_XL_sim_null_rep2_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_null_rep3_SE(metadata = FALSE)
Weber_BCR_XL_sim_null_rep3_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_random_seeds_rep1_SE(metadata = FALSE)
Weber_BCR_XL_sim_random_seeds_rep1_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_random_seeds_rep2_SE(metadata = FALSE)
Weber_BCR_XL_sim_random_seeds_rep2_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_random_seeds_rep3_SE(metadata = FALSE)
Weber_BCR_XL_sim_random_seeds_rep3_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_less_distinct_less_50pc_SE(metadata = FALSE)
Weber_BCR_XL_sim_less_distinct_less_50pc_flowSet(metadata = FALSE)
Weber_BCR_XL_sim_less_distinct_less_75pc_SE(metadata = FALSE)
Weber_BCR_XL_sim_less_distinct_less_75pc_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Weber_BCR_XL_sim_main_SE()
Weber_BCR_XL_sim_main_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Bodenmiller et al. (2012). "Multiplexed mass cytometry profiling of cellular states perturbed 
by small-molecule regulators." Nature Biotechnology, 30(9), 858-867: https://www.ncbi.nlm.nih.gov/pubmed/22902532

Nowicka et al. (2017). "CyTOF workflow: differential discovery in high-throughput high-dimensional 
cytometry datasets." F1000Research, v2: https://f1000research.com/articles/6-748/v2

Weber et al. (2019). "diffcyt: Differential discovery in high-dimensional cytometry via high-resolution 
clustering." Communications Biology, 2:183: https://www.ncbi.nlm.nih.gov/pubmed/31098416
}
\keyword{datasets}
