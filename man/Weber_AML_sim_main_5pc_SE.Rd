\name{Weber_AML_sim}
\docType{data}
\alias{Weber_AML_sim}
\alias{Weber_AML_sim_main}
\alias{Weber_AML_sim_main_5pc}
\alias{Weber_AML_sim_main_5pc_SE}
\alias{Weber_AML_sim_main_5pc_flowSet}
\alias{Weber_AML_sim_main_1pc}
\alias{Weber_AML_sim_main_1pc_SE}
\alias{Weber_AML_sim_main_1pc_flowSet}
\alias{Weber_AML_sim_main_0.1pc}
\alias{Weber_AML_sim_main_0.1pc_SE}
\alias{Weber_AML_sim_main_0.1pc_flowSet}
\alias{Weber_AML_sim_main_blasts_all}
\alias{Weber_AML_sim_main_blasts_all_SE}
\alias{Weber_AML_sim_main_blasts_all_flowSet}
\alias{Weber_AML_sim_null}
\alias{Weber_AML_sim_null_SE}
\alias{Weber_AML_sim_null_flowSet}
\alias{Weber_AML_sim_random_seeds}
\alias{Weber_AML_sim_random_seeds_5pc}
\alias{Weber_AML_sim_random_seeds_5pc_SE}
\alias{Weber_AML_sim_random_seeds_5pc_flowSet}
\alias{Weber_AML_sim_random_seeds_1pc}
\alias{Weber_AML_sim_random_seeds_1pc_SE}
\alias{Weber_AML_sim_random_seeds_1pc_flowSet}
\alias{Weber_AML_sim_random_seeds_0.1pc}
\alias{Weber_AML_sim_random_seeds_0.1pc_SE}
\alias{Weber_AML_sim_random_seeds_0.1pc_flowSet}
\alias{Weber_AML_sim_less_distinct}
\alias{Weber_AML_sim_less_distinct_5pc}
\alias{Weber_AML_sim_less_distinct_5pc_SE}
\alias{Weber_AML_sim_less_distinct_5pc_flowSet}
\alias{Weber_AML_sim_less_distinct_1pc}
\alias{Weber_AML_sim_less_distinct_1pc_SE}
\alias{Weber_AML_sim_less_distinct_1pc_flowSet}
\alias{Weber_AML_sim_less_distinct_0.1pc}
\alias{Weber_AML_sim_less_distinct_0.1pc_SE}
\alias{Weber_AML_sim_less_distinct_0.1pc_flowSet}
\title{
'Weber_AML_sim' semi-simulated datasets
}
\description{
Semi-simulated mass cytometry (CyTOF) datasets from Weber et al. (2019), constructed by 
computationally 'spiking in' small percentages of AML (acute myeloid leukemia) blast cells 
into samples of healthy BMMCs (bone marrow mononuclear cells), simulating the phenotype of 
minimal residual disease (MRD) in AML patients. These datasets can be used to benchmark 
differential analysis algorithms used to test for differentially abundant rare cell populations. 
Raw data sourced from Levine et al. (2015), and data generation strategy modified from 
Arvaniti et al. (2017). See Weber et al. (2019), Supplementary Note 1, for more details.
}
\details{
This is a set of semi-simulated mass cytometry (CyTOF) datasets, generated for benchmarking 
purposes in our paper introducing the 'diffcyt' framework (Weber et al., 2019).

The datasets are constructed by computationally 'spiking in' small percentages of AML (acute 
myeloid leukemia) blast cells into samples of healthy BMMCs (bone marrow mononuclear cells), 
simulating the phenotype of minimal residual disease (MRD) in AML patients. Blast cells are 
spiked in at 3 different thresholds of abundance (5\%, 1\%, and 0.1\%), to create multiple 
simulations with varying levels of difficulty.

These datasets can be used to benchmark differential analysis algorithms used to test for 
differentially abundant rare cell populations.

The raw data consists of 5 healthy samples (H1-H5), 1 AML (diseased) sample from condition CN 
(cytogenetically normal), and 1 AML sample from condition CBF (core-binding factor translocation). 
The dataset is constructed by splitting the healthy samples into 3 parts; one part is kept as 
the healthy reference condition, and small proportions of either CN or CBF blast cells are spiked 
into the other two parts to create the semi-simulated MRD conditions CN and CBF. We are then 
interested in detecting the differentially abundant rare population of either CN or CBF blasts 
in differential comparisons between CN and healthy, or CBF and healthy.

The dataset contains 16 surface protein markers used to define cell populations. For additional 
details, including numbers of cells per sample, see Weber et al. (2019), Supplementary Note 1 
(in particular Supplementary Tables 1 and 2).

Multiple simulations are available, as described in our paper (Weber et al., 2019). These are 
stored in the objects listed below.

In each case, the objects are available in both \code{\link{SummarizedExperiment}} and 
\code{\link{flowSet}} formats, with cells stored in rows, and protein markers in columns (i.e. 
the usual format for cytometry data). After loading the datasets, they can be inspected using the 
standard accessor functions for either \code{SummarizedExperiments} or \code{flowSets} (e.g. for 
\code{SummarizedExperiments}: \code{\link{rowData}}, \code{\link{colData}}, \code{\link{assays}}, 
and \code{\link{metadata}}).

For the \code{SummarizedExperiments}: \code{assays} contain tables of expression values 
(multiple \code{assays} for datasets with multiple replicates); 
\code{rowData} contains group IDs, patient IDs, sample IDs, and a column identifying spike-in cells; 
\code{colData} contains channel names, marker names, and marker classes; and 
\code{metadata} contains experiment information and number of cells.

For the \code{flowSets}: individual \code{flowFrames} within the \code{flowSet} contain tables 
of expression values (multiple \code{flowFrames} for datasets with multiple replicates); 
row data is stored as additional columns of numeric values within the expression tables; 
column data is stored in the \code{pData(parameters())} slot of the individual \code{flowFrames}; and 
additional information (e.g. experiment information, marker information, replicate information, 
and lookup tables to identify row data values) is stored in the \code{description()} slot of 
the \code{flowFrames}.

\strong{Main simulations}
\emph{
Separate files for each threshold of abundance (5\%, 1\%, and 0.1\%), as well additional objects 
containing all blast cells.
}
\itemize{
\item Weber_AML_sim_main_5pc_SE
\item Weber_AML_sim_main_5pc_flowSet
\item Weber_AML_sim_main_1pc_SE
\item Weber_AML_sim_main_1pc_flowSet
\item Weber_AML_sim_main_0.1pc_SE
\item Weber_AML_sim_main_0.1pc_flowSet
\item Weber_AML_sim_main_blasts_all_SE
\item Weber_AML_sim_main_blasts_all_flowSet
}

\strong{Additional simulations: null simulations}
\itemize{
\item Weber_AML_sim_null_SE
\item Weber_AML_sim_null_flowSet
}

\strong{Additional simulations: modified random seeds}
\emph{
Separate files for each threshold of abundance (5\%, 1\%, and 0.1\%).
}
\itemize{
\item Weber_AML_sim_random_seeds_5pc_SE
\item Weber_AML_sim_random_seeds_5pc_flowSet
\item Weber_AML_sim_random_seeds_1pc_SE
\item Weber_AML_sim_random_seeds_1pc_flowSet
\item Weber_AML_sim_random_seeds_0.1pc_SE
\item Weber_AML_sim_random_seeds_0.1pc_flowSet
}

\strong{Additional simulations: 'less distinct' spike-in cells}
\emph{
Separate files for each threshold of abundance (5\%, 1\%, and 0.1\%).
}
\itemize{
\item Weber_AML_sim_less_distinct_5pc_SE
\item Weber_AML_sim_less_distinct_5pc_flowSet
\item Weber_AML_sim_less_distinct_1pc_SE
\item Weber_AML_sim_less_distinct_1pc_flowSet
\item Weber_AML_sim_less_distinct_0.1pc_SE
\item Weber_AML_sim_less_distinct_0.1pc_flowSet
}

Note that prior to performing any downstream analyses, the expression values should be 
transformed. A standard transformation used for mass cytometry data is the \code{\link{asinh}} 
with \code{cofactor = 5}.

The raw data is sourced from Levine et al. (2015), and the data generation strategy is modified 
from Arvaniti et al. (2017). See Weber et al. (2019), Supplementary Note 1, for more details.

Original links to raw data from Cytobank:
- all cells (also contains gating scheme for CD34+CD45mid cells, i.e. blasts): 
https://community.cytobank.org/cytobank/experiments/46098/illustrations/121588
- blasts (repository cloned from the one for 'all cells' above, using the gating scheme 
for CD34+CD45mid cells; this allows .fcs files for the subset to be exported): 
https://community.cytobank.org/cytobank/experiments/63534/illustrations/125318

Data files are also available from FlowRepository (FR-FCM-ZYL8): http://flowrepository.org/id/FR-FCM-ZYL8
}
\usage{
Weber_AML_sim_main_5pc_SE(metadata = FALSE)
Weber_AML_sim_main_5pc_flowSet(metadata = FALSE)
Weber_AML_sim_main_1pc_SE(metadata = FALSE)
Weber_AML_sim_main_1pc_flowSet(metadata = FALSE)
Weber_AML_sim_main_0.1pc_SE(metadata = FALSE)
Weber_AML_sim_main_0.1pc_flowSet(metadata = FALSE)
Weber_AML_sim_main_blasts_all_SE(metadata = FALSE)
Weber_AML_sim_main_blasts_all_flowSet(metadata = FALSE)
Weber_AML_sim_null_SE(metadata = FALSE)
Weber_AML_sim_null_flowSet(metadata = FALSE)
Weber_AML_sim_random_seeds_5pc_SE(metadata = FALSE)
Weber_AML_sim_random_seeds_5pc_flowSet(metadata = FALSE)
Weber_AML_sim_random_seeds_1pc_SE(metadata = FALSE)
Weber_AML_sim_random_seeds_1pc_flowSet(metadata = FALSE)
Weber_AML_sim_random_seeds_0.1pc_SE(metadata = FALSE)
Weber_AML_sim_random_seeds_0.1pc_flowSet(metadata = FALSE)
Weber_AML_sim_less_distinct_5pc_SE(metadata = FALSE)
Weber_AML_sim_less_distinct_5pc_flowSet(metadata = FALSE)
Weber_AML_sim_less_distinct_1pc_SE(metadata = FALSE)
Weber_AML_sim_less_distinct_1pc_flowSet(metadata = FALSE)
Weber_AML_sim_less_distinct_0.1pc_SE(metadata = FALSE)
Weber_AML_sim_less_distinct_0.1pc_flowSet(metadata = FALSE)
}
\arguments{
\item{metadata}{\code{logical} value indicating whether ExperimentHub metadata (describing the 
overall dataset) should be returned only, or if the whole dataset should be loaded. Default = FALSE, 
which loads the whole dataset.}
}
\examples{
Weber_AML_sim_main_5pc_SE()
Weber_AML_sim_main_5pc_flowSet()
}
\value{Returns a \code{\link{SummarizedExperiment}} or \code{\link{flowSet}} object.}
\references{
Arvaniti and Claassen (2017), "Sensitive detection of rare disease-associated cell subsets via 
representation learning", Nature Communications, 8:14825: https://www.ncbi.nlm.nih.gov/pubmed/28382969

Levine et al. (2015), "Data-driven phenotypic dissection of AML reveals progenitor-like cells 
that correlate with prognosis", Cell, 162, 184-197: https://www.ncbi.nlm.nih.gov/pubmed/26095251

Weber et al. (2019). "diffcyt: Differential discovery in high-dimensional cytometry via high-resolution 
clustering." Communications Biology, 2:183: https://www.ncbi.nlm.nih.gov/pubmed/31098416
}
\keyword{datasets}
