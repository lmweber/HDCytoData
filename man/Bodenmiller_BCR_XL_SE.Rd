\name{Bodenmiller_BCR_XL_SE}
\docType{data}
\alias{Bodenmiller_BCR_XL_SE}
\alias{Bodenmiller_BCR_XL_flowSet}
\title{'Bodenmiller_BCR_XL' data set}
\description{
Mass cytometry data set from Bodenmiller et al. (2012), consisting of paired samples of healthy 
PBMCs in unstimulated and BCR-XL stimulated condition.
}
\details{
Mass cytometry data set from Bodenmiller et al. (2012), consisting of paired samples of healthy 
peripheral blood mononuclear cells (PBMCs), where one sample from each pair was stimulated with 
B cell receptor / Fc receptor cross-linker (BCR-XL). This induces a strong differential signal 
for several differential state markers. The strongest observed signal is for phosphorylated 
S6 (pS6) in B cells.

Contains 16 samples (8 paired samples); 24 protein markers, consisting of 10 'cell type' 
markers and 14 'cell state' markers; 172,791 cells in total.

This data set has previously been used for benchmarking evaluations by Bruggner et al. (2014) 
and Nowicka et al. (2017).

Raw data files are available from Cytobank, experiment 15713 
(https://community.cytobank.org/cytobank/experiments/15713/download_files; see also Citrus wiki, 
section 'PBMC Example 1': https://github.com/nolanlab/citrus/wiki/PBMC-Example-1).

The population labels (which can be used to compare with clustering results) are sourced from 
Nowicka et al. (2017), where they were generated using a strategy of expert-guided manual merging 
of automatically generated clusters from the FlowSOM algorithm (Van Gassen et al., 2015).

The data set is provided in two formats: 'SummarizedExperiment' and 'flowSet'. In each case, 
cells are stored in rows, and protein markers in columns (this is the usual format used for flow 
cytometry data).

For the 'SummarizedExperiment', row and column meta-data can be accessed with 
\code{rowData(Bodenmiller_BCR_XL_SE)} and \code{colData(Bodenmiller_BCR_XL_SE)}. The row meta-data 
contains group IDs, patient IDs, sample IDs, and population IDs. The column meta-data contains 
protein marker names, and logical vectors indicating all markers, 'cell type' markers, and 'cell 
state' markers. The expression values are formatted as a single table (i.e. the 16 samples are 
stacked together into one table of values).

For the 'flowSet', the row meta-data is stored as additional columns within the 'flowFrame' 
object for each sample. The protein marker names are stored in the column names of the tables of values 
(it is not possible to include additional information in the column names to identify the 
'cell type' and 'cell state' markers). The expression values are formatted into one table per 
sample, where each table represents a 'flowFrame' within the overall 'flowSet' object.
}
\usage{
Bodenmiller_BCR_XL_SE(metadata = FALSE)
Bodenmiller_BCR_XL_flowSet(metadata = FALSE)
}
\arguments{
\itemize{
\item{metadata}{\code{logical} value indicating whether metadata only should be returned, 
or if the resource (data set) should be loaded. Default behavior (\code{metadata = FALSE}) 
loads the data set.}
}
}
\examples{
Bodenmiller_BCR_XL_SE()
Bodenmiller_BCR_XL_flowSet()
}
\value{A 'SummarizedExperiment' or 'flowSet' object.}
\references{
Bodenmiller et al. (2012); Nowicka et al. (2017)
}
\keyword{datasets}
