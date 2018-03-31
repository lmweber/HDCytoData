# HDCytoData

## Overview

This package contains a set of publicly available high-dimensional cytometry data sets, formatted into `SummarizedExperiment` and `flowSet` objects. The objects contain expression values as well as meta-data including sample IDs, group IDs, population labels, and protein marker names.

These data sets have been used in our previous publications for benchmarking evaluations. They are provided here as `SummarizedExperiments` and `flowSets` to make them easier to access for other method developers.

The raw data sets are also available from FlowRepository (see help files for each data set).


## How to load data

The data objects are stored on Bioconductor's ExperimentHub. To load the data objects, first install the `HDCytoData` package:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("HDCytoData")
```

Then, the data sets can be loaded by referring to the object names:

```{r}
library("HDCytoData")

# 'SummarizedExperiment' format
Bodenmiller_BCR_XL_SE()

# 'flowSet' format
Bodenmiller_BCR_XL_flowSet()
```

Alternatively, the data sets can be loaded via ExperimentHub. See the package vignette for details.

