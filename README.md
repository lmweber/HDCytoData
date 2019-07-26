# HDCytoData

[![Build Status](https://travis-ci.org/lmweber/HDCytoData.svg?branch=master)](https://travis-ci.org/lmweber/HDCytoData)


## Summary

The `HDCytoData` package is an extensible resource containing a set of publicly available high-dimensional flow cytometry and mass cytometry (CyTOF) benchmark datasets, which have been formatted into `SummarizedExperiment` and `flowSet` Bioconductor object formats. The data objects are hosted on Bioconductor's `ExperimentHub` platform.

The objects each contain one or more tables of cell-level expression values, as well as all required metadata. Row metadata includes sample IDs, group IDs, patient IDs, reference cell population or cluster labels (where available), and labels identifying 'spiked in' cells (where available). Column metadata includes channel names, protein marker names, and protein marker classes (cell type or cell state).

Note that raw expression values should be transformed prior to any downstream analyses.

Currently, the package includes benchmark datasets used in our previous work to evaluate methods for clustering and differential analyses. The datasets are provided here in `SummarizedExperiment` and `flowSet` formats in order to make them easier to access and integrate into R/Bioconductor workflows.


## Details

For additional details and an example showing how to load the datasets, see the package vignette available from [Bioconductor](http://bioconductor.org/packages/HDCytoData).

For details on the individual datasets, see the help files within the package.


## Availability and installation

The `HDCytoData` package is freely available from [Bioconductor](http://bioconductor.org/packages/HDCytoData), and can be installed by following standard Bioconductor package installation procedures:

```{r}
# install BiocManager (if not already installed)
install.packages("BiocManager")

# install HDCytoData package
BiocManager::install("HDCytoData")
```

