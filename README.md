# HDCytoData

[![Build Status](https://travis-ci.org/lmweber/HDCytoData.svg?branch=master)](https://travis-ci.org/lmweber/HDCytoData)


## Summary

Data package containing high-dimensional cytometry datasets saved in Bioconductor object formats, hosted on Bioconductor ExperimentHub.

This package contains a set of publicly available high-dimensional flow cytometry and mass cytometry (CyTOF) datasets, which have been formatted into the `SummarizedExperiment` and `flowSet` Bioconductor object formats. The objects contain the cell-level expression values, as well as row and column metadata, including sample IDs, group IDs, true cell population labels or cluster labels (where available), channel names, protein marker names, and protein marker classes (cell type or cell state).

These datasets have been used in our previous work and publications for benchmarking purposes, e.g. to benchmark clustering algorithms or methods for differential analysis. They are provided here in the `SummarizedExperiment` and `flowSet` formats to make them easier to access for ourselves and other method developers.


## Details

For additional details, including references and raw data sources, see the help files for each dataset.


## Tutorial

A short tutorial showing how to load the data objects is included in the package vignette.


## Availability and installation

The `HDCytoData` package is available as an experiment data package from [Bioconductor](http://bioconductor.org/packages/HDCytoData). It can be installed using `BiocManager`:

```{r}
# install BiocManager from CRAN (if not already installed)
install.packages("BiocManager")

# install HDCytoData package
BiocManager::install("HDCytoData")
```

