# HDCytoData

[![Build Status](https://travis-ci.org/lmweber/HDCytoData.svg?branch=master)](https://travis-ci.org/lmweber/HDCytoData)


## Summary

Data package containing high-dimensional cytometry data sets saved in Bioconductor object formats, hosted on Bioconductor ExperimentHub.

This package contains a set of publicly available high-dimensional flow cytometry and mass cytometry (CyTOF) data sets, which have been formatted into the `SummarizedExperiment` and `flowSet` Bioconductor object formats. The objects contain the cell-level expression values, as well as row and column meta-data, including sample IDs, group IDs, true cell population or cluster labels (where available), and protein marker names.

These data sets have been used in our previous work and publications for benchmarking purposes, e.g. to evaluate the performance of clustering algorithms. They are provided here in the `SummarizedExperiment` and `flowSet` formats to make them easier to access for ourselves and other method developers.


## Details

For additional details, including references and raw data sources, see the help files for each data set.


## Tutorial and example workflow

For a short example workflow demonstrating how to load the data objects and use them in an analysis workflow, see the package vignette.


## Availability and installation

The `HDCytoData` package is available as an experiment data package from Bioconductor. It can be installed using the Bioconductor installer (`biocLite`):

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("HDCytoData")
```


