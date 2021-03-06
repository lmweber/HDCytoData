# HDCytoData

[![Build Status](https://travis-ci.org/lmweber/HDCytoData.svg?branch=master)](https://travis-ci.org/lmweber/HDCytoData)


## Summary

The `HDCytoData` package is an extensible resource containing a set of publicly available high-dimensional flow cytometry and mass cytometry (CyTOF) benchmark datasets, which have been formatted into `SummarizedExperiment` and `flowSet` Bioconductor object formats. The data objects are hosted on Bioconductor's `ExperimentHub` platform.

The objects each contain one or more tables of cell-level expression values, as well as all required metadata. Row metadata includes sample IDs, group IDs, patient IDs, reference cell population or cluster labels (where available), and labels identifying 'spiked in' cells (where available). Column metadata includes channel names, protein marker names, and protein marker classes (cell type, cell state, as well as non protein marker columns).

Note that raw expression values should be transformed prior to any downstream analyses.

Currently, the package includes benchmark datasets used in our previous work to evaluate methods for clustering and differential analyses. The datasets are provided here in `SummarizedExperiment` and `flowSet` formats in order to make them easier to access and integrate into R/Bioconductor workflows.

For more details, see our paper describing the `HDCytoData` package:

- [Weber L.M. and Soneson C. (2019), *HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats*, F1000Research, 8:1459, v2.](https://f1000research.com/articles/8-1459)


## Vignettes

Additional details are provided in the following vignettes, available from the [Bioconductor](http://bioconductor.org/packages/HDCytoData) website:

- [HDCytoData package](http://bioconductor.org/packages/release/data/experiment/vignettes/HDCytoData/inst/doc/HDCytoData_package.html): Overview of the package and example showing how to load the datasets
- [Examples and use cases](http://bioconductor.org/packages/release/data/experiment/vignettes/HDCytoData/inst/doc/Examples_and_use_cases.html): Examples and use cases
- [Contribution guidelines](http://bioconductor.org/packages/release/data/experiment/vignettes/HDCytoData/inst/doc/Contribution_guidelines.html): Guidelines for contributing new datasets

For details on the datasets, see the help files for each dataset available within the package, or the metadata from the `ExperimentHub` database (see "HDCytoData package" vignette).


## Availability and installation

The `HDCytoData` package is freely available from [Bioconductor](http://bioconductor.org/packages/HDCytoData), and can be installed by following standard Bioconductor package installation procedures:

```{r}
# install BiocManager (if not already installed)
install.packages("BiocManager")

# install HDCytoData package
BiocManager::install("HDCytoData")
```


## Citation

Please cite the following paper:

- [Weber L.M. and Soneson C. (2019), *HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats*, F1000Research, 8:1459, v2.](https://f1000research.com/articles/8-1459)


