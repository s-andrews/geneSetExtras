# Gene Set Extras

This package contains a number of functions which supplement the functionality found in the ```clusterProfiler``` package.  These add some useful plots and functionality for visualising and simplifying gene set results.

In addition the package provides functions for performing a differential enrichment analysis of gene set data containing multiple sets.

## Installation

The package can be installed using the ```remotes``` or ```devtools``` packages.

Either:

```
install.packages("remotes")

remotes::install_github("s-andrews/geneSetExtras")
```

or

```
install.packages("devtools")

devtools::install_github("s-andrews/geneSetExtras")
```

## Input Requirements
The package currently supports gene set enrichment analysis from named lists of genes.  The inputs will be a set of gene set enrichment results from ```clusterProfiler``` calculated with functions such as ```enrichGO```.  The input for these will simply be a list of 'hit' genes as well as a 'background' list of all genes measured.

## Usage

The package is still evolving, so the easiest way to see the current capabilities and usage is to look at the [geneSetExtras Vignette](https://html-preview.github.io/?url=https://raw.githubusercontent.com/s-andrews/geneSetExtras/refs/heads/main/doc/gene_set_extras_usage.html) for details of how to use the package.