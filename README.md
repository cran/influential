
<!-- README.md is generated from README.Rmd. Please edit that file -->

# influential

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/asalavaty/influential?branch=master&svg=true)](https://ci.appveyor.com/project/asalavaty/influential)
[![](https://www.r-pkg.org/badges/version/influential?color=blue)](https://cran.r-project.org/package=influential)
[![](http://cranlogs.r-pkg.org/badges/grand-total/influential?color=green)](https://cran.r-project.org/package=influential)
[![Rdoc](http://www.rdocumentation.org/badges/version/influential)](http://www.rdocumentation.org/packages/influential)
[![](https://img.shields.io/badge/Integrated%20Value%20of%20Influence-IVI-blue.svg)](https://doi.org/10.1016/j.patter.2020.100052)
[![](https://img.shields.io/badge/SIR--based%20Influence%20Ranking-SIRIR-green.svg)](https://doi.org/10.1016/j.patter.2020.100052)
[![](https://img.shields.io/badge/Experimental--data--based%20Integrative%20Ranking-ExIR-blue.svg)](https://github.com/asalavaty/influential)
<!-- badges: end -->

![The influential R package
logo](https://github.com/asalavaty/influential/blob/master/logo.png)

## Overview

The goal of `influential` is to help identification of the most
`influential` nodes in a network as well as the classification and
ranking of top candidate features. This package contains functions for
the classification and ranking of features, reconstruction of networks
from adjacency matrices and data frames, analysis of the topology of the
network and calculation of centrality measures as well as a novel and
powerful `influential` node ranking. The **Experimental-data-based
Integrative Ranking (ExIR)** is a sophisticated model for classification
and ranking of the top candidate features based on only the experimental
data. The first integrative method, namely the **Integrated Value of
Influence (IVI)**, that captures all topological dimensions of the
network for the identification of network most `influential` nodes is
also provided as a function. Also, neighborhood connectivity, H-index,
local H-index, and collective influence (CI), all of which required
centrality measures for the calculation of **IVI**, are for the first
time provided in an R package. Additionally, a function is provided for
running **SIRIR** model, which is the combination of leave-one-out cross
validation technique and the conventional SIR model, on a network to
unsupervisedly rank the true influence of vertices. Furthermore, some
functions have been provided for the assessment of dependence and
correlation of two network centrality measures as well as the
conditional probability of deviation from their corresponding means in
opposite directions.

Check out [**our paper**](https://doi.org/10.1016/j.patter.2020.100052)
for a more complete description of the IVI formula and all of its
underpinning methods and analyses.

## Author

The `influential` package was written by [Adrian (Abbas)
Salavaty](https://www.AbbasSalavaty.com)

## Advisors

Mirana Ramialison and Peter D. Currie

## How to Install

You can install the official [CRAN
release](https://cran.r-project.org/package=influential) of the
`influential` with the following code:

``` r
install.packages("influential")
```

Or the development version from GitHub:

``` r
## install.packages("devtools")
devtools::install_github("asalavaty/influential")
```

## Vignettes

[Detailed description of the functions and their
outputs](https://github.com/asalavaty/influential/blob/master/vignettes/Vignettes.md)

You may browse Vignettes from within R using the following code.

``` r
browseVignettes("influential")
```

## An Example for Calculation of IVI

This is a basic example which shows you how to solve a common problem:

``` r
library(influential)

MyData <- centrality.measures # A data frame of centrality measures

# This function calculates the Integrated Value of Influence (IVI)
My.vertices.IVI <- ivi.from.indices(DC = centrality.measures$DC,       # Calculation of IVI
                                   CR = centrality.measures$CR,
                                   NC = centrality.measures$NC,
                                   LH_index = centrality.measures$LH_index,
                                   BC = centrality.measures$BC,
                                   CI = centrality.measures$CI)

print(head(My.vertices.IVI))
#> [1] 24.670056  8.344337 18.621049  1.017768 29.437028 33.512598
```

## How to cite `influential`

To cite `influential`, please cite the [**associated
paper**](https://doi.org/10.1016/j.patter.2020.100052). You can also
refer to the package’s citation information using the citation()
function.

``` r
citation("influential")
```

## How to contribute

Please don’t hesitate to report any bugs/issues and request for
enhancement or any other contributions. To submit a bug report or
enhancement request, please use the [`influential` GitHub issues
tracker](https://github.com/asalavaty/influential/issues).
