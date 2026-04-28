# StatMatch R pakage
Marcello D’Orazio

# StatMatch <img src="https://www.r-project.org/logo/Rlogo.png" align="right" width="80"/>

> **Statistical Matching or Data Fusion**

[![CRAN
version](https://www.r-pkg.org/badges/version/StatMatch.png)](https://CRAN.R-project.org/package=StatMatch)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/StatMatch.png)](https://CRAN.R-project.org/package=StatMatch)
[![Mentioned in Awesome Official
Statistics](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)

------------------------------------------------------------------------

## Overview

**StatMatch** is an R package for integrating two data sources referring
to the same target population that share a number of common variables.
Integration can be performed at:

-   **micro level** — producing a synthetic (fused) file
-   **macro level** — estimating correlation coefficients, regression
    coefficients, contingency tables, etc.

Some functions can also be used to **impute missing values** through hot
deck imputation methods.

------------------------------------------------------------------------

## Installation

Install the stable version from CRAN:

``` r
install.packages("StatMatch")
```

------------------------------------------------------------------------

## Main Features

-   **Non-parametric hot deck imputation**: random, rank, and
    nearest-neighbour methods
-   **Predictive mean matching**: a mixed parametric-nonparametric
    procedure
-   **Complex survey designs**: methods for data from surveys with
    unequal weights
-   **Uncertainty analysis**: tools to explore the uncertainty inherent
    in the statistical matching framework

------------------------------------------------------------------------

## Documentation & Tutorials

Additional material including vignettes, tutorials, and presentations is
available in the
[`Tutorials_Vignette_OtherDocs`](https://github.com/marcellodo/StatMatch/tree/master/Tutorials_Vignette_OtherDocs)
folder.

------------------------------------------------------------------------

## References

-   D’Orazio M., Di Zio M., Scanu M. (2006) *Statistical Matching,
    Theory and Practice*. Wiley, Chichester.
-   D’Orazio M. (2015) Integration and imputation of survey data in R:
    the StatMatch package, *Romanian Statistical Review*, 2/2015,
    pp. 57–68.

------------------------------------------------------------------------

## Links

-   📦 [CRAN page](https://CRAN.R-project.org/package=StatMatch)
-   📄 [Reference manual
    (PDF)](https://cran.r-project.org/web/packages/StatMatch/StatMatch.pdf)
-   📋 [Official Statistics & Survey Statistics task
    view](https://cran.r-project.org/web/views/OfficialStatistics.html)
