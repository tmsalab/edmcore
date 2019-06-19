
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ecdmcore

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/tmsalab/ecdmcore.svg?branch=master)](https://travis-ci.com/tmsalab/ecdmcore)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=2\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/ecdmcore)](https://cran.r-project.org/package=ecdmcore)
[![CRAN
status](https://www.r-pkg.org/badges/version/ecdmcore)](https://cran.r-project.org/package=ecdmcore)
[![RStudio CRAN Mirror’s Monthly
Downloads](http://cranlogs.r-pkg.org/badges/ecdmcore?color=brightgreen)](http://www.r-pkg.org/pkg/ecdmcore)
[![RStudio CRAN Mirror’s Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/ecdmcore?color=brightgreen)](http://www.r-pkg.org/pkg/ecdmcore)
[![Codecov test
coverage](https://codecov.io/gh/tmsalab/ecdmcore/branch/master/graph/badge.svg)](https://codecov.io/gh/tmsalab/ecdmcore?branch=master)
<!-- badges: end -->

The goal of `ecdmcore` is to house a set of functions shared by many
packages within the exploratory cognitive diagnostic modeling framework.

## Installation

You can install the released version of ecdmcore from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ecdmcore")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tmsalab/ecdmcore")
```

## Usage

To use `ecdmcore`, load the package using:

``` r
library("ecdmcore")
```

## Overview

The package contains class structure shared between different estimation
units.

In particular, we have:

  - Attributes: `attribute_classes()`, `attribute_bijection()`,
    `attribute_inv_bijection()`.
  - Matrix: `q_matrix()`

## Authors

James Joseph Balamuta, Steven Andrew Culpepper, and Jeffrey Douglas

## Citing the `ecdmcore` package

To ensure future development of the package, please cite `ecdmcore`
package if used during the analysis or simulations. Citation information
for the package may be acquired by using in *R*:

``` r
citation("ecdmcore")
```

## License

GPL (\>= 2)
