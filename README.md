# `bivnegbincp`

> **Bayesian Multiple Change Point Detection for Bivariate Negative Binomial Data in R**

<p align="center">

**Bayesian Inference • Multiple Change Points • Binary Segmentation • Genomics • Rcpp Optimized**

</p>

---

## Overview

`bivnegbincp` is an R package for Bayesian multiple change point detection in sequences of **bivariate Negative Binomial (BNB)** random vectors.

The package implements the methodology developed in the accompanying dissertation and supports Bayesian inference under both

* Constant prior
* Dirichlet prior

Multiple change points are detected recursively using **Binary Segmentation**.

The methodology is designed for sequences of correlated count observations where abrupt structural changes may occur over time or along an ordered sequence. The package is applicable to a wide range of problems involving bivariate count data, including epidemiology, genomics, finance, environmental science, reliability engineering, and other longitudinal or sequential studies.

---

## Key Features

* Bayesian single change point detection
* Multiple change point detection via Binary Segmentation
* Constant prior formulation
* Dirichlet prior formulation
* Posterior probability estimation
* Interactive 3D visualization using Plotly
* Utilities for index correction
* Optimized computational backend using Rcpp
* Numerically stable implementation using `lgamma()`

---

## Installation

Install the development version directly from GitHub.

```r
# install.packages("remotes")

remotes::install_github("asinha639/bivnegbincp")

library(bivnegbincp)
```

---

## Methodology

The package implements Bayesian change point detection for observations

```text
(X₁,Y₁), (X₂,Y₂), … , (Xₙ,Yₙ)
```

where each observation follows a bivariate Negative Binomial distribution.

Posterior probabilities are computed for every admissible change point location.

Multiple change points are then obtained using recursive Binary Segmentation.

The underlying Bayesian methodology follows the dissertation exactly.

The optimized implementation changes **only the numerical evaluation** of the posterior and **does not modify the statistical model**.

---

## Available Functions

| Function               | Description                            |
| ---------------------- | -------------------------------------- |
| `bivnegbin_singlecp()` | Bayesian single change point detection |
| `BinSeg_BNB()`         | Multiple change point detection        |
| `plot_cp_3d()`         | Interactive Plotly visualization       |
| `index_correction()`   | Convert local to global indices        |

---

# Example 1 — Constant Prior

```r
library(bivnegbincp)

x1 <- c(12,13,15,18,17,65,70,69,72)
x2 <- c( 8,10,11,12,10,40,42,39,43)

result <- bivnegbin_singlecp(
    x1,
    x2,
    prior = "constant"
)

result
```

---

# Example 2 — Dirichlet Prior

```r
result <- bivnegbin_singlecp(

    x1,
    x2,

    prior = "dirichlet"

)

result
```

---

# Example 3 — Multiple Change Points

```r
cp <- BinSeg_BNB(

    x1,
    x2,

    prior = "constant"

)

cp
```

The returned object contains

* posterior probability
* local change point index
* corrected global index
* recursion level
* segment information

---

# Example 4 — Interactive Visualization

```r
plot_cp_3d(

    x1,
    x2,
    cp

)
```

The interactive Plotly visualization displays

* First count variable
* Second count variable
* Detected change points
* Posterior probabilities

Users can freely rotate, zoom and inspect detected change points in three dimensions.

---

# General Workflow

```r
library(bivnegbincp)

cp <- BinSeg_BNB(

    x1,
    x2,

    prior="dirichlet"

)

plot_cp_3d(

    x1,
    x2,
    cp

)
```

---

# Returned Output

Typical output contains

| Variable              | Description                                       |
| --------------------- | ------------------------------------------------- |
| posterior_probability | Posterior probability of a change point           |
| cp_local_index        | Estimated change point within the current segment |
| cp_global_index       | Estimated change point in the original sequence   |
| cp_status             | Indicator of detected change point                |
| additional metadata   | User-supplied indexing variables, if available    |

---

# Computational Optimization

Beginning with Version 2, the computational core has been substantially optimized while preserving the Bayesian methodology.

The following numerical improvements are incorporated:

* Log-scale computation
* `lgamma()` reformulation of Gamma products
* Prefix sums (`cumsum`)
* Removal of repeated summations
* Rcpp implementation of computational kernels
* Reduced memory allocation
* Improved numerical stability

These optimizations **do not change the underlying Bayesian model**.

Instead, they provide a mathematically equivalent evaluation of the posterior while substantially reducing execution time for large genomic datasets.

---

# Applications

The package is intended for Bayesian change point analysis of sequential bivariate count data arising in a wide range of scientific disciplines, including

* epidemiology
* genomics
* finance
* environmental monitoring
* industrial process control
* reliability engineering
* ecology
* quality control

Any ordered sequence of paired count observations satisfying the model assumptions can be analyzed using the implemented methodology.

---

# Package Structure

```
R/
    bivar_negbin_cp.R
    binseg_bnb.R
    plot_cp_3d.R
    index_correction.R

src/
    bnb_singlecp.cpp
    bnb_helpers.cpp
```

# Model Assumptions

The implemented methodology assumes

- observations are ordered in sequence,
- observations within each segment arise from a common bivariate Negative Binomial model,
- changes occur at unknown locations separating homogeneous segments,
- Bayesian inference is performed independently within each recursively identified segment during Binary Segmentation.

Users should ensure these assumptions are appropriate for their application before interpreting detected change points.

---

# Citation

If you use this package in published work, please cite the accompanying dissertation and any subsequent journal publication describing the Bayesian methodology.

---

# License

GPL (>=3)

---

# Development

Source code

https://github.com/asinha639/bivnegbincp

Issues and feature requests are welcome through the GitHub issue tracker.

---

## Author

**Arnoneel Sinha, PhD**

Center for Biotechnology & Genomic Medicine
Augusta University

---

## Future Development

Planned enhancements include

* additional visualization tools
* simulation utilities
* expanded genomic workflows
* benchmarking studies
* CRAN-ready release
* extended documentation and vignettes

---

**bivnegbincp** provides a reproducible, statistically rigorous, and computationally efficient framework for Bayesian multiple change point detection in bivariate Negative Binomial data.
