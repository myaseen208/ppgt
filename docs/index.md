# ppgt: Pooled Group Testing Design Toolkit

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/myaseen208/ppgt)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://myaseen208.github.io/ppgt)

## Overview

**ppgt** is an R package for constructing, analyzing, and comparing
pooled group testing designs. It brings together a set of classical and
modern designs for research and experimentation, including polynomial
pooling, PBEST, HYPER, separable and disjunct constructions, and
finite-field utilities.

The package is intended for reproducible method development and
documentation-driven workflows. The code, tests, and documentation
follow the package notation defined in `Notations.tex` and are validated
with
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html),
[`devtools::test()`](https://devtools.r-lib.org/reference/test.html),
and
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html).
Use
[`ppgt_notations()`](https://myaseen208.github.io/ppgt/reference/ppgt_notations.md)
to inspect the package notation registry directly.

## Features

- Polynomial pooling design construction and decoding
- PP-based P-BEST design support
- HYPER-style design support
- Separable and disjunct matrix constructors
- GF(q) field construction and arithmetic utilities
- Design comparison helpers
- Roxygen2 documentation and pkgdown-ready articles

## Installation

### From GitHub

``` r

devtools::install_github("myaseen208/ppgt")
```

### From source

``` r

install.packages("ppgt_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

``` r

library(ppgt)

# Construct a PP matrix
M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
M

# Simulate positives
x_true <- rep(0L, 384)
x_true[c(72, 142)] <- 1L
x_true[c(72, 142)]

# Generate pooled outcomes
z <- as.integer(as.vector(M %*% x_true) > 0)
z

# Decode
result <- pp_decode(M, z)
result
result$positives
```

## Local development

You can work with the package locally without installing it:

``` r

devtools::load_all("ppgt")
```

Then run tests and checks:

``` r

devtools::test()
devtools::check()
```

## Package contents

### Core designs

- [`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
  and
  [`pp_decode()`](https://myaseen208.github.io/ppgt/reference/pp_decode.md)
- [`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
- [`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md)
  and
  [`hyper_ec_decode()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_decode.md)
- [`SeparableMatrix()`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md)
  and
  [`DisjunctMatrix()`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md)

### Finite-field utilities

- [`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md)
- [`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md)
- [`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md)
- [`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md)
- [`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md)
- [`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md)

### Comparison tools

- [`compare_all_designs()`](https://myaseen208.github.io/ppgt/reference/compare_all_designs.md)
- [`compare_all_designs_honest()`](https://myaseen208.github.io/ppgt/reference/compare_all_designs_honest.md)
- `compare_*()` helper functions

## Design families

### Polynomial Pooling (PP)

Finite-field-based pooling designs with decoding support and clear
parameterization in terms of `q`, `d`, `nl`, and `N`.

### PBEST

PBEST is supported in the package through the PP-based 48-pool design
configuration and related package documentation.

### HYPER

HYPER is supported through an upstream-compatible matrix generator and
related package documentation.

### HyperEC

[`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md)
and
[`hyper_ec_decode()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_decode.md)
are documented as a **package-specific experimental parity-augmented
extension**. They are not presented as a direct implementation of the
Hong et al. HYPER procedure.

### Separable and disjunct designs

The package includes deterministic and randomized separable/disjunct
constructions together with verification helpers.

## Documentation

- Package reference: pkgdown site
- Notation registry:
  [`ppgt_notations()`](https://myaseen208.github.io/ppgt/reference/ppgt_notations.md)
- Intro article: `ppgt_theory`
- Topic articles:
  - `ppgt_pp`
  - `ppgt_pbest`
  - `ppgt_hyper`
  - `ppgt_hyperec`
  - `ppgt_galois_field`
  - `ppgt_separable_disjunct`
  - `ppgt_comparisons`
  - `ppgt_diagnostics`
- Workflow guide: package documentation and roxygen comments

## Testing and validation

Current package validation includes:

- [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
- [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
- [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)

The package is maintained so that these commands run cleanly in the
release workflow.

## Release notes

### ppgt 0.1.0

- Version updated to `0.1.0`
- Package metadata and documentation prepared for a professional release
- HYPER-EC material re-scoped as an experimental package extension
- Finite-field helpers and core design constructors documented and
  tested

## References

- Dorfman, R. (1943). The detection of defective members of large
  populations.
- Hong, F. et al. (2022). HYPER pooled testing results.
- Zismanov, V. et al. (2024). PBEST design paper.
- Additional references are listed in the package articles and reference
  manual.

## License

GPL-3

## Author

Muhammad Yaseen
