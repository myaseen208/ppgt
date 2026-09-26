# ppgt: Pooled Group Testing Design Toolkit

[![](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![CRAN
status](https://www.r-pkg.org/badges/version/ppgt)](https://CRAN.R-project.org/package=ppgt)
[![](https://cranlogs.r-pkg.org/badges/grand-total/ppgt)](https://cran.r-project.org/package=ppgt)

## Overview

**ppgt** is an R package for constructing, analyzing, and comparing
pooled group testing designs. It brings together a set of classical and
modern designs for research and experimentation, including polynomial
pooling, PBEST, HYPER, Kirkman triple systems, other combinatorial
designs (BIBDs, hypercubes, arrays, projective planes), and
separable/disjunct constructions, together with finite-field (Galois
field) utilities and COMP/GPSR decoders.

## Installation

`ppgt` is not yet on CRAN. Install the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("myaseen208/ppgt")

# or with remotes
# install.packages("remotes")
remotes::install_github("myaseen208/ppgt")
```

## Quick start

``` r

library(ppgt)

# Construct a Polynomial Pools (PP) matrix -- the standard PP-based
# P-BEST configuration (q = 8, d = 3, nl = 6, N = 384)
M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
dim(M)
#> [1]  48 384

# Simulate two positive individuals
x_true <- rep(0L, 384)
x_true[c(72, 142)] <- 1L

# Generate pooled outcomes and decode
z <- as.integer(as.vector(M %*% x_true) > 0)
result <- pp_decode(M, z)
#> 
#> === PP Decode ===
#> Matrix: 48 pools x 384 samples
#> Positive pools: 11/48
#> COMP candidates: 2
#> Result: 72, 142
#> Method: COMP (exact)
result$positives
#> [1]  72 142
```

## Kirkman triple system (KTS) helpers

`kirkman_matrix(v)` constructs a resolvable Steiner triple system on `v`
points – a Kirkman triple system exists exactly when `v %% 6 == 3`.
Every returned matrix is independently re-verified with
[`is_kts()`](https://myaseen208.github.io/ppgt/reference/is_kts.md)
before it is returned, checking all four properties that define a KTS:
block size 3, correct block count, every pair of points covered exactly
once, and resolvability into `(v-1)/2` parallel classes that each
partition the points:

``` r

# v = 15 is the classic Kirkman schoolgirl problem
M15 <- kirkman_matrix(15)
dim(M15)                                  # 35 triples x 15 points
#> [1] 35 15
is_kts(M15)                               # verified: TRUE
#> [1] TRUE
#> attr(,"reasons")
#> character(0)
attr(M15, "design")$parallel_class        # which of the 7 "days" each triple is in
#>  [1] 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7
```

Any `v` that is a power of 3 (`3, 9, 27, 81, ...`, via the affine
geometry AG(n,3)) and `v = 15` (the classical solution built from the
lines of PG(3,2)) construct quickly and deterministically. Other
admissible `v` (e.g. `21`, `33`, `63`) fall back to a general randomized
search that is not tuned for larger `v` – it is not guaranteed to
succeed even for `v = 21` – in which case
[`kirkman_matrix()`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md)
raises an informative error rather than ever returning an unverified
result.

## Package contents

### Core designs

- [`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
  /
  [`pp_decode()`](https://myaseen208.github.io/ppgt/reference/pp_decode.md)
  /
  [`pp_design()`](https://myaseen208.github.io/ppgt/reference/pp_design.md)
  /
  [`pp_verify()`](https://myaseen208.github.io/ppgt/reference/pp_verify.md)
  – Polynomial Pools construction, decoding, and design search
- [`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
  /
  [`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md)
  – HYPER-style pooling matrices
- [`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md)
  /
  [`hyper_ec_decode()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_decode.md)
  – package-specific experimental parity-augmented extension of HYPER
- [`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)
  – Tapestry-style pooling matrix
- [`SeparableMatrix()`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md)
  /
  [`DisjunctMatrix()`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md)
  – d-separable / d-disjunct constructions, self-verified against their
  defining properties
- [`dorfman_matrix()`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md)
  /
  [`dorfman_optimal_g()`](https://myaseen208.github.io/ppgt/reference/dorfman_optimal_g.md)
  – classic two-stage Dorfman pooling
- [`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md)
  /
  [`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md)
  – row/column and 3D array designs
- [`bibd_matrix()`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md)
  – balanced incomplete block designs (built-in for specific
  `(v, k, lambda)`)
- [`hypercube_matrix()`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md)
  – d-dimensional hypercube designs
- [`kirkman_matrix()`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md)
  / [`is_kts()`](https://myaseen208.github.io/ppgt/reference/is_kts.md)
  – Kirkman triple systems, see above
- [`pg_matrix()`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md)
  – projective plane PG(2, q) designs
- [`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md)
  – clinical 384-to-94 P-BEST heuristic layout
- [`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md)
  /
  [`score_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/score_pooling_matrix.md)
  – prescribed-marginal regular pooling matrices

### Finite-field utilities

- [`gf()`](https://myaseen208.github.io/ppgt/reference/gf.md),
  [`gf_add()`](https://myaseen208.github.io/ppgt/reference/gf_add.md),
  [`gf_mult()`](https://myaseen208.github.io/ppgt/reference/gf_mult.md),
  [`gf_pow()`](https://myaseen208.github.io/ppgt/reference/gf_pow.md) –
  Galois field construction and arithmetic
- [`is_prime()`](https://myaseen208.github.io/ppgt/reference/is_prime.md),
  [`is_prime_power()`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md)

### Workflow and comparison

- [`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
  – unified constructor across all design families
- [`simulate_group_testing()`](https://myaseen208.github.io/ppgt/reference/simulate_group_testing.md),
  [`run_testing_workflow()`](https://myaseen208.github.io/ppgt/reference/run_testing_workflow.md),
  [`protocol_summary()`](https://myaseen208.github.io/ppgt/reference/protocol_summary.md)
- [`compare_designs()`](https://myaseen208.github.io/ppgt/reference/compare_designs.md),
  [`list_designs()`](https://myaseen208.github.io/ppgt/reference/list_designs.md),
  [`compare_all_designs()`](https://myaseen208.github.io/ppgt/reference/compare_all_designs.md),
  [`compare_all_designs_honest()`](https://myaseen208.github.io/ppgt/reference/compare_all_designs_honest.md)

## Documentation

- Package reference:
  [github.com/myaseen208/ppgt](https://github.com/myaseen208/ppgt)
- Notation registry:
  [`ppgt_notations()`](https://myaseen208.github.io/ppgt/reference/ppgt_notations.md)
- Vignettes: `ppgt_theory` (overview), and topic articles `ppgt_pp`,
  `ppgt_pp_188x282` (worked 188 x 282 pooling matrix example),
  `ppgt_pbest`, `ppgt_hyper`, `ppgt_hyperec`, `ppgt_galois_field`,
  `ppgt_separable_disjunct`, `ppgt_comparisons`, `ppgt_diagnostics`

## License

GPL-3

## Author and citation

Package author and maintainer: **Muhammad Yaseen**
(<myaseen208@gmail.com>).

**Authors/Contributors**: Christopher McMahan, Christopher Bilder,
Joshua Tebbs, and Pranta Das.

``` R
Yaseen, M., McMahan, C.,  Bilder, C., Tebbs, J. and Das, P. (2026). ppgt: Pooled Group Testing Design Toolkit.
R package version 0.1.2. https://github.com/myaseen208/ppgt
```
