# Polynomial Pools Algorithm for Group Testing

Implementation of the Polynomial Pools (PP) algorithm for constructing
efficient group testing pooling matrices.

## Details

Group testing pools multiple individuals to reduce the number of assays.
The Polynomial Pools construction uses arithmetic over
\\\mathrm{GF}(q)\\ to build sparse binary pooling matrices with
controlled overlap.

The package uses PP both as a general design family and as the
matrix-level source of the standard PP-based P-BEST configuration \\(q,
d, n_l, N) = (8, 3, 6, 384)\\.

A central structural property is the intersection bound: any two
distinct individuals share at most \\d - 1\\ pools. This is the matrix
property behind the standard PP COMP detection bound.

## Relationship To Other Components

- Finite fields:

  [`gf`](https://myaseen208.github.io/ppgt/reference/gf.md) and related
  helpers provide the arithmetic used by the PP construction.

- PP matrices:

  [`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
  builds the incidence matrix \\\mathbf{M}\\ for the chosen PP family.

- PP-based P-BEST:

  The standard PP-based P-BEST family is one specific PP configuration,
  documented separately in the package vignettes.

## References

Tan, Y. H. I. (2020). Pooling Matrix Designs for Group Testing. *SIAM
Undergraduate Research Online*.

Brust, D., & Brust, J. J. (2023). Effective Matrix Designs for COVID-19
Group Testing. *BMC Bioinformatics*, 24, 26.

## See also

[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md) for creating
Galois Fields,
[`pp_decode`](https://myaseen208.github.io/ppgt/reference/pp_decode.md)
for decoding test results

Other polynomial_pools:
[`pp_design()`](https://myaseen208.github.io/ppgt/reference/pp_design.md),
[`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md),
[`pp_verify()`](https://myaseen208.github.io/ppgt/reference/pp_verify.md)
