# Compute An Approximate Optimal Dorfman Pool Size

Choose the pool size \\g\\ that minimizes the expected number of tests
per individual under the homogeneous Dorfman model.

## Usage

``` r
dorfman_optimal_g(p)
```

## Arguments

- p:

  Numeric scalar giving the common individual prevalence \\p_i = p \in
  (0,1)\\. It must satisfy \\0 \< p \< 1\\.

## Value

An integer scalar giving the pool size \\g\\ that minimizes the
objective evaluated by this function over a local search range around
\\1 / \sqrt{p}\\.

## Details

In two-stage Dorfman testing, a pool of size \\g\\ requires one stage-1
pooled test and then \\g\\ follow-up individual tests whenever the
pooled result is positive. Under a homogeneous prevalence model
\\\Pr(\widetilde{y}\_i = 1) = p\\ with independent individual statuses
and perfect testing, the expected number of tests per individual is
\$\$\psi(g; p) = \frac{1}{g} + 1 - (1-p)^g.\$\$ This function
numerically minimizes \\\psi(g; p)\\ over integer values in a finite
window centered near the classical approximation \$\$g\_{\mathrm{opt}}
\approx \frac{1}{\sqrt{p}}.\$\$ The returned value can be passed
directly to
[`dorfman_matrix`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md).

## References

Dorfman, R. (1943). The detection of defective members of large
populations. *The Annals of Mathematical Statistics*, 14(4), 436-440.

Sobel, M., & Groll, P. A. (1959). Group testing to eliminate efficiently
all defectives in a binomial sample. *Bell System Technical Journal*,
38(5), 1179-1252.

## See also

[`dorfman_matrix`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md),
[`compare_designs`](https://myaseen208.github.io/ppgt/reference/compare_designs.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
[`bibd_matrix()`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
[`compare_designs()`](https://myaseen208.github.io/ppgt/reference/compare_designs.md),
[`dorfman_matrix()`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md),
[`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md),
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`hypercube_matrix()`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md),
[`is_kts()`](https://myaseen208.github.io/ppgt/reference/is_kts.md),
[`kirkman_matrix()`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md),
[`list_designs()`](https://myaseen208.github.io/ppgt/reference/list_designs.md),
[`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md),
[`pg_matrix()`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
g1 <- dorfman_optimal_g(0.01)
g1
#> [1] 11
g2 <- dorfman_optimal_g(0.05)
g2
#> [1] 5
g3 <- dorfman_optimal_g(0.001)
g3
#> [1] 32
```
