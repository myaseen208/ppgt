# Construct a Smart Regular Pooling Matrix

Construct a sparse binary pooling matrix with prescribed constant row
and column sums. By default, an exact polynomial-pools construction is
returned when available; otherwise an overlap-aware regular construction
is used.

## Usage

``` r
regular_pooling_matrix(
  n_pools,
  n_samples,
  pool_size,
  pools_per_sample,
  seed = NULL,
  max_tries = 1000,
  method = c("auto", "pp", "regular", "overlap_optimized"),
  prefer_pp = TRUE,
  d_values = 2:8
)
```

## Arguments

- n_pools:

  Integer number of pools (rows).

- n_samples:

  Integer number of samples (columns).

- pool_size:

  Required number of samples in every pool.

- pools_per_sample:

  Required number of pools containing every sample.

- seed:

  Optional integer seed for randomized fallback construction.

- max_tries:

  Maximum randomized candidates for the optimized fallback.

- method:

  One of `"auto"`, `"pp"`, `"regular"`, or `"overlap_optimized"`.
  Automatic mode tries PP first. PP mode requires an exact PP
  construction. Regular mode uses one generic realization.

- prefer_pp:

  In automatic mode, try an exact PP construction first?

- d_values:

  Polynomial dimensions tried, in order, during PP detection.

## Value

A sparse binary matrix with the requested dimensions and marginals. An
exact PP result is returned directly and retains its `"pp"` attribute; a
generic result has a `"design"` attribute.

## Details

Feasibility requires \$\$n\_{pools} pool\\size = n\_{samples}
pools\\per\\sample,\$\$ `pool_size <= n_samples`, and
`pools_per_sample <= n_pools`.

PP detection derives `q = n_pools / pools_per_sample` and
`nl = pools_per_sample`, calls
[`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
for each requested `d`, and verifies all dimensions and row and column
sums. The generic fallback keeps exact marginals and ranks candidates by
maximum specimen overlap, number of specimen pairs with overlap greater
than one, and maximum pool overlap.

Generic fallback matrices do not automatically claim PP, P-BEST,
disjunctness, separability, or protocol-specific decoding guarantees.
Use
[`score_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/score_pooling_matrix.md)
or other diagnostics to evaluate overlaps.

## See also

[`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md),
[`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
[`bibd_matrix()`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
[`compare_designs()`](https://myaseen208.github.io/ppgt/reference/compare_designs.md),
[`dorfman_matrix()`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md),
[`dorfman_optimal_g()`](https://myaseen208.github.io/ppgt/reference/dorfman_optimal_g.md),
[`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md),
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`hypercube_matrix()`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md),
[`kirkman_matrix()`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md),
[`list_designs()`](https://myaseen208.github.io/ppgt/reference/list_designs.md),
[`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md),
[`pg_matrix()`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
M1 <- regular_pooling_matrix(94, 188, 8, 4, seed = 4L)
M2 <- regular_pooling_matrix(186, 372, 8, 4, seed = 4L)
M3 <- regular_pooling_matrix(188, 282, 6, 4, seed = 4L)
attr(M1, "pp")
#> NULL
attr(M2, "pp")
#> NULL
attr(M3, "pp")
#> $q
#> [1] 47
#> 
#> $d
#> [1] 3
#> 
#> $nl
#> [1] 4
#> 
#> $N
#> [1] 282
#> 
#> $M
#> [1] 188
#> 
#> $pool_size
#> [1] 6
#> 
#> $k_max
#> [1] 1
#> 
M3_pp <- pp_matrix(q = 47, d = 3, nl = 4, N = 282)
identical(as.matrix(M3), as.matrix(M3_pp))
#> [1] TRUE
score_pooling_matrix(M1)
#> $dim
#> [1]  94 188
#> 
#> $row_sum_table
#> rs
#>  8 
#> 94 
#> 
#> $col_sum_table
#> cs
#>   4 
#> 188 
#> 
#> $max_col_overlap
#> [1] 1
#> 
#> $n_col_pairs_overlap_gt_1
#> [1] 0
#> 
#> $col_overlap_table
#> col_positive
#>    1 
#> 5264 
#> 
#> $max_pool_overlap
#> [1] 1
#> 
#> $pool_overlap_table
#> pool_positive
#>    1 
#> 2256 
#> 
#> $density
#> [1] 0.04255319
#> 
#> $has_pp_attr
#> [1] FALSE
#> 
score_pooling_matrix(M2)
#> $dim
#> [1] 186 372
#> 
#> $row_sum_table
#> rs
#>   8 
#> 186 
#> 
#> $col_sum_table
#> cs
#>   4 
#> 372 
#> 
#> $max_col_overlap
#> [1] 1
#> 
#> $n_col_pairs_overlap_gt_1
#> [1] 0
#> 
#> $col_overlap_table
#> col_positive
#>     1 
#> 10416 
#> 
#> $max_pool_overlap
#> [1] 1
#> 
#> $pool_overlap_table
#> pool_positive
#>    1 
#> 4464 
#> 
#> $density
#> [1] 0.02150538
#> 
#> $has_pp_attr
#> [1] FALSE
#> 
score_pooling_matrix(M3)
#> $dim
#> [1] 188 282
#> 
#> $row_sum_table
#> rs
#>   6 
#> 188 
#> 
#> $col_sum_table
#> cs
#>   4 
#> 282 
#> 
#> $max_col_overlap
#> [1] 1
#> 
#> $n_col_pairs_overlap_gt_1
#> [1] 0
#> 
#> $col_overlap_table
#> col_positive
#>    1 
#> 5640 
#> 
#> $max_pool_overlap
#> [1] 1
#> 
#> $pool_overlap_table
#> pool_positive
#>    1 
#> 3384 
#> 
#> $density
#> [1] 0.0212766
#> 
#> $has_pp_attr
#> [1] TRUE
#> 
```
