# HyperDesign Compatibility Wrapper

Provide a compatibility wrapper with the same argument interface as the
former `HyperDesign::HyperDesign()` constructor. This function forwards
directly to
[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
and returns the same matrix object.

## Usage

``` r
HyperDesign(n, m, q, reorder = TRUE)
```

## Arguments

- n:

  Integer scalar giving the total number of individuals. This becomes
  the column dimension \\N\\ of the returned matrix.

- m:

  Integer scalar giving the total number of pools. This becomes the row
  dimension \\J\\ of the returned matrix.

- q:

  Integer scalar giving the number of pools containing each individual.
  Allowed values are `1L`, `2L`, and `3L`.

- reorder:

  Logical scalar controlling whether the columns are reordered by the
  upstream HYPER grouping convention. The default is `TRUE`.

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\. The object is exactly the output of
[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)`(n = n, m = m, q = q, reorder = reorder)`
and carries the same `design` attribute.

## Details

This wrapper exists for interface compatibility only. If \\\mathbf{M}
\in \\0,1\\^{J \times N}\\ denotes the returned matrix, then \\J = m\\,
\\N = n\\, and every individual \\i\\ satisfies \$\$w_i =
\|\mathcal{J}\_i\| = q.\$\$ All arithmetic constraints on `m` and `q`
are inherited unchanged from
[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md):
\$\$q \in \\1,2,3\\,\$\$ with \\m\\ even for \\q = 2\\, and \\m \equiv 0
\pmod 6\\ together with primality of \\m-1\\ for \\q = 3\\. No
additional computation is carried out by this wrapper beyond forwarding
the call.

## References

Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
method. *Nature Communications*, 13, Article 3626.

## See also

[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`compare_all_designs`](https://myaseen208.github.io/ppgt/reference/compare_all_designs.md),
[`compare_all_designs_honest`](https://myaseen208.github.io/ppgt/reference/compare_all_designs_honest.md)

Other pooling_designs:
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
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
design1 <- HyperDesign(n = 10, m = 5, q = 1)
design1
#> 5 x 10 sparse Matrix of class "dgCMatrix"
#>                         
#> [1,] 1 1 . . . . . . . .
#> [2,] . . 1 1 . . . . . .
#> [3,] . . . . 1 1 . . . .
#> [4,] . . . . . . 1 1 . .
#> [5,] . . . . . . . . 1 1

design2 <- HyperDesign(n = 12, m = 6, q = 2)
design2
#> 6 x 12 sparse Matrix of class "dgCMatrix"
#>                             
#> [1,] 1 . . . 1 . . . 1 . . 1
#> [2,] . 1 . 1 . . . 1 . . . 1
#> [3,] . . 1 . 1 . 1 . . . 1 .
#> [4,] . . 1 . . 1 . 1 . 1 . .
#> [5,] . 1 . . . 1 . . 1 . 1 .
#> [6,] 1 . . 1 . . 1 . . 1 . .
Matrix::colSums(design2)
#>  [1] 2 2 2 2 2 2 2 2 2 2 2 2
```
