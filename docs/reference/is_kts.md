# Verify the Defining Properties of a Kirkman Triple System

Independently checks whether a putative Kirkman triple system (KTS)
incidence matrix, such as one returned by
[`kirkman_matrix`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md),
actually satisfies all four properties that define a KTS: block size
three, existence condition on \\v\\, every pair of points covered
exactly once, and resolvability into parallel classes that each
partition the point set.

## Usage

``` r
is_kts(M, verbose = FALSE)
```

## Arguments

- M:

  A candidate KTS incidence matrix: a \\J \times N\\ binary matrix
  (`matrix`, `dgCMatrix`, or other `Matrix` class) whose rows are blocks
  and whose columns are the \\N\\ points. To check resolvability, `M`
  should carry a `design` attribute with a `parallel_class` component (a
  length-\\J\\ vector of class labels), exactly as attached by
  [`kirkman_matrix`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md);
  without it, resolvability is reported as a failure rather than
  skipped, since an unresolved block set is not a Kirkman triple system.

- verbose:

  Logical; if `TRUE` and the check fails, print the failure reasons via
  [`message()`](https://rdrr.io/r/base/message.html). Default `FALSE`.

## Value

A single `logical` scalar, `TRUE` only if every one of the four
properties holds. On `FALSE`, the result carries a character-vector
attribute `"reasons"` describing every violation found; pass
`verbose = TRUE` to also print them.

## Details

Writing \\v = N\\ for the number of points and \\b = J\\ for the number
of blocks, this function checks, directly against the incidence matrix
(not against any internal bookkeeping the constructor may have used):

1.  \\v \equiv 3 \pmod 6\\ and \\b = v(v-1)/6\\;

2.  every row has exactly 3 nonzero entries;

3.  every unordered pair of columns co-occurs (both 1) in exactly one
    row;

4.  the rows, grouped by `attr(M, "design")$parallel_class`, form
    exactly \\(v-1)/2\\ classes of \\v/3\\ rows each, and within each
    class the rows partition \\\\1,\ldots,v\\\\ (every point covered
    exactly once, none missing, none repeated).

## References

Ray-Chaudhuri, D. K., & Wilson, R. M. (1971). Solution of Kirkman's
schoolgirl problem. *Proceedings of Symposia in Pure Mathematics*, 19,
187-203.

## See also

[`kirkman_matrix`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md)

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
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
M15 <- kirkman_matrix(15)
is_kts(M15)
#> [1] TRUE
#> attr(,"reasons")
#> character(0)

# A deliberately broken KTS(15) (the pre-fix package bug): a hand-built
# block set whose last three parallel classes repeat some points and
# omit others, so several pairs are covered zero or more than one time.
broken_blocks <- list(
  c(1,2,3), c(4,5,6), c(7,8,9), c(10,11,12), c(13,14,15),
  c(1,4,7), c(2,5,8), c(3,6,9), c(10,13,14), c(11,12,15),
  c(1,5,9), c(2,6,7), c(3,4,8), c(10,12,14), c(11,13,15),
  c(1,6,8), c(2,4,9), c(3,5,7), c(10,11,14), c(12,13,15),
  c(1,10,15), c(2,11,13), c(3,12,14), c(4,8,13), c(5,7,12),
  c(1,11,14), c(2,10,12), c(3,13,15), c(4,7,15), c(5,8,10),
  c(1,12,13), c(2,14,15), c(3,10,11), c(4,9,14), c(5,6,13)
)
rows_i <- unlist(lapply(seq_along(broken_blocks), function(i) rep(i, 3)))
cols_j <- unlist(broken_blocks)
M_broken <- Matrix::sparseMatrix(
  i = rows_i, j = cols_j, x = rep(1L, length(rows_i)), dims = c(35, 15)
)
attr(M_broken, "design") <- list(parallel_class = rep(1:7, each = 5))
is_kts(M_broken)
#> [1] FALSE
#> attr(,"reasons")
#> [1] "43 pair(s) of points are covered 0 or >1 times (each pair must occur in exactly one block); 24 pair(s) never covered, 19 pair(s) covered more than once"
#> [2] "parallel class '5' does not partition the 15 points exactly once"                                                                                       
#> [3] "parallel class '6' does not partition the 15 points exactly once"                                                                                       
#> [4] "parallel class '7' does not partition the 15 points exactly once"                                                                                       
attr(is_kts(M_broken), "reasons")
#> [1] "43 pair(s) of points are covered 0 or >1 times (each pair must occur in exactly one block); 24 pair(s) never covered, 19 pair(s) covered more than once"
#> [2] "parallel class '5' does not partition the 15 points exactly once"                                                                                       
#> [3] "parallel class '6' does not partition the 15 points exactly once"                                                                                       
#> [4] "parallel class '7' does not partition the 15 points exactly once"                                                                                       
```
