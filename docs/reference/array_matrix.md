# Generate A Two-Dimensional Array Pooling Matrix

Construct a row-column array design in which the \\N = rc\\ individuals
are placed on an \\r \times c\\ grid and pooled by grid rows and
columns.

## Usage

``` r
array_matrix(r, c)
```

## Arguments

- r:

  Integer scalar giving the number of grid rows. It must satisfy \\r \ge
  2\\. This parameter determines the first family of pools and
  contributes \\r\\ rows to \\\mathbf{M}\\.

- c:

  Integer scalar giving the number of grid columns. It must satisfy \\c
  \ge 2\\. This parameter determines the second family of pools and
  contributes \\c\\ rows to \\\mathbf{M}\\.

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\, where \\J = r + c\\ and \\N = rc\\. Every individual belongs to
exactly one row pool and one column pool, so \\w_i = \|\mathcal{J}\_i\|
= 2\\. The attached `design` attribute is a named list with fields
`type`, `N`, `M`, `rows`, `cols`, `pool_size_row`, `pool_size_col`, and
`pools_per_sample`, where `M` is a legacy display label for the number
of pools \\J\\.

## Details

Index the individuals by ordered pairs \\(u,v)\\ with \\u \in
\\1,\ldots,r\\\\ and \\v \in \\1,\ldots,c\\\\. After linearizing the
grid into column indices \\i \in \\1,\ldots,N\\\\, the returned design
matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\ has two pool families:
\$\$\mathcal{P}^{\mathrm{row}}\_u = \\i : \text{individual } i \text{
lies in row } u\\, \qquad u = 1,\ldots,r,\$\$
\$\$\mathcal{P}^{\mathrm{col}}\_v = \\i : \text{individual } i \text{
lies in column } v\\, \qquad v = 1,\ldots,c.\$\$ Hence \$\$J = r + c,
\qquad n_j \in \\r, c\\, \qquad w_i = 2.\$\$ For a single positive
individual, the unique positive row and unique positive column identify
that individual by intersection. With multiple positives, the stage-1
positive pattern may correspond to several candidate intersections, so
the design is primarily useful as a low-weight screening matrix rather
than an exact combinatorial code.

## References

Phatarfod, R. M., & Sudbury, A. (1994). The use of a square array scheme
in blood testing. *Statistics in Medicine*, 13(22), 2337-2343.

Kim, H. Y., Hudgens, M. G., Dreyfuss, J. M., Westreich, D. J., &
Pilcher, C. D. (2007). Comparison of group testing algorithms for case
identification in the presence of test error. *Biometrics*, 63(4),
1152-1163.

## See also

[`array_matrix_3d`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
[`hypercube_matrix`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md),
[`compare_designs`](https://myaseen208.github.io/ppgt/reference/compare_designs.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
[`bibd_matrix()`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
[`compare_designs()`](https://myaseen208.github.io/ppgt/reference/compare_designs.md),
[`dorfman_matrix()`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md),
[`dorfman_optimal_g()`](https://myaseen208.github.io/ppgt/reference/dorfman_optimal_g.md),
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
M <- array_matrix(r = 4, c = 5)
M
#> 9 x 20 sparse Matrix of class "dgCMatrix"
#>                                              
#>  [1,] 1 1 1 1 1 . . . . . . . . . . . . . . .
#>  [2,] . . . . . 1 1 1 1 1 . . . . . . . . . .
#>  [3,] . . . . . . . . . . 1 1 1 1 1 . . . . .
#>  [4,] . . . . . . . . . . . . . . . 1 1 1 1 1
#>  [5,] 1 . . . . 1 . . . . 1 . . . . 1 . . . .
#>  [6,] . 1 . . . . 1 . . . . 1 . . . . 1 . . .
#>  [7,] . . 1 . . . . 1 . . . . 1 . . . . 1 . .
#>  [8,] . . . 1 . . . . 1 . . . . 1 . . . . 1 .
#>  [9,] . . . . 1 . . . . 1 . . . . 1 . . . . 1
Matrix::rowSums(M)
#> [1] 5 5 5 5 4 4 4 4 4
Matrix::colSums(M)
#>  [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

M96 <- array_matrix(r = 8, c = 12)
M96
#> 20 x 96 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,] 1 1 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [2,] . . . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . . .
#>  [3,] . . . . . . . . . . . . . . . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 .
#>  [4,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 1
#>  [5,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [6,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [7,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [8,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [9,] 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1
#> [10,] . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . .
#> [11,] . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . .
#> [12,] . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . .
#> [13,] . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . .
#> [14,] . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . .
#> [15,] . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . .
#> [16,] . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . .
#> [17,] . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . .
#> [18,] . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . .
#> [19,] . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . .
#> [20,] . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 .
#>                                                                                
#>  [1,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [2,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [3,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [4,] 1 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [5,] . . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . . . .
#>  [6,] . . . . . . . . . . . . . . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1 1 . .
#>  [7,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 1 1
#>  [8,] . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#>  [9,] . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 .
#> [10,] 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1
#> [11,] . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . .
#> [12,] . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . .
#> [13,] . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . .
#> [14,] . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . .
#> [15,] . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . .
#> [16,] . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . .
#> [17,] . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . .
#> [18,] . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . .
#> [19,] . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . . .
#> [20,] . . . . . . . . . . 1 . . . . . . . . . . . 1 . . . . . . . . . . . 1 . .
#>                                                  
#>  [1,] . . . . . . . . . . . . . . . . . . . . . .
#>  [2,] . . . . . . . . . . . . . . . . . . . . . .
#>  [3,] . . . . . . . . . . . . . . . . . . . . . .
#>  [4,] . . . . . . . . . . . . . . . . . . . . . .
#>  [5,] . . . . . . . . . . . . . . . . . . . . . .
#>  [6,] . . . . . . . . . . . . . . . . . . . . . .
#>  [7,] 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . .
#>  [8,] . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1 1
#>  [9,] . . . . . . . . . . 1 . . . . . . . . . . .
#> [10,] . . . . . . . . . . . 1 . . . . . . . . . .
#> [11,] 1 . . . . . . . . . . . 1 . . . . . . . . .
#> [12,] . 1 . . . . . . . . . . . 1 . . . . . . . .
#> [13,] . . 1 . . . . . . . . . . . 1 . . . . . . .
#> [14,] . . . 1 . . . . . . . . . . . 1 . . . . . .
#> [15,] . . . . 1 . . . . . . . . . . . 1 . . . . .
#> [16,] . . . . . 1 . . . . . . . . . . . 1 . . . .
#> [17,] . . . . . . 1 . . . . . . . . . . . 1 . . .
#> [18,] . . . . . . . 1 . . . . . . . . . . . 1 . .
#> [19,] . . . . . . . . 1 . . . . . . . . . . . 1 .
#> [20,] . . . . . . . . . 1 . . . . . . . . . . . 1
```
