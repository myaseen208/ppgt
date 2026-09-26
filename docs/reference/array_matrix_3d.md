# Generate A Three-Dimensional Array Pooling Matrix

Construct a three-dimensional array design in which the \\N = d_1 d_2
d_3\\ individuals are arranged on a rectangular lattice and pooled by
the three coordinate directions.

## Usage

``` r
array_matrix_3d(d1, d2, d3)
```

## Arguments

- d1:

  Integer scalar giving the size of the first coordinate dimension. It
  must satisfy \\d_1 \ge 2\\. This contributes \\d_1\\ pools to the
  returned matrix.

- d2:

  Integer scalar giving the size of the second coordinate dimension. It
  must satisfy \\d_2 \ge 2\\. This contributes \\d_2\\ pools to the
  returned matrix.

- d3:

  Integer scalar giving the size of the third coordinate dimension. It
  must satisfy \\d_3 \ge 2\\. This contributes \\d_3\\ pools to the
  returned matrix.

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\, where \\J = d_1 + d_2 + d_3\\ and \\N = d_1 d_2 d_3\\. Each
individual lies in exactly one pool from each coordinate family, so
\\w_i = \|\mathcal{J}\_i\| = 3\\. The attached `design` attribute is a
named list with fields `type`, `N`, `M`, `dimensions`, and
`pools_per_sample`, where `M` is a legacy display label for the number
of pools \\J\\.

## Details

Label the individuals by triples \\(u,v,w)\\ with \\u \in
\\1,\ldots,d_1\\\\, \\v \in \\1,\ldots,d_2\\\\, and \\w \in
\\1,\ldots,d_3\\\\. The design matrix \\\mathbf{M} \in \\0,1\\^{J \times
N}\\ is formed from three parallel pool families:
\$\$\mathcal{P}^{(1)}\_u = \\i : \text{the first coordinate of } i
\text{ is } u\\,\$\$ \$\$\mathcal{P}^{(2)}\_v = \\i : \text{the second
coordinate of } i \text{ is } v\\,\$\$ \$\$\mathcal{P}^{(3)}\_w = \\i :
\text{the third coordinate of } i \text{ is } w\\.\$\$ Therefore \$\$J =
d_1 + d_2 + d_3\$\$ and every individual participates in exactly three
pools. The family sizes are \$\$\|\mathcal{P}^{(1)}\_u\| = d_2 d_3,
\qquad \|\mathcal{P}^{(2)}\_v\| = d_1 d_3, \qquad
\|\mathcal{P}^{(3)}\_w\| = d_1 d_2.\$\$ When there is a single positive
individual, the unique positive pool in each family identifies its three
coordinates.

## References

Li, X., & Ying, K. (2013). Two- and three-dimensional pooling algorithms
for high-throughput screening assays. *Technometrics*, 55(3), 285-296.

Kim, H. Y., Hudgens, M. G., Dreyfuss, J. M., Westreich, D. J., &
Pilcher, C. D. (2007). Comparison of group testing algorithms for case
identification in the presence of test error. *Biometrics*, 63(4),
1152-1163.

## See also

[`array_matrix`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`hypercube_matrix`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md),
[`compare_designs`](https://myaseen208.github.io/ppgt/reference/compare_designs.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
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
M <- array_matrix_3d(3, 3, 2)
M
#> 8 x 18 sparse Matrix of class "dgCMatrix"
#>                                         
#> [1,] 1 1 1 1 1 1 . . . . . . . . . . . .
#> [2,] . . . . . . 1 1 1 1 1 1 . . . . . .
#> [3,] . . . . . . . . . . . . 1 1 1 1 1 1
#> [4,] 1 1 . . . . 1 1 . . . . 1 1 . . . .
#> [5,] . . 1 1 . . . . 1 1 . . . . 1 1 . .
#> [6,] . . . . 1 1 . . . . 1 1 . . . . 1 1
#> [7,] 1 . 1 . 1 . 1 . 1 . 1 . 1 . 1 . 1 .
#> [8,] . 1 . 1 . 1 . 1 . 1 . 1 . 1 . 1 . 1
Matrix::rowSums(M)
#> [1] 6 6 6 6 6 6 9 9
Matrix::colSums(M)
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
```
