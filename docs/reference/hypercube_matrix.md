# Generate A Hypercube Pooling Matrix

Construct the incidence matrix of a \\d\\-dimensional hypercube design
by arranging the individuals on a \\q \times \cdots \times q\\ lattice
and pooling along coordinate hyperplanes.

## Usage

``` r
hypercube_matrix(q, d)
```

## Arguments

- q:

  Integer scalar giving the side length of each coordinate dimension. It
  must satisfy \\q \ge 2\\. The returned design contains \\q\\ pools in
  each coordinate family.

- d:

  Integer scalar giving the number of dimensions. It must satisfy \\d
  \ge 2\\. Each individual belongs to exactly \\d\\ pools, one from each
  coordinate family.

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\, where \\N = q^d\\ and \\J = dq\\. Each row has constant size
\\q^{d-1}\\ and each column has constant sum \\d\\. The attached
`design` attribute is a named list with fields `type`, `N`, `M`, `q`,
`d`, `pool_size`, and `pools_per_sample`, where `M` is a legacy display
label for the number of pools \\J\\.

## Details

Label each individual by a coordinate vector \\(x_1,\ldots,x_d) \in
\\0,\ldots,q-1\\^d\\. For each coordinate index \\\ell \in
\\1,\ldots,d\\\\ and each value \\a \in \\0,\ldots,q-1\\\\, define the
pool \$\$\mathcal{P}\_{\ell,a} = \\i : x\_{\ell}(i) = a\\.\$\$ The
resulting matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\ has \$\$J =
dq\$\$ rows and \$\$N = q^d\$\$ columns. Every pool fixes one coordinate
and lets the remaining \\d-1\\ coordinates vary freely, so \$\$n_j =
\|\mathcal{P}\_j\| = q^{d-1}.\$\$ Every individual has one coordinate
value in each dimension, hence \$\$w_i = \|\mathcal{J}\_i\| = d.\$\$ For
\\d = 2\\, this recovers the square array design with equally sized row
and column families.

## References

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

Du, D.-Z., & Hwang, F. K. (2000). *Combinatorial Group Testing and Its
Applications* (2nd ed.). World Scientific.

## See also

[`array_matrix`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`array_matrix_3d`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
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
[`dorfman_optimal_g()`](https://myaseen208.github.io/ppgt/reference/dorfman_optimal_g.md),
[`hyper_ec_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md),
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`is_kts()`](https://myaseen208.github.io/ppgt/reference/is_kts.md),
[`kirkman_matrix()`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md),
[`list_designs()`](https://myaseen208.github.io/ppgt/reference/list_designs.md),
[`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md),
[`pg_matrix()`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
M2 <- hypercube_matrix(q = 4, d = 2)
M2
#> 8 x 16 sparse Matrix of class "dgCMatrix"
#>                                     
#> [1,] 1 . . . 1 . . . 1 . . . 1 . . .
#> [2,] . 1 . . . 1 . . . 1 . . . 1 . .
#> [3,] . . 1 . . . 1 . . . 1 . . . 1 .
#> [4,] . . . 1 . . . 1 . . . 1 . . . 1
#> [5,] 1 1 1 1 . . . . . . . . . . . .
#> [6,] . . . . 1 1 1 1 . . . . . . . .
#> [7,] . . . . . . . . 1 1 1 1 . . . .
#> [8,] . . . . . . . . . . . . 1 1 1 1
Matrix::rowSums(M2)
#> [1] 4 4 4 4 4 4 4 4
Matrix::colSums(M2)
#>  [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

M3 <- hypercube_matrix(q = 3, d = 3)
M3
#> 9 x 27 sparse Matrix of class "dgCMatrix"
#>                                                            
#>  [1,] 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . .
#>  [2,] . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 .
#>  [3,] . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1 . . 1
#>  [4,] 1 1 1 . . . . . . 1 1 1 . . . . . . 1 1 1 . . . . . .
#>  [5,] . . . 1 1 1 . . . . . . 1 1 1 . . . . . . 1 1 1 . . .
#>  [6,] . . . . . . 1 1 1 . . . . . . 1 1 1 . . . . . . 1 1 1
#>  [7,] 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . . . . . . . .
#>  [8,] . . . . . . . . . 1 1 1 1 1 1 1 1 1 . . . . . . . . .
#>  [9,] . . . . . . . . . . . . . . . . . . 1 1 1 1 1 1 1 1 1
```
