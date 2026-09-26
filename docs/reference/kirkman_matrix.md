# Generate A Kirkman Triple System Pooling Matrix

Construct the incidence matrix of a built-in Kirkman triple system,
viewed as a resolvable Steiner triple system for pooled testing.

## Usage

``` r
kirkman_matrix(v)
```

## Arguments

- v:

  Integer scalar giving the number of individuals. A Kirkman triple
  system exists only when \\v \equiv 3 \pmod 6\\; this function returns
  `NULL` with a warning when that condition fails or when no built-in
  construction is available.

## Value

Either a sparse binary matrix of class `dgCMatrix` with dimensions \\J
\times N\\, where \\N = v\\ and \\J = v(v-1)/6\\, or `NULL`. When a
matrix is returned, each row has size three and each column sum is
\\(v-1)/2\\. The attached `design` attribute stores `type`, `v`, `b`,
`r`, and `k`.

## Details

A Kirkman triple system is a Steiner triple system whose blocks can be
partitioned into parallel classes. In the returned incidence matrix
\\\mathbf{M} \in \\0,1\\^{J \times N}\\, the columns index the \\N=v\\
individuals and the rows index the triples. The defining combinatorial
relations are \$\$n_j = \|\mathcal{P}\_j\| = 3 \quad \text{for all }
j,\$\$ \$\$w_i = \|\mathcal{J}\_i\| = \frac{v-1}{2} \quad \text{for all
} i,\$\$ and every unordered pair of distinct individuals appears
together in exactly one pool. Consequently, \$\$J =
\frac{v(v-1)}{6}.\$\$ The necessary and sufficient existence condition
for a Kirkman triple system is \$\$v \equiv 3 \pmod 6.\$\$ This
implementation contains explicit block sets for a small number of values
of \\v\\; other admissible values may still return `NULL`.

## References

Colbourn, C. J., & Rosa, A. (1999). *Triple Systems*. Oxford University
Press.

Colbourn, C. J., & Dinitz, J. H. (Eds.). (2007). *Handbook of
Combinatorial Designs* (2nd ed.). Chapman & Hall/CRC.

## See also

[`bibd_matrix`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
[`pg_matrix`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
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
[`hypercube_matrix()`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md),
[`list_designs()`](https://myaseen208.github.io/ppgt/reference/list_designs.md),
[`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md),
[`pg_matrix()`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
M9 <- kirkman_matrix(9)
M9
#> 12 x 9 sparse Matrix of class "dgCMatrix"
#>                        
#>  [1,] 1 1 1 . . . . . .
#>  [2,] . . . 1 1 1 . . .
#>  [3,] . . . . . . 1 1 1
#>  [4,] 1 . . 1 . . 1 . .
#>  [5,] . 1 . . 1 . . 1 .
#>  [6,] . . 1 . . 1 . . 1
#>  [7,] 1 . . . 1 . . . 1
#>  [8,] . 1 . . . 1 1 . .
#>  [9,] . . 1 1 . . . 1 .
#> [10,] 1 . . . . 1 . 1 .
#> [11,] . 1 . 1 . . . . 1
#> [12,] . . 1 . 1 . 1 . .
Matrix::rowSums(M9)
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3
Matrix::colSums(M9)
#> [1] 4 4 4 4 4 4 4 4 4

M10 <- kirkman_matrix(10)
#> Warning: KTS(v) requires v == 3 (mod 6). Got v =10
M10
#> NULL
```
