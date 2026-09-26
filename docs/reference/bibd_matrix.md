# Generate A BIBD Pooling Matrix

Construct the incidence matrix of a balanced incomplete block design
(BIBD) when one of the built-in block constructions is available.

## Usage

``` r
bibd_matrix(v, k, lambda = 1L)
```

## Arguments

- v:

  Integer scalar giving the number of treatments in the underlying BIBD.
  In the pooling interpretation this is the total number of individuals,
  so the returned matrix has \\N = v\\ columns.

- k:

  Integer scalar giving the block size of the BIBD. In pooling terms,
  each pool \\\mathcal{P}\_j\\ has size \\k\\. It must satisfy \\2 \le k
  \< v\\.

- lambda:

  Integer scalar giving the pairwise concurrence parameter \\\lambda \ge
  1\\. Every unordered pair of distinct individuals appears together in
  exactly \\\lambda\\ pools when the design exists. The default is `1L`.

## Value

Either a sparse binary matrix of class `dgCMatrix` with dimensions \\J
\times N\\, where \\N = v\\ and \\J = b\\, or `NULL` if the required
arithmetic conditions fail or no built-in construction is available.
When a matrix is returned, the `design` attribute stores `type`, `v`,
`b`, `r`, `k`, and `lambda`.

## Details

A \\(v, b, r, k, \lambda)\\ balanced incomplete block design is
represented by an incidence matrix \\\mathbf{M} \in \\0,1\\^{b \times
v}\\. In the pooling interpretation, \\N = v\\ and \\J = b\\. The
defining relations are \$\$n_j = \|\mathcal{P}\_j\| = k \quad \text{for
all } j,\$\$ \$\$w_i = \|\mathcal{J}\_i\| = r \quad \text{for all }
i,\$\$ and every unordered pair of distinct individuals occurs together
in exactly \\\lambda\\ pools. The standard parameter identities are
\$\$bk = vr\$\$ and \$\$\lambda(v-1) = r(k-1).\$\$ This function first
computes \\r\\ and \\b\\ from these identities and then dispatches to a
small collection of built-in constructions such as the Fano plane. If no
matching construction is coded, the function returns `NULL`.

## References

Beth, T., Jungnickel, D., & Lenz, H. (1999). *Design Theory*. Cambridge
University Press.

Colbourn, C. J., & Dinitz, J. H. (Eds.). (2007). *Handbook of
Combinatorial Designs* (2nd ed.). Chapman & Hall/CRC.

## See also

[`kirkman_matrix`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md),
[`pg_matrix`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`compare_designs`](https://myaseen208.github.io/ppgt/reference/compare_designs.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
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
M <- bibd_matrix(v = 7, k = 3)
M
#> 7 x 7 sparse Matrix of class "dgCMatrix"
#>                   
#> [1,] 1 1 . 1 . . .
#> [2,] . 1 1 . 1 . .
#> [3,] . . 1 1 . 1 .
#> [4,] . . . 1 1 . 1
#> [5,] 1 . . . 1 1 .
#> [6,] . 1 . . . 1 1
#> [7,] 1 . 1 . . . 1
Matrix::rowSums(M)
#> [1] 3 3 3 3 3 3 3
Matrix::colSums(M)
#> [1] 3 3 3 3 3 3 3

M_fail <- bibd_matrix(v = 8, k = 3)
#> Warning: No BIBD exists: lambda(v-1) not divisible by (k-1)
M_fail
#> NULL
```
