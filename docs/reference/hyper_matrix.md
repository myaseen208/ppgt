# Generate A HYPER-Style Pooling Matrix

Construct a binary HYPER-style pooling matrix by cyclically developing
combinatorial starter blocks for \\q = 1\\, \\q = 2\\, or \\q = 3\\. The
returned design matches the package's upstream-compatible HYPER
generator and produces a matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\
with exactly \\q\\ ones in each column.

## Usage

``` r
hyper_matrix(n, m, q = 1L, reorder = TRUE)
```

## Arguments

- n:

  Integer scalar giving the total number of individuals. This is the
  column dimension \\N\\ of the returned design matrix and must satisfy
  \\N \ge 1\\.

- m:

  Integer scalar giving the total number of pools. This is the row
  dimension \\J\\ of the returned design matrix and must satisfy \\J \ge
  1\\. Additional arithmetic constraints depend on `q`.

- q:

  Integer scalar giving the number of pools containing each individual.
  Allowed values are `1L`, `2L`, and `3L`; the default is `1L`. In the
  returned matrix, \\w_i = \|\mathcal{J}\_i\| = q\\ for every individual
  \\i\\.

- reorder:

  Logical scalar controlling whether the generated column blocks are
  reordered using the upstream HYPER grouping pattern. The default is
  `TRUE`. This changes column order only; it does not change the row set
  or the column weights.

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\, where \\J = m\\ and \\N = n\\. Every column sum equals \\q\\. The
matrix carries `attr(., "design")` with named fields `type`, `n`, `m`,
`q`, and `pools_per_sample`.

## Details

Let \\\mathbf{M} \in \\0,1\\^{J \times N}\\ denote the returned pooling
matrix, with \\J = m\\ pools and \\N = n\\ individuals. The construction
assigns each individual \\i\\ to a set \\\mathcal{J}\_i\\ of exactly
\\q\\ pools, so that \$\$w_i = \|\mathcal{J}\_i\| = q \qquad \text{for
all } i \in \\1,\ldots,N\\.\$\$ The design is generated from a canonical
list of \\q\\-subsets of the pool index set, and the columns cycle
through this factor list when \\N\\ exceeds the number of available
canonical factors.

For \\q = 1\\, the factors are singleton pools and the design cycles
through the \\J\\ pools. For \\q = 2\\, the construction uses cyclic
development of starter pairs over a projective line of size \\m-1\\
augmented by an infinite point, which requires \$\$m \equiv 0 \pmod
2.\$\$ For \\q = 3\\, the construction uses developed projective orbits
and requires \$\$m \equiv 0 \pmod 6\$\$ together with the primality
condition \$\$m - 1 \text{ is prime}.\$\$ Under all three constructions,
the row sets \\\mathcal{P}\_j = \\i : M\_{ji} = 1\\\\ are intended to be
balanced as evenly as permitted by the cyclic factor development and by
the finite truncation to \\N\\ columns. When `reorder = TRUE`, the
columns are permuted by factor blocks to match the upstream HYPER
ordering convention.

## References

Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
method. *Nature Communications*, 13, Article 3626.

McMahan, C. S., Tebbs, J. M., & Bilder, C. R. (2012). Informative
Dorfman screening. *Biometrics*, 68(1), 287-296.

## See also

[`HyperDesign`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`hyper_ec_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md),
[`compare_all_designs`](https://myaseen208.github.io/ppgt/reference/compare_all_designs.md)

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
M1 <- hyper_matrix(n = 10, m = 5, q = 1)
M1
#> 5 x 10 sparse Matrix of class "dgCMatrix"
#>                         
#> [1,] 1 1 . . . . . . . .
#> [2,] . . 1 1 . . . . . .
#> [3,] . . . . 1 1 . . . .
#> [4,] . . . . . . 1 1 . .
#> [5,] . . . . . . . . 1 1
Matrix::colSums(M1)
#>  [1] 1 1 1 1 1 1 1 1 1 1

M2 <- hyper_matrix(n = 12, m = 6, q = 2)
M2
#> 6 x 12 sparse Matrix of class "dgCMatrix"
#>                             
#> [1,] 1 . . . 1 . . . 1 . . 1
#> [2,] . 1 . 1 . . . 1 . . . 1
#> [3,] . . 1 . 1 . 1 . . . 1 .
#> [4,] . . 1 . . 1 . 1 . 1 . .
#> [5,] . 1 . . . 1 . . 1 . 1 .
#> [6,] 1 . . 1 . . 1 . . 1 . .
Matrix::colSums(M2)
#>  [1] 2 2 2 2 2 2 2 2 2 2 2 2

M3 <- hyper_matrix(n = 24, m = 12, q = 3)
M3
#> 12 x 24 sparse Matrix of class "dgCMatrix"
#>                                                      
#>  [1,] 1 . . . . 1 . . 1 . . . . . 1 . . 1 . . . . . 1
#>  [2,] . . . 1 1 . . . . 1 . . 1 . . . . . 1 . . 1 . .
#>  [3,] . 1 . . . . . 1 1 . . . . 1 . . 1 . . . . . 1 .
#>  [4,] . . 1 . . 1 . . . . . 1 1 . . . . 1 . . 1 . . .
#>  [5,] . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . 1 . .
#>  [6,] . . . 1 . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
#>  [7,] . . . 1 . . . 1 . . 1 . . . 1 . . 1 . . . . . 1
#>  [8,] . 1 . . . . . 1 . . . 1 . . 1 . . . 1 . . 1 . .
#>  [9,] . . 1 . . 1 . . . . . 1 . . . 1 . . 1 . . . 1 .
#> [10,] 1 . . . . . 1 . . 1 . . . . . 1 . . . 1 . . 1 .
#> [11,] . 1 . . 1 . . . . . 1 . . 1 . . . . . 1 . . . 1
#> [12,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
Matrix::colSums(M3)
#>  [1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
```
