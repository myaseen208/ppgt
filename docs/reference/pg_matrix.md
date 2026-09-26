# Generate A Projective Plane Pooling Matrix

Construct the incidence matrix of the projective plane
\\\mathrm{PG}(2,q)\\ over the finite field \\\mathrm{GF}(q)\\, where
points act as individuals and lines act as pools.

## Usage

``` r
pg_matrix(q)
```

## Arguments

- q:

  Integer scalar giving the field order. It must be a prime power so
  that \\\mathrm{GF}(q)\\ exists. Supported values are those handled by
  [`gf`](https://myaseen208.github.io/ppgt/reference/gf.md) and
  [`is_prime_power`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md).

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\, where \\N = J = q^2 + q + 1\\. Every row sum and every column sum
equals \\q+1\\. The attached `design` attribute is a named list with
fields `type`, `q`, `N`, `M`, `points_per_line`, and `lines_per_point`.

## Details

The projective plane \\\mathrm{PG}(2,q)\\ consists of the
one-dimensional subspaces of \\\mathrm{GF}(q)^3\\. This implementation
uses canonical homogeneous representatives of the form \$\$(1, a, b),
\qquad (0, 1, c), \qquad (0, 0, 1),\$\$ where \\a,b,c \in
\mathrm{GF}(q)\\. These representatives index both the \\N = q^2 + q +
1\\ points and, by duality, the \\J = q^2 + q + 1\\ lines. A point \\P =
\[x:y:z\]\\ is incident with a line \\L = \[\alpha:\beta:\gamma\]\\
precisely when \$\$\alpha x + \beta y + \gamma z = 0 \quad \text{in }
\mathrm{GF}(q).\$\$ The returned incidence matrix \\\mathbf{M} \in
\\0,1\\^{J \times N}\\ therefore satisfies \$\$M\_{ji} = 1 \iff P_i
\text{ lies on line } j.\$\$ The standard projective-plane identities
become \$\$n_j = \|\mathcal{P}\_j\| = q + 1, \qquad w_i =
\|\mathcal{J}\_i\| = q + 1,\$\$ and every pair of distinct individuals
belongs to exactly one common pool.

## References

Hirschfeld, J. W. P. (1998). *Projective Geometries over Finite Fields*
(2nd ed.). Oxford University Press.

Colbourn, C. J., & Dinitz, J. H. (Eds.). (2007). *Handbook of
Combinatorial Designs* (2nd ed.). Chapman & Hall/CRC.

## See also

[`gf`](https://myaseen208.github.io/ppgt/reference/gf.md),
[`is_prime_power`](https://myaseen208.github.io/ppgt/reference/is_prime_power.md),
[`bibd_matrix`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md)

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
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
M2 <- pg_matrix(2)
M2
#> 7 x 7 sparse Matrix of class "dgCMatrix"
#>                   
#> [1,] . . . . 1 1 1
#> [2,] . 1 . 1 1 . .
#> [3,] . . 1 1 . . 1
#> [4,] . 1 1 . . 1 .
#> [5,] 1 1 . . . . 1
#> [6,] 1 . . 1 . 1 .
#> [7,] 1 . 1 . 1 . .
Matrix::rowSums(M2)
#> [1] 3 3 3 3 3 3 3
Matrix::colSums(M2)
#> [1] 3 3 3 3 3 3 3

M3 <- pg_matrix(3)
M3
#> 13 x 13 sparse Matrix of class "dgCMatrix"
#>                                
#>  [1,] . . . . . . . . . 1 1 1 1
#>  [2,] . . 1 . . 1 . . 1 1 . . .
#>  [3,] . 1 . . 1 . . 1 . 1 . . .
#>  [4,] . . . . . . 1 1 1 . . . 1
#>  [5,] . . 1 . 1 . 1 . . . . 1 .
#>  [6,] . 1 . . . 1 1 . . . 1 . .
#>  [7,] . . . 1 1 1 . . . . . . 1
#>  [8,] . . 1 1 . . . 1 . . 1 . .
#>  [9,] . 1 . 1 . . . . 1 . . 1 .
#> [10,] 1 1 1 . . . . . . . . . 1
#> [11,] 1 . . . . 1 . 1 . . . 1 .
#> [12,] 1 . . . 1 . . . 1 . 1 . .
#> [13,] 1 . . 1 . . 1 . . 1 . . .
```
