# Generate A Dorfman Pooling Matrix

Construct the binary pooling design matrix for classical two-stage
Dorfman testing with non-overlapping pools of common target size `g`.

## Usage

``` r
dorfman_matrix(N, g)
```

## Arguments

- N:

  Integer scalar giving the total number of individuals. This is the
  column dimension \\N\\ of the returned design matrix \\\mathbf{M} \in
  \\0,1\\^{J \times N}\\ and must satisfy \\N \ge 1\\.

- g:

  Integer scalar giving the target pool size. It must satisfy \\2 \le g
  \le N\\. In the returned design, each pool \\\mathcal{P}\_j\\ has size
  at most \\g\\, and all but possibly the last pool have size exactly
  \\g\\.

## Value

A sparse binary matrix of class `dgCMatrix` with dimensions \\J \times
N\\, where \\J = \lceil N / g \rceil\\. Entry \\M\_{ji} = 1\\ if and
only if individual \\i\\ belongs to pool \\j\\. Every column sum is one,
so \\w_i = \|\mathcal{J}\_i\| = 1\\ for every individual \\i\\. The
matrix carries `attr(., "design")` with named fields `type`, `N`, `M`,
`pool_size`, and `pools_per_sample`, where `M` is a legacy display label
for the number of pools \\J\\.

## Details

Let \\N\\ denote the total number of individuals. Dorfman testing
partitions the individuals into \$\$J = \left\lceil \frac{N}{g}
\right\rceil\$\$ disjoint pools. The returned matrix \\\mathbf{M} \in
\\0,1\\^{J \times N}\\ satisfies \$\$\mathcal{P}\_j = \\i : M\_{ji} =
1\\, \qquad \mathcal{J}\_i = \\j : M\_{ji} = 1\\,\$\$ with
\$\$\|\mathcal{J}\_i\| = 1 \quad \text{for all } i \in
\\1,\ldots,N\\.\$\$ Hence the stage-1 pools form a partition of the
individual set. If \\\widetilde{y}\_i \in \\0,1\\\\ denotes the latent
status of individual \\i\\, then the latent status of pool \\j\\ is
\$\$\widetilde{z}\_j = \mathbb{I}\left(\sum\_{i \in \mathcal{P}\_j}
\widetilde{y}\_i \> 0\right).\$\$ Under a homogeneous prevalence model
\\\Pr(\widetilde{y}\_i = 1) = p\\ and perfect testing, the expected
number of tests per individual in the two-stage procedure is
\$\$\frac{\mathbb{E}\[T\]}{N} = \frac{1}{g} + 1 - (1-p)^g.\$\$ The
heuristic choice \$\$g\_{\mathrm{opt}} \approx p^{-1/2}\$\$ motivates
[`dorfman_optimal_g`](https://myaseen208.github.io/ppgt/reference/dorfman_optimal_g.md).

## References

Dorfman, R. (1943). The detection of defective members of large
populations. *The Annals of Mathematical Statistics*, 14(4), 436-440.

Sobel, M., & Groll, P. A. (1959). Group testing to eliminate efficiently
all defectives in a binomial sample. *Bell System Technical Journal*,
38(5), 1179-1252.

## See also

[`dorfman_optimal_g`](https://myaseen208.github.io/ppgt/reference/dorfman_optimal_g.md),
[`array_matrix`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`compare_designs`](https://myaseen208.github.io/ppgt/reference/compare_designs.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
[`bibd_matrix()`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
[`compare_designs()`](https://myaseen208.github.io/ppgt/reference/compare_designs.md),
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
M <- dorfman_matrix(N = 12, g = 3)
M
#> 4 x 12 sparse Matrix of class "dgCMatrix"
#>                             
#> [1,] 1 1 1 . . . . . . . . .
#> [2,] . . . 1 1 1 . . . . . .
#> [3,] . . . . . . 1 1 1 . . .
#> [4,] . . . . . . . . . 1 1 1
Matrix::rowSums(M)
#> [1] 3 3 3 3
Matrix::colSums(M)
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1

g_opt <- dorfman_optimal_g(0.01)
g_opt
#> [1] 11
M_opt <- dorfman_matrix(N = 33, g = g_opt)
M_opt
#> 3 x 33 sparse Matrix of class "dgCMatrix"
#>                                                                       
#> [1,] 1 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . . . . . . . . . . . . .
#> [2,] . . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1 . . . . . . . . . . .
#> [3,] . . . . . . . . . . . . . . . . . . . . . . 1 1 1 1 1 1 1 1 1 1 1
```
