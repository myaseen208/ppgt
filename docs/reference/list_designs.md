# List Available Pooling Design Families

Return a compact catalog of the pooling design families exposed by this
package layer, together with their structural interpretations.

## Usage

``` r
list_designs()
```

## Value

A data frame with columns `type`, `function_name`, `description`,
`pools_per_sample`, and `detection_guarantee`. Each row describes one
design family and its associated constructor or parameterization.

## Details

The returned table is a descriptive index rather than a design matrix.
Its entries summarize how each family maps individuals to pools in a
binary matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\. In particular,
the `pools_per_sample` column records the intended column weight \\w_i =
\|\mathcal{J}\_i\|\\, while `detection_guarantee` records a brief
family-specific statement about the positive-pattern resolution
guaranteed by the corresponding construction.

## References

Du, D.-Z., & Hwang, F. K. (2000). *Combinatorial Group Testing and Its
Applications* (2nd ed.). World Scientific.

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

## See also

[`compare_designs`](https://myaseen208.github.io/ppgt/reference/compare_designs.md),
[`dorfman_matrix`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md),
[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)

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
[`is_kts()`](https://myaseen208.github.io/ppgt/reference/is_kts.md),
[`kirkman_matrix()`](https://myaseen208.github.io/ppgt/reference/kirkman_matrix.md),
[`pbest_clinical_matrix()`](https://myaseen208.github.io/ppgt/reference/pbest_clinical_matrix.md),
[`pg_matrix()`](https://myaseen208.github.io/ppgt/reference/pg_matrix.md),
[`regular_pooling_matrix()`](https://myaseen208.github.io/ppgt/reference/regular_pooling_matrix.md),
[`tapestry_matrix()`](https://myaseen208.github.io/ppgt/reference/tapestry_matrix.md)

## Examples

``` r
design_catalog <- list_designs()
design_catalog
#>        type                 function_name
#> 1   dorfman                dorfman_matrix
#> 2  array_2d                  array_matrix
#> 3  array_3d               array_matrix_3d
#> 4      bibd                   bibd_matrix
#> 5 hypercube              hypercube_matrix
#> 6   kirkman                kirkman_matrix
#> 7        pg                     pg_matrix
#> 8        pp                     pp_matrix
#> 9     pbest pp_matrix(q=8,d=3,nl=6,N=384)
#>                                   description pools_per_sample
#> 1 Dorfman two-stage pooling (non-overlapping)                1
#> 2                2D array (row/column) design                2
#> 3                             3D array design                3
#> 4            Balanced Incomplete Block Design                r
#> 5                     d-dimensional hypercube                d
#> 6                       Kirkman Triple System          (v-1)/2
#> 7                 Projective Geometry PG(2,q)              q+1
#> 8                  Polynomial Pools algorithm               nl
#> 9             P-BEST: 384 samples in 48 pools                6
#>   detection_guarantee
#> 1                   1
#> 2                   1
#> 3                   1
#> 4             depends
#> 5                   1
#> 6                   1
#> 7             depends
#> 8 floor((nl-1)/(d-1))
#> 9                   2
```
