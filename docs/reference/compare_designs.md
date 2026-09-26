# Compare Candidate Pooling Designs

Summarize several built-in pooling design families for a target number
of individuals by reporting the achieved matrix dimensions and simple
structural metrics.

## Usage

``` r
compare_designs(N, designs = c("dorfman", "array", "hypercube", "pp"))
```

## Arguments

- N:

  Integer scalar giving the target number of individuals. This is the
  desired column dimension \\N\\ for the candidate design matrices.

- designs:

  Character vector selecting which design families to evaluate. Allowed
  values are drawn from `"dorfman"`, `"array"`, `"hypercube"`, and
  `"pp"`. The default evaluates all four.

## Value

A data frame with one row per retained design candidate. The columns are
`type`, `N`, `M`, `pool_size`, `pools_per_sample`, and `compression`.
Here `N` is the achieved number of individuals, `M` is a legacy display
label for the number of pools \\J\\, and `compression` is the ratio \\N
/ J\\. The function returns `NULL` if no candidate design is available.

## Details

For each selected design family, this function computes a candidate
pooling matrix \\\mathbf{M} \in \\0,1\\^{J \times N}\\ or a parameter
proxy and then reports the basic design dimensions. The reported
compression metric is \$\$\mathrm{compression} = \frac{N}{J},\$\$ where
\\N\\ is the achieved number of individuals and \\J\\ is the number of
pools in the candidate design. Some families, such as array and
hypercube designs, may adjust the achieved \\N\\ upward from the target
value in order to satisfy structural constraints like \\N = rc\\ or \\N
= q^d\\. The summary is therefore comparative rather than an exact
design optimization routine.

## References

Dorfman, R. (1943). The detection of defective members of large
populations. *The Annals of Mathematical Statistics*, 14(4), 436-440.

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

## See also

[`list_designs`](https://myaseen208.github.io/ppgt/reference/list_designs.md),
[`dorfman_matrix`](https://myaseen208.github.io/ppgt/reference/dorfman_matrix.md),
[`hypercube_matrix`](https://myaseen208.github.io/ppgt/reference/hypercube_matrix.md),
[`pp_design`](https://myaseen208.github.io/ppgt/reference/pp_design.md)

Other pooling_designs:
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md),
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`additional_designs`](https://myaseen208.github.io/ppgt/reference/additional_designs.md),
[`array_matrix()`](https://myaseen208.github.io/ppgt/reference/array_matrix.md),
[`array_matrix_3d()`](https://myaseen208.github.io/ppgt/reference/array_matrix_3d.md),
[`bibd_matrix()`](https://myaseen208.github.io/ppgt/reference/bibd_matrix.md),
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
comparison_100 <- compare_designs(100)
comparison_100
#>                type   N  M pool_size pools_per_sample compression
#> dorfman     dorfman 100 10        11                1       10.00
#> array      array_2d 100 20     10/10                2        5.00
#> hypercube hypercube 125 15        25                3        8.33
#> pp        pp(9,2,3)  81 27         9                3        3.00
comparison_384 <- compare_designs(384, designs = c("dorfman", "pp"))
comparison_384
#>               type   N  M pool_size pools_per_sample compression
#> dorfman    dorfman 384 35        11                1       10.97
#> pp      pp(19,2,3) 361 57        19                3        6.33
```
