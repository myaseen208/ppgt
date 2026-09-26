# Score a Pooling Matrix

Summarize marginal and pairwise-overlap properties of a pooling matrix.

## Usage

``` r
score_pooling_matrix(A)
```

## Arguments

- A:

  A binary pooling matrix with pools in rows and specimens in columns.

## Value

A list containing dimensions, marginal tables, overlap summaries,
density, and an indicator for PP metadata.

## Examples

``` r
A <- regular_pooling_matrix(12, 18, 6, 4, seed = 4L)
score_pooling_matrix(A)
#> $dim
#> [1] 12 18
#> 
#> $row_sum_table
#> rs
#>  6 
#> 12 
#> 
#> $col_sum_table
#> cs
#>  4 
#> 18 
#> 
#> $max_col_overlap
#> [1] 3
#> 
#> $n_col_pairs_overlap_gt_1
#> [1] 54
#> 
#> $col_overlap_table
#> col_positive
#>   1   2   3 
#> 108  72  36 
#> 
#> $max_pool_overlap
#> [1] 6
#> 
#> $pool_overlap_table
#> pool_positive
#>  2  6 
#> 90  6 
#> 
#> $density
#> [1] 0.3333333
#> 
#> $has_pp_attr
#> [1] TRUE
#> 
```
