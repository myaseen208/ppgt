# Validate Decoding Results

Compares decoded results against known true positives and computes
performance metrics.

## Usage

``` r
pp_validate(decoded, true_pos)
```

## Arguments

- decoded:

  Integer vector containing the decoded positive individual indices,
  that is, a subset of \\\\1,\ldots,N\\\\ returned by a decoder.

- true_pos:

  Integer vector containing the true positive individual indices, that
  is, the ground-truth active set \\\mathcal{S}\\.

## Value

A list with components:

- TP:

  Integer; true positives

- FP:

  Integer; false positives

- FN:

  Integer; false negatives

- sensitivity:

  Numeric; TP / (TP + FN), in \\\[0,1\]\\

- precision:

  Numeric; TP / (TP + FP), in \\\[0,1\]\\

- exact_match:

  Logical; TRUE if perfect identification

## Metrics Computed

- **True Positives (TP)**: Correctly identified positives

- **False Positives (FP)**: Incorrectly identified as positive

- **False Negatives (FN)**: Missed positives

- **Sensitivity**: TP / (TP + FN) = fraction of positives found

- **Precision**: TP / (TP + FP) = fraction of identified that are truly
  positive

- **Exact Match**: TRUE if perfect identification (TP = all true, FP =
  0)

## COMP Properties

For COMP decoding with PP matrices:

- **FN is always 0** (no false negatives guaranteed)

- **Sensitivity is always 1** (100%)

- **Exact match** when \\k \leq k\_{max}\\

## References

Altman, D. G., & Bland, J. M. (1994). Diagnostic tests 1: Sensitivity
and specificity. *BMJ*, 308(6943), 1552.

Zhou, X.-H., Obuchowski, N. A., & McClish, D. K. (2011). *Statistical
Methods in Diagnostic Medicine* (2nd ed.). Wiley.

## See also

[`pp_decode`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
[`pp_comp`](https://myaseen208.github.io/ppgt/reference/pp_comp.md)

Other decoding:
[`decoding_algorithms`](https://myaseen208.github.io/ppgt/reference/decoding_algorithms.md),
[`pp_comp()`](https://myaseen208.github.io/ppgt/reference/pp_comp.md),
[`pp_decode()`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
[`pp_gpsr()`](https://myaseen208.github.io/ppgt/reference/pp_gpsr.md)

## Examples

``` r
# Perfect match
pp_validate(c(72, 142), c(72, 142))
#> $TP
#> [1] 2
#> 
#> $FP
#> [1] 0
#> 
#> $FN
#> [1] 0
#> 
#> $sensitivity
#> [1] 1
#> 
#> $precision
#> [1] 1
#> 
#> $exact_match
#> [1] TRUE
#> 

# One false positive
pp_validate(c(72, 142, 50), c(72, 142))
#> $TP
#> [1] 2
#> 
#> $FP
#> [1] 1
#> 
#> $FN
#> [1] 0
#> 
#> $sensitivity
#> [1] 1
#> 
#> $precision
#> [1] 0.6666667
#> 
#> $exact_match
#> [1] FALSE
#> 

# One false negative (unusual for COMP)
pp_validate(c(72), c(72, 142))
#> $TP
#> [1] 1
#> 
#> $FP
#> [1] 0
#> 
#> $FN
#> [1] 1
#> 
#> $sensitivity
#> [1] 0.5
#> 
#> $precision
#> [1] 1
#> 
#> $exact_match
#> [1] FALSE
#> 

# Complete example
M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
true_pos <- c(72, 142)
x <- rep(0L, 384)
x[true_pos] <- 1L
y <- as.integer(as.vector(M %*% x) > 0)

result <- pp_decode(M, y, verbose = FALSE)
metrics <- pp_validate(result$positives, true_pos)
metrics$exact_match  # TRUE
#> [1] TRUE
metrics$sensitivity  # 1.0
#> [1] 1
```
