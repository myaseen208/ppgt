# Decode Experimental Parity-Augmented HYPER Results

Apply the package-specific decoder paired with
[`hyper_ec_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md).
The decoder performs a bounded syndrome search under the package's
experimental mod-2 parity model and then applies the stage-1 HYPER
conservative rule to the corrected base rows.

## Usage

``` r
hyper_ec_decode(y_obs, design, method = c("syndrome", "none"))
```

## Arguments

- y_obs:

  Numeric or integer vector of length \\J = J_0 + K\\ containing the
  observed binary pool outcomes for the combined matrix returned by
  [`hyper_ec_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md).

- design:

  Design object returned by
  [`hyper_ec_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md).
  It must contain at least `matrix_base`, `matrix_parity`, `m`,
  `m_base`, `k_parity`, and `n`.

- method:

  Character string specifying the decoder option: `"syndrome"` for
  bounded syndrome search or `"none"` to skip the search and decode from
  the uncorrected base rows.

## Value

A named list with components:

- x_decoded:

  Integer indicator vector of length \\N\\; entry \\i\\ equals one when
  individual \\i\\ remains a stage-1 candidate

- syndrome:

  Computed syndrome vector of length \\K\\

- errors_detected:

  Syndrome weight \\\sum\_\ell s\_\ell\\

- errors_corrected:

  Weight of the base-pool flip pattern returned by the bounded search

- y_corrected:

  Corrected base-pool outcome vector of length \\J_0\\

## Details

This function is not a faithful implementation of the Hong et al. (2022)
HYPER workflow. In paper HYPER, stage 1 identifies putative positives
and stage 2 individually retests them. Here, the output `x_decoded` is
the set of stage-1 candidates that remain after optional bounded
syndrome search.

Under the current implementation:

1.  split the observation vector into base and parity parts,

2.  compute the syndrome \$\$\mathbf{s} = \mathbf{y}\_{parity} \oplus
    (\mathbf{H}\mathbf{y}\_{base}),\$\$

3.  optionally search for a low-weight base-pool flip pattern, and

4.  mark individual \\i\\ as a stage-1 candidate when every base pool in
    \\\mathcal{J}\_i\\ is positive after correction.

No deterministic error-correction guarantee is implied by this
interface.

## References

Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
method. *Nature Communications*, 13, Article 3626.

The decoder documented here is a package-specific experimental extension
and is not the stage-1 plus stage-2 workflow from the paper above.

## See also

[`hyper_ec_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_ec_matrix.md),
[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`pp_decode`](https://myaseen208.github.io/ppgt/reference/pp_decode.md)

## Examples

``` r
# Build the matching experimental design
design <- hyper_ec_matrix(n = 100, k = 2, q = 2, alpha = 0.20)

# Simulate pooled outcomes for a sparse binary signal
x <- integer(100)
x[c(10, 50)] <- 1
y_obs <- as.integer((design$matrix %*% x) > 0)

# Decode to obtain stage-1 candidates
result <- hyper_ec_decode(y_obs, design)
which(result$x_decoded == 1)
#> [1] 10 50
result$errors_detected
#> [1] 0

# Flip a couple of base-pool outcomes and rerun the bounded search
y_error <- y_obs
y_error[c(5, 15)] <- 1 - y_error[c(5, 15)]
result2 <- hyper_ec_decode(y_error, design, method = "syndrome")
result2$errors_detected
#> [1] 1
result2$errors_corrected
#> [1] 1
which(result2$x_decoded == 1)
#> [1] 10 50
```
