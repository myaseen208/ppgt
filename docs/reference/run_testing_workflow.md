# Run a Group-Testing Workflow Above the Constructor Layer

Coordinate stage-1 decoding and optional stage-2 individual retesting
above a fixed pooling design. This function belongs to the workflow
layer. It does not construct the pooling matrix \\\mathbf{M}\\; that
remains the role of
[`PoolMatrix`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
and the underlying matrix constructors.

## Usage

``` r
run_testing_workflow(
  design,
  z,
  workflow = c("non_adaptive", "adaptive_retest"),
  decoder,
  decoder_args = list(),
  stage2 = c("none", "individual_retest"),
  individual_results = NULL,
  finalize = TRUE
)
```

## Arguments

- design:

  A pooling design supplied either as a binary matrix \\\mathbf{M} \in
  \\0,1\\^{J \times N}\\ or as a design object containing a component
  named `matrix`. The function extracts \\\mathbf{M}\\ for structural
  bookkeeping, but preserves the original object for decoder calls when
  required.

- z:

  Integer, logical, or numeric vector of length \\J\\ containing the
  observed stage-1 pool outcomes \\\mathbf{z}\\.

- workflow:

  Character scalar selecting the workflow layer. The current
  implementation supports exactly two labels: `"non_adaptive"` and
  `"adaptive_retest"`.

- decoder:

  Either a function or a character scalar naming the stage-1 decoder.
  Supported character values in the current implementation are
  `"pp_decode"` and `"hyper_ec_decode"`.

- decoder_args:

  Named list of additional arguments passed to the selected decoder.

- stage2:

  Character scalar describing the stage-2 action. The current
  implementation supports exactly two values: `"none"` and
  `"individual_retest"`.

- individual_results:

  Optional vector of individual retest outcomes used only when
  `stage2 = "individual_retest"`. This may be:

  a length-\\N\\ vector

  : indexed by individual position;

  a named vector

  : whose names are retested individual indices;

  a vector of length equal to the planned retest set

  : matched in the planned retest order.

- finalize:

  Logical scalar; if `TRUE`, compute final calls from the available
  stage-1 and optional stage-2 information. If `FALSE`, return the
  workflow state without final call consolidation.

## Value

A named list with class `"ppgt_workflow"` and components:

- `workflow`:

  Workflow label.

- `matrix`:

  The pooling design matrix \\\mathbf{M}\\.

- `z`:

  Observed stage-1 pool outcomes \\\mathbf{z}\\.

- `stage1`:

  A named list describing the stage-1 decoder, its raw output, and
  standardized stage-1 calls.

- `stage2`:

  A named list describing planned retests, supplied individual results,
  and stage-2 completion status.

- `final_calls`:

  A named list with integer vectors `positive`, `negative`, and
  `unresolved`.

## Details

Let \\\mathbf{M} \in \\0,1\\^{J \times N}\\ denote the stage-1 pooling
design matrix and let \\\mathbf{z} = (z_1, \ldots, z_J)^\top\\ denote
the observed pool outcomes. This function organizes workflow execution
into two explicit layers:

1.  stage 1: decode \\\mathbf{z}\\ relative to \\\mathbf{M}\\, and

2.  stage 2: optionally retest selected individuals and consolidate
    final calls.

The function is intentionally generic and architecture-oriented. It does
not claim a single decoder or protocol is correct for every design
family. Instead, the decoder is explicit, and the stage-2 policy is
explicit.

In the current minimal implementation, the supported combinations are:

- `workflow = "non_adaptive"`, `stage2 = "none"`:

  Stage 1 only.

- `workflow = "adaptive_retest"`, `stage2 = "individual_retest"`:

  Stage-1 decoding followed by individual retesting of the stage-1
  non-negative set.

Standardized stage-1 calls are recorded as integer index sets named
`positive`, `suspected`, `negative`, and `candidates`. The consolidation
rule for the current skeleton is conservative:

- stage-1 negatives remain negative;

- stage-1 positives remain positive only in the non-adaptive
  stage-1-only case;

- when stage 2 is requested, all stage-1 non-negative individuals are
  placed into the planned retest set;

- any planned retest without an observed individual result remains
  unresolved.

This interface is stable as a layering boundary, but it is intentionally
narrow. It should be read as a workflow coordinator for the currently
implemented combinations above, not as a claim that every design family
or every protocol variant is already implemented here.

## References

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
method. *Nature Communications*, 13, Article 3626.

## See also

[`PoolMatrix`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`simulate_group_testing`](https://myaseen208.github.io/ppgt/reference/simulate_group_testing.md),
[`protocol_summary`](https://myaseen208.github.io/ppgt/reference/protocol_summary.md),
[`pp_decode`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
[`hyper_ec_decode`](https://myaseen208.github.io/ppgt/reference/hyper_ec_decode.md)

## Examples

``` r
M <- pp_matrix(q = 4, d = 3, nl = 5)
M
#> 20 x 64 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
#>  [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
#>  [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
#>  [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
#>  [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
#> [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
#> [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
#> [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
#> [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
#> [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
#> [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
#> [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
#> [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#> [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>                                                            
#>  [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
#>  [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
#>  [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
#>  [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
#>  [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
#> [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
#> [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
#> [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
#> [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
#> [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
#> [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
#> [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
#> [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#> [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1

Y_tilde <- rep(0L, 64L)
Y_tilde[c(2L, 13L)] <- 1L
Y_tilde
#>  [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

sim <- simulate_group_testing(M, Y_tilde)
sim
#> $matrix
#> 20 x 64 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
#>  [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
#>  [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
#>  [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
#>  [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
#> [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
#> [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
#> [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
#> [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
#> [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
#> [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
#> [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
#> [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#> [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>                                                            
#>  [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
#>  [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
#>  [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
#>  [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
#>  [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
#> [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
#> [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
#> [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
#> [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
#> [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
#> [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
#> [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
#> [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#> [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#> 
#> $Y_tilde
#>  [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $z_tilde
#>  [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
#> 
#> $z
#>  [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
#> 
#> $s_e
#> NULL
#> 
#> $s_p
#> NULL
#> 
#> $mode
#> [1] "noiseless"
#> 

wf1 <- run_testing_workflow(
  design = M,
  z = sim$z,
  workflow = "non_adaptive",
  decoder = "pp_decode",
  decoder_args = list(verbose = FALSE),
  stage2 = "none"
)
wf1
#> $workflow
#> [1] "non_adaptive"
#> 
#> $matrix
#> 20 x 64 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
#>  [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
#>  [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
#>  [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
#>  [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
#> [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
#> [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
#> [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
#> [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
#> [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
#> [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
#> [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
#> [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#> [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>                                                            
#>  [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
#>  [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
#>  [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
#>  [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
#>  [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
#> [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
#> [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
#> [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
#> [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
#> [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
#> [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
#> [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
#> [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#> [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#> 
#> $z
#>  [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
#> 
#> $stage1
#> $stage1$decoder
#> [1] "pp_decode"
#> 
#> $stage1$raw
#> $stage1$raw$positives
#> [1]  2 13
#> 
#> $stage1$raw$n_positives
#> [1] 2
#> 
#> $stage1$raw$suspected
#> integer(0)
#> 
#> $stage1$raw$n_suspected
#> [1] 0
#> 
#> $stage1$raw$candidates
#> [1]  2 13
#> 
#> $stage1$raw$n_candidates
#> [1] 2
#> 
#> $stage1$raw$error
#> [1] 0
#> 
#> $stage1$raw$method
#> [1] "COMP (exact)"
#> 
#> 
#> $stage1$positive
#> [1]  2 13
#> 
#> $stage1$suspected
#> integer(0)
#> 
#> $stage1$negative
#> integer(0)
#> 
#> $stage1$candidates
#> [1]  2 13
#> 
#> 
#> $stage2
#> $stage2$method
#> [1] "none"
#> 
#> $stage2$planned
#> integer(0)
#> 
#> $stage2$results
#>  [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 
#> $stage2$completed
#> [1] TRUE
#> 
#> 
#> $final_calls
#> $final_calls$positive
#> [1]  2 13
#> 
#> $final_calls$negative
#> integer(0)
#> 
#> $final_calls$unresolved
#>  [1]  1  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
#> [26] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
#> [51] 53 54 55 56 57 58 59 60 61 62 63 64
#> 
#> 
#> attr(,"class")
#> [1] "ppgt_workflow"
wf1$final_calls
#> $positive
#> [1]  2 13
#> 
#> $negative
#> integer(0)
#> 
#> $unresolved
#>  [1]  1  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
#> [26] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
#> [51] 53 54 55 56 57 58 59 60 61 62 63 64
#> 

retest_results <- integer(ncol(M))
retest_results[c(2L, 13L)] <- 1L
retest_results
#>  [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

wf2 <- run_testing_workflow(
  design = M,
  z = sim$z,
  workflow = "adaptive_retest",
  decoder = "pp_decode",
  decoder_args = list(verbose = FALSE),
  stage2 = "individual_retest",
  individual_results = retest_results
)
wf2
#> $workflow
#> [1] "adaptive_retest"
#> 
#> $matrix
#> 20 x 64 sparse Matrix of class "dgCMatrix"
#>                                                                                
#>  [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
#>  [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
#>  [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
#>  [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
#>  [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
#> [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
#> [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
#> [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
#> [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
#> [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
#> [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
#> [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
#> [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#> [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>                                                            
#>  [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#>  [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#>  [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#>  [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#>  [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
#>  [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
#>  [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
#>  [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
#>  [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
#> [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
#> [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
#> [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
#> [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
#> [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
#> [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
#> [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
#> [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
#> [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
#> [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
#> [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
#> 
#> $z
#>  [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
#> 
#> $stage1
#> $stage1$decoder
#> [1] "pp_decode"
#> 
#> $stage1$raw
#> $stage1$raw$positives
#> [1]  2 13
#> 
#> $stage1$raw$n_positives
#> [1] 2
#> 
#> $stage1$raw$suspected
#> integer(0)
#> 
#> $stage1$raw$n_suspected
#> [1] 0
#> 
#> $stage1$raw$candidates
#> [1]  2 13
#> 
#> $stage1$raw$n_candidates
#> [1] 2
#> 
#> $stage1$raw$error
#> [1] 0
#> 
#> $stage1$raw$method
#> [1] "COMP (exact)"
#> 
#> 
#> $stage1$positive
#> [1]  2 13
#> 
#> $stage1$suspected
#> integer(0)
#> 
#> $stage1$negative
#> integer(0)
#> 
#> $stage1$candidates
#> [1]  2 13
#> 
#> 
#> $stage2
#> $stage2$method
#> [1] "individual_retest"
#> 
#> $stage2$planned
#> [1]  2 13
#> 
#> $stage2$results
#>  [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $stage2$completed
#> [1] TRUE
#> 
#> 
#> $final_calls
#> $final_calls$positive
#> [1]  2 13
#> 
#> $final_calls$negative
#> integer(0)
#> 
#> $final_calls$unresolved
#>  [1]  1  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
#> [26] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
#> [51] 53 54 55 56 57 58 59 60 61 62 63 64
#> 
#> 
#> attr(,"class")
#> [1] "ppgt_workflow"
wf2$stage2
#> $method
#> [1] "individual_retest"
#> 
#> $planned
#> [1]  2 13
#> 
#> $results
#>  [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $completed
#> [1] TRUE
#> 
wf2$final_calls
#> $positive
#> [1]  2 13
#> 
#> $negative
#> integer(0)
#> 
#> $unresolved
#>  [1]  1  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
#> [26] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
#> [51] 53 54 55 56 57 58 59 60 61 62 63 64
#> 
```
