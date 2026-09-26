# Simulate Pool Outcomes From a Pooling Design

Simulate latent and optionally noisy observed pool outcomes from a
pooling design matrix \\\mathbf{M}\\ and a latent individual-status
vector \\\widetilde{\mathbf{y}}\\. This function belongs to the
simulation layer; it does not construct \\\mathbf{M}\\, choose a
decoder, or execute a testing protocol.

## Usage

``` r
simulate_group_testing(design, Y_tilde, s_e = NULL, s_p = NULL, seed = NULL)
```

## Arguments

- design:

  A pooling design supplied either as a binary matrix \\\mathbf{M} \in
  \\0,1\\^{J \times N}\\ or as a design object containing a component
  named `matrix`. When a design object is supplied, the simulation uses
  its matrix component as the pooling design matrix.

- Y_tilde:

  Integer or logical vector of length \\N\\ giving the latent individual
  statuses \\\widetilde{\mathbf{y}} = (\widetilde{y}\_1, \ldots,
  \widetilde{y}\_N)^\top\\. Entries are interpreted as binary
  indicators.

- s_e:

  Optional numeric scalar or numeric vector of length \\J\\ specifying
  the pool-level sensitivities \\\mathbf{s}\_e = (S\_{e_1}, \ldots,
  S\_{e_J})^\top\\. If `NULL`, the simulation is noiseless on truly
  positive pools.

- s_p:

  Optional numeric scalar or numeric vector of length \\J\\ specifying
  the pool-level specificities \\\mathbf{s}\_p = (S\_{p_1}, \ldots,
  S\_{p_J})^\top\\. If `NULL`, the simulation is noiseless on truly
  negative pools.

- seed:

  Optional integer scalar used to initialize the random-number generator
  when noisy outcomes are simulated.

## Value

A named list with components:

- `matrix`:

  The pooling design matrix \\\mathbf{M}\\ used by the simulation.

- `Y_tilde`:

  The latent individual-status vector \\\widetilde{\mathbf{y}}\\.

- `z_tilde`:

  The latent pool-status vector \\\widetilde{\mathbf{z}}\\.

- `z`:

  The observed pool-outcome vector \\\mathbf{z}\\.

- `s_e`:

  The realized pool-level sensitivity vector used in the simulation, or
  `NULL`.

- `s_p`:

  The realized pool-level specificity vector used in the simulation, or
  `NULL`.

- `mode`:

  Either `"noiseless"` or `"noisy"`.

## Details

Let \\\mathbf{M} \in \\0,1\\^{J \times N}\\ denote the pooling design
matrix and let \\\widetilde{\mathbf{y}} = (\widetilde{y}\_1, \ldots,
\widetilde{y}\_N)^\top\\ denote the latent individual-status vector. The
latent pool statuses are computed as \$\$ \widetilde{z}\_j =
\mathbb{I}\left( \sum\_{i = 1}^N M\_{ji}\widetilde{y}\_i \> 0 \right),
\qquad j \in \\1, \ldots, J\\. \$\$ Equivalently, \$\$
\widetilde{\mathbf{z}} = \left(
\mathbb{I}\left((\mathbf{M}\widetilde{\mathbf{y}})\_1 \> 0\right),
\ldots, \mathbb{I}\left((\mathbf{M}\widetilde{\mathbf{y}})\_J \>
0\right) \right)^\top. \$\$

If `s_e` and `s_p` are both `NULL`, then the observed outcome vector
equals the latent pool-status vector, \\\mathbf{z} =
\widetilde{\mathbf{z}}\\. Otherwise the function samples \\z_j \mid
\widetilde{z}\_j\\ using the supplied pool-level sensitivity and
specificity parameters.

## References

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

## See also

[`PoolMatrix`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`run_testing_workflow`](https://myaseen208.github.io/ppgt/reference/run_testing_workflow.md),
[`protocol_summary`](https://myaseen208.github.io/ppgt/reference/protocol_summary.md)

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
sim$z_tilde
#>  [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
sim$z
#>  [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0

sim_noisy <- simulate_group_testing(M, Y_tilde, s_e = 0.95, s_p = 0.99, seed = 42L)
sim_noisy
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
#>  [1] 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95
#> [16] 0.95 0.95 0.95 0.95 0.95
#> 
#> $s_p
#>  [1] 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99
#> [16] 0.99 0.99 0.99 0.99 0.99
#> 
#> $mode
#> [1] "noisy"
#> 
sim_noisy$z
#>  [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
```
