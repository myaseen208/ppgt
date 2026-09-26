# Polynomial Pools: Construction, Truncation, and Decoding

## 1 Abstract

This vignette develops the Polynomial Pools (PP) construction
implemented by
[`pp_matrix()`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md).
It describes the full design with \\N = q^d\\, the package handling of
truncated designs with \\N \< q^d\\, and the public decoding helpers
[`pp_comp()`](https://myaseen208.github.io/ppgt/reference/pp_comp.md),
[`pp_gpsr()`](https://myaseen208.github.io/ppgt/reference/pp_gpsr.md),
[`pp_decode()`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
and
[`pp_validate()`](https://myaseen208.github.io/ppgt/reference/pp_validate.md).

## 2 PP Construction

Let \\\mathbf{M} \in \\0,1\\^{J \times N}\\ denote a PP matrix. In the
full untruncated construction:

- \\N = q^d\\
- \\J = n_l q\\
- \\w_i = \|\mathcal{J}\_i\| = n_l\\ for every individual \\i\\
- \\n_j = \|\mathcal{P}\_j\| = q^{d-1}\\ for every pool \\j\\

The construction indexes individuals by field-valued coordinates and
assigns them to pools by polynomial evaluation over \\GF(q)\\. In the
full design, the number of pools is \\J = n_l q\\, and the standard COMP
detection bound is \\k\_{\max} = \left\lfloor (n_l - 1) / (d - 1)
\right\rfloor\\.

More explicitly, a PP design with field order \\q\\ and dimension \\d\\
starts from the full ambient index set of size \\ N\_{\max} = q^d. \\
Individuals are represented by field-valued coordinate tuples, and each
layer corresponds to a family of polynomial evaluations over
\\\mathrm{GF}(q)\\. In the full design, \\ J = n_l q, \qquad w_i = n_l,
\qquad n_j = q^{d-1}, \qquad k\_{\max} = \left\lfloor \frac{n_l - 1}{d -
1} \right\rfloor. \\ The PP intersection bound implies that any two
distinct individuals share at most \\d - 1\\ pools.

### 2.1 Why The COMP Bound Follows

The PP guarantee used throughout the package is a direct consequence of
the intersection bound and the uniform column weight \\w_i = n_l\\ in
the full design.

Let \\i^\star\\ be an uninfected individual and let \\\mathcal{S}
\subseteq \\1, \ldots, N\\\\ be the active set of truly infected
individuals, with \\\|\mathcal{S}\| = k\\. Since any two distinct
columns of a PP matrix meet in at most \\d - 1\\ pools, the pools
containing \\i^\star\\ that are also contaminated by positives are
bounded by \$\$ \| *{i^}* {i } \_i \|

k(d - 1). \\ If \\ n_l \> k(d - 1), \\ then at least one pool containing
\$i^\star\$ remains negative, so COMP removes \$i^\star\$ from the
candidate set. Solving this inequality for the largest guaranteed \$k\$
yields \\ k\_{} = . \$\$ That is the package-wide reason the PP help
pages state exact COMP recovery for active sets up to \\k\_{\max}\\.

## 3 Building A PP Matrix

``` r
M_pp <- pp_matrix(q = 4, d = 3, nl = 5)
M_pp
  20 x 64 sparse Matrix of class "dgCMatrix"
                                                                                 
   [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
   [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
   [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
   [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
   [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
  [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
  [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
  [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
  [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
  [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
  [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
  [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
  [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
                                                             
   [1,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [2,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [3,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [4,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [5,] . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . .
   [6,] . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . .
   [7,] 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 .
   [8,] . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1
   [9,] . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1
  [10,] . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 .
  [11,] 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . .
  [12,] . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . .
  [13,] . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1
  [14,] 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 .
  [15,] . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . .
  [16,] . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
  [17,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [18,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [19,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
  [20,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1

dim(M_pp)
  [1] 20 64
attr(M_pp, "pp")
  $q
  [1] 4
  
  $d
  [1] 3
  
  $nl
  [1] 5
  
  $N
  [1] 64
  
  $M
  [1] 20
  
  $pool_size
  [1] 16
  
  $k_max
  [1] 2
```

This small example matches the standard `q = 4`, `d = 3`, `k = 2` PP
family used throughout the package documentation.

## 4 Verifying Structure

``` r
pp_summary <- pp_verify(M_pp)
pp_summary
  $n_pools
  [1] 20
  
  $n_samples
  [1] 64
  
  $pool_size
  [1] 16
  
  $sample_coverage
  [1] 5
  
  $uniform_pools
  [1] TRUE
  
  $uniform_coverage
  [1] TRUE
```

The helper checks the realized row sizes \\n_j = \|\mathcal{P}\_j\|\\
and column weights \\w_i = \|\mathcal{J}\_i\|\\.

## 5 Truncation Behavior

The constructor also supports \\N \< q^d\\. The column weights stay
fixed at \\w_i = n_l\\, while row sizes near the truncation boundary
need not stay equal.

Operationally, truncation keeps the first \\N\\ columns of the full PP
design and therefore preserves the layer structure at the individual
level. The cost is that the uniform row identity \\n_j = q^{d-1}\\ is no
longer guaranteed once some columns are removed.

``` r
M_trunc <- pp_matrix(q = 4, d = 3, nl = 5, N = 48)
M_trunc
  20 x 48 sparse Matrix of class "dgCMatrix"
                                                                                 
   [1,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
   [2,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
   [3,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
   [4,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
   [5,] 1 . . . . 1 . . . . 1 . . . . 1 . 1 . . 1 . . . . . . 1 . . 1 . . . 1 . .
   [6,] . 1 . . 1 . . . . . . 1 . . 1 . 1 . . . . 1 . . . . 1 . . . . 1 . . . 1 .
   [7,] . . 1 . . . . 1 1 . . . . 1 . . . . . 1 . . 1 . . 1 . . 1 . . . 1 . . . .
   [8,] . . . 1 . . 1 . . 1 . . 1 . . . . . 1 . . . . 1 1 . . . . 1 . . . 1 . . 1
   [9,] 1 . . . . . 1 . . . . 1 . 1 . . . . . 1 . 1 . . 1 . . . . . 1 . . 1 . . .
  [10,] . 1 . . . . . 1 . . 1 . 1 . . . . . 1 . 1 . . . . 1 . . . . . 1 1 . . . .
  [11,] . . 1 . 1 . . . . 1 . . . . . 1 . 1 . . . . . 1 . . 1 . 1 . . . . . . 1 .
  [12,] . . . 1 . 1 . . 1 . . . . . 1 . 1 . . . . . 1 . . . . 1 . 1 . . . . 1 . 1
  [13,] 1 . . . . . . 1 . 1 . . . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . . . 1 1
  [14,] . 1 . . . . 1 . 1 . . . . . . 1 . . . 1 1 . . . . . 1 . . 1 . . . . 1 . .
  [15,] . . 1 . . 1 . . . . . 1 1 . . . 1 . . . . . . 1 . 1 . . . . 1 . . 1 . . .
  [16,] . . . 1 1 . . . . . 1 . . 1 . . . 1 . . . . 1 . 1 . . . . . . 1 1 . . . .
  [17,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1
  [18,] . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .
  [19,] . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . .
  [20,] . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 .
                             
   [1,] . . . 1 . . . 1 . . .
   [2,] 1 . . . 1 . . . 1 . .
   [3,] . 1 . . . 1 . . . 1 .
   [4,] . . 1 . . . 1 . . . 1
   [5,] . . 1 1 . . . . 1 . .
   [6,] . 1 . . 1 . . 1 . . .
   [7,] 1 . . . . 1 . . . . 1
   [8,] . . . . . . 1 . . 1 .
   [9,] . . 1 . . 1 . 1 . . .
  [10,] . 1 . . . . 1 . 1 . .
  [11,] 1 . . 1 . . . . . 1 .
  [12,] . . . . 1 . . . . . 1
  [13,] . . . . . 1 . . 1 . .
  [14,] 1 . . . . . 1 1 . . .
  [15,] . 1 . 1 . . . . . . 1
  [16,] . . 1 . 1 . . . . 1 .
  [17,] . . . 1 . . . 1 . . .
  [18,] 1 . . . 1 . . . 1 . .
  [19,] . 1 . . . 1 . . . 1 .
  [20,] . . 1 . . . 1 . . . 1

dim(M_trunc)
  [1] 20 48
table(Matrix::rowSums(M_trunc))
  
  12 
  20
table(Matrix::colSums(M_trunc))
  
   5 
  48
M_trunc[1:8, 1:16]
  8 x 16 sparse Matrix of class "dgCMatrix"
                                      
  [1,] 1 . . . 1 . . . 1 . . . 1 . . .
  [2,] . 1 . . . 1 . . . 1 . . . 1 . .
  [3,] . . 1 . . . 1 . . . 1 . . . 1 .
  [4,] . . . 1 . . . 1 . . . 1 . . . 1
  [5,] 1 . . . . 1 . . . . 1 . . . . 1
  [6,] . 1 . . 1 . . . . . . 1 . . 1 .
  [7,] . . 1 . . . . 1 1 . . . . 1 . .
  [8,] . . . 1 . . 1 . . 1 . . 1 . . .
```

This is the key structural distinction between the full PP design and a
truncated one.

### 5.1 Matrix Products And PP Structure

Two Gramian-style products are especially informative for PP matrices:
\\ \mathbf{M}\mathbf{M}^\top \qquad\text{and}\qquad
\mathbf{M}^\top\mathbf{M}. \\ For the full design, their diagonals
recover the canonical row and column weights, \\
(\mathbf{M}\mathbf{M}^\top)\_{jj} = n_j = q^{d-1}, \qquad
(\mathbf{M}^\top\mathbf{M})\_{ii} = w_i = n_l, \\ while the off-diagonal
entries record pairwise overlaps. The PP intersection bound is a
statement about the off-diagonal entries of
\\\mathbf{M}^\top\mathbf{M}\\: \\ (\mathbf{M}^\top\mathbf{M})\_{ii'} =
\|\mathcal{J}\_i \cap \mathcal{J}\_{i'}\| \le d - 1 \qquad \text{for } i
\ne i'. \\ This is the matrix-product view of the same combinatorial
structure used by the COMP recovery argument.

``` r
pp_pool_overlap <- Matrix::tcrossprod(M_pp)
pp_pool_overlap
  20 x 20 sparse Matrix of class "dsCMatrix"
                                                                   
   [1,] 16  .  .  .  4  4  4  4  4  4  4  4  4  4  4  4 16  .  .  .
   [2,]  . 16  .  .  4  4  4  4  4  4  4  4  4  4  4  4  . 16  .  .
   [3,]  .  . 16  .  4  4  4  4  4  4  4  4  4  4  4  4  .  . 16  .
   [4,]  .  .  . 16  4  4  4  4  4  4  4  4  4  4  4  4  .  .  . 16
   [5,]  4  4  4  4 16  .  .  .  4  4  4  4  4  4  4  4  4  4  4  4
   [6,]  4  4  4  4  . 16  .  .  4  4  4  4  4  4  4  4  4  4  4  4
   [7,]  4  4  4  4  .  . 16  .  4  4  4  4  4  4  4  4  4  4  4  4
   [8,]  4  4  4  4  .  .  . 16  4  4  4  4  4  4  4  4  4  4  4  4
   [9,]  4  4  4  4  4  4  4  4 16  .  .  .  4  4  4  4  4  4  4  4
  [10,]  4  4  4  4  4  4  4  4  . 16  .  .  4  4  4  4  4  4  4  4
  [11,]  4  4  4  4  4  4  4  4  .  . 16  .  4  4  4  4  4  4  4  4
  [12,]  4  4  4  4  4  4  4  4  .  .  . 16  4  4  4  4  4  4  4  4
  [13,]  4  4  4  4  4  4  4  4  4  4  4  4 16  .  .  .  4  4  4  4
  [14,]  4  4  4  4  4  4  4  4  4  4  4  4  . 16  .  .  4  4  4  4
  [15,]  4  4  4  4  4  4  4  4  4  4  4  4  .  . 16  .  4  4  4  4
  [16,]  4  4  4  4  4  4  4  4  4  4  4  4  .  .  . 16  4  4  4  4
  [17,] 16  .  .  .  4  4  4  4  4  4  4  4  4  4  4  4 16  .  .  .
  [18,]  . 16  .  .  4  4  4  4  4  4  4  4  4  4  4  4  . 16  .  .
  [19,]  .  . 16  .  4  4  4  4  4  4  4  4  4  4  4  4  .  . 16  .
  [20,]  .  .  . 16  4  4  4  4  4  4  4  4  4  4  4  4  .  .  . 16

pp_sample_overlap <- Matrix::crossprod(M_pp)
pp_sample_overlap[1:12, 1:12]
  12 x 12 sparse Matrix of class "dsCMatrix"
                               
   [1,] 5 . . . 2 1 1 1 2 1 1 1
   [2,] . 5 . . 1 2 1 1 1 2 1 1
   [3,] . . 5 . 1 1 2 1 1 1 2 1
   [4,] . . . 5 1 1 1 2 1 1 1 2
   [5,] 2 1 1 1 5 . . . 2 1 1 1
   [6,] 1 2 1 1 . 5 . . 1 2 1 1
   [7,] 1 1 2 1 . . 5 . 1 1 2 1
   [8,] 1 1 1 2 . . . 5 1 1 1 2
   [9,] 2 1 1 1 2 1 1 1 5 . . .
  [10,] 1 2 1 1 1 2 1 1 . 5 . .
  [11,] 1 1 2 1 1 1 2 1 . . 5 .
  [12,] 1 1 1 2 1 1 1 2 . . . 5
```

## 6 Latent And Observed Pool Status

Given latent individual statuses \\\widetilde{\mathbf{y}} \in
\\0,1\\^N\\, the latent pool status is \\\widetilde{z}\_j =
\mathbb{I}\left(\sum\_{i = 1}^N M\_{ji}\widetilde{y}\_i \> 0\right)\\.

In the simple noiseless examples below, we set the observed outcomes
equal to the latent pool statuses.

Equivalently, \\ \widetilde{\mathbf{z}} = \left(
\mathbb{I}\left((\mathbf{M}\widetilde{\mathbf{y}})\_1 \> 0\right),
\ldots, \mathbb{I}\left((\mathbf{M}\widetilde{\mathbf{y}})\_J \>
0\right) \right)^\top, \\ and the noiseless examples use \\\mathbf{z} =
\widetilde{\mathbf{z}}\\.

``` r
Y_tilde <- rep(0L, 64)
Y_tilde[c(2, 13)] <- 1L
Y_tilde
   [1] 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
  [39] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

z <- as.integer(as.vector(M_pp %*% Y_tilde) > 0)
z
   [1] 1 1 0 0 0 1 0 1 0 1 0 0 0 1 1 0 1 1 0 0
```

## 7 COMP Decoding

``` r
comp_result <- pp_comp(M_pp, z)
comp_result
  $candidates
  [1]  2 13
  
  $negatives
   [1]  1  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27
  [26] 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52
  [51] 53 54 55 56 57 58 59 60 61 62 63 64
  
  $n_positive_pools
  [1] 9
```

For this small example,
[`pp_comp()`](https://myaseen208.github.io/ppgt/reference/pp_comp.md)
recovers the two positives exactly. The COMP rule marks an individual as
negative whenever that individual appears in at least one observed
negative pool. If \\\mathcal{N}\_0 = \\j : Z_j = 0\\\\ is the set of
negative pools, then the COMP candidate set is \\
\mathcal{C}\_{\mathrm{COMP}} = \left\\ i \in \\1, \ldots, N\\ : M\_{ji}
= 0 \text{ for all } j \in \mathcal{N}\_0 \right\\. \\

## 8 GPSR Candidate Scoring

[`pp_gpsr()`](https://myaseen208.github.io/ppgt/reference/pp_gpsr.md)
solves the nonnegative \\\ell_1\\-regularized least-squares problem used
by the package sparse-scoring layer.

In this decoder-local subsection, the symbols \\\mathbf{A}\\,
\\\mathbf{x}\\, and \\\mathbf{y}\\ follow the sparse-reconstruction
literature and are local to the decoder. They do not replace the
package-wide canonical symbols \\\mathbf{M}\\,
\\\widetilde{\mathbf{y}}\\, and \\\mathbf{z}\\. In that decoder-local
notation, the scoring problem is \\ \min\_{\mathbf{x} \ge 0}
\frac{1}{2}\\\mathbf{y} - \mathbf{A}\mathbf{x}\\\_2^2 + \tau
\\\mathbf{x}\\\_1, \\ where \\\mathbf{A}\\ is typically the
positive-pool submatrix restricted to the COMP candidate set. The
solution is interpreted as a sparse score vector over those candidates
rather than as a final binary call by itself.

``` r
A_sub <- M_pp[which(z == 1), comp_result$candidates, drop = FALSE]
y_sub <- rep(1, nrow(A_sub))
A_sub
  9 x 2 sparse Matrix of class "dgCMatrix"
           
   [1,] . 1
   [2,] 1 .
   [3,] 1 .
   [4,] . 1
   [5,] 1 1
   [6,] 1 .
   [7,] . 1
   [8,] . 1
   [9,] 1 .
y_sub
  [1] 1 1 1 1 1 1 1 1 1

gpsr_result <- pp_gpsr(A_sub, y_sub)
gpsr_result
  $x
  [1] 0.8295815 0.8295815
  
  $iterations
  [1] 5
  
  $converged
  [1] TRUE
```

## 9 Complete Package Decoding Pipeline

[`pp_decode()`](https://myaseen208.github.io/ppgt/reference/pp_decode.md)
is the package pipeline combining COMP, verification, and GPSR
refinement when needed. It is package functionality built around the PP
matrix interface; it should not be described as a blanket claim of exact
equivalence to every bundled upstream decoding path.

The implemented logic is:

1.  threshold or binarize the observed pool outcomes \\\mathbf{z}\\
2.  run COMP to eliminate obvious negatives
3.  if the candidate set is already small enough, return a direct
    stage-1 call
4.  otherwise score the candidate set with GPSR
5.  search a bounded sparse subset and report package decoding outputs
    such as `positive`, `suspected`, and `negative`

``` r
decode_result <- pp_decode(M_pp, z, verbose = FALSE)
decode_result
  $positives
  [1]  2 13
  
  $n_positives
  [1] 2
  
  $suspected
  integer(0)
  
  $n_suspected
  [1] 0
  
  $candidates
  [1]  2 13
  
  $n_candidates
  [1] 2
  
  $error
  [1] 0
  
  $method
  [1] "COMP (exact)"

decode_result$positives
  [1]  2 13
decode_result$suspected
  integer(0)
decode_result$negatives[1:12]
  NULL
```

## 10 Validation Against Known Positives

``` r
validation <- pp_validate(decode_result$positives, c(2, 13))
validation
  $TP
  [1] 2
  
  $FP
  [1] 0
  
  $FN
  [1] 0
  
  $sensitivity
  [1] 1
  
  $precision
  [1] 1
  
  $exact_match
  [1] TRUE
```

## 11 Exploring Candidate Configurations

``` r
design_grid <- pp_design(N_min = 50, N_max = 150, M_max = 40, k_min = 2, pool_size_max = 32)
design_grid
     q d nl   N  M pool_size k_max compression
  1  5 3  5 125 25        25     2        5.00
  2  5 3  6 125 30        25     2        4.17
  3 11 2  3 121 33        11     2        3.67
  4  4 3  5  64 20        16     2        3.20
  5  9 2  3  81 27         9     2        3.00
  6  8 2  3  64 24         8     2        2.67
  7  9 2  4  81 36         9     3        2.25
  8  8 2  4  64 32         8     3        2.00
  9  8 2  5  64 40         8     4        1.60
head(design_grid, 6)
     q d nl   N  M pool_size k_max compression
  1  5 3  5 125 25        25     2        5.00
  2  5 3  6 125 30        25     2        4.17
  3 11 2  3 121 33        11     2        3.67
  4  4 3  5  64 20        16     2        3.20
  5  9 2  3  81 27         9     2        3.00
  6  8 2  3  64 24         8     2        2.67
```

## 12 References

- Tan, Y. H. I. (2020). Pooling matrix designs for group testing. *SIAM
  Undergraduate Research Online*, 13, 1-21.
- PP reference algorithms as discussed in the package documentation and
  source comparisons.
