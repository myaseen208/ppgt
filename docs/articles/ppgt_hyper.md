# HYPER: Matrix Generation and Interpretation

## 1 Abstract

This vignette covers
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
and
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md).
In the package, the HYPER implementation is documented as an
upstream-compatible matrix generator. It should be read as a
construction for the stage-1 design matrix \\\mathbf{M}\\, not as a
claim that the package implements the full two-stage protocol from the
HYPER paper.

## 2 Construction Summary

The package constructor returns a sparse matrix \\\mathbf{M} \in
\\0,1\\^{J \times N}\\ in which each individual belongs to exactly \\q\\
pools, so \\w_i = \|\mathcal{J}\_i\| = q\\.

The current implementation supports \\q = 1\\, \\2\\, and \\3\\, with
the same admissibility conditions documented in the exported help page.

At the design level, the HYPER generator is a balanced low-weight
incidence construction. The key exported claim is about the matrix
\\\mathbf{M}\\ itself: the package generates the stage-1 design and
preserves the upstream-compatible combinatorial factorization pattern
for the supported values of \\q\\.

## 3 A Small HYPER Example

``` r
M_hyper <- hyper_matrix(n = 12, m = 6, q = 2, reorder = FALSE)
M_hyper
  6 x 12 sparse Matrix of class "dgCMatrix"
                              
  [1,] 1 . . . 1 . . . 1 . . 1
  [2,] . 1 . 1 . . . 1 . . . 1
  [3,] . . 1 . 1 . 1 . . . 1 .
  [4,] . . 1 . . 1 . 1 . 1 . .
  [5,] . 1 . . . 1 . . 1 . 1 .
  [6,] 1 . . 1 . . 1 . . 1 . .

dim(M_hyper)
  [1]  6 12
M_hyper[1:6, 1:12]
  6 x 12 sparse Matrix of class "dgCMatrix"
                              
  [1,] 1 . . . 1 . . . 1 . . 1
  [2,] . 1 . 1 . . . 1 . . . 1
  [3,] . . 1 . 1 . 1 . . . 1 .
  [4,] . . 1 . . 1 . 1 . 1 . .
  [5,] . 1 . . . 1 . . 1 . 1 .
  [6,] 1 . . 1 . . 1 . . 1 . .
Matrix::colSums(M_hyper)
   [1] 2 2 2 2 2 2 2 2 2 2 2 2
Matrix::rowSums(M_hyper)
  [1] 4 4 4 4 4 4
```

For this example, every individual belongs to exactly `q = 2` pools, so
`w_i = |\mathcal{J}_i| = 2`.

## 4 Reordering

``` r
M_hyper_reordered <- hyper_matrix(n = 12, m = 6, q = 2, reorder = TRUE)
M_hyper_reordered
  6 x 12 sparse Matrix of class "dgCMatrix"
                              
  [1,] 1 . . . 1 . . . 1 . . 1
  [2,] . 1 . 1 . . . 1 . . . 1
  [3,] . . 1 . 1 . 1 . . . 1 .
  [4,] . . 1 . . 1 . 1 . 1 . .
  [5,] . 1 . . . 1 . . 1 . 1 .
  [6,] 1 . . 1 . . 1 . . 1 . .
```

The `reorder` argument changes only the column order of the generated
factors. It does not change the row space, the individual coverages, or
the underlying set of incidence relations.

## 5 Compatibility Wrapper

``` r
M_wrapper <- HyperDesign(n = 12, m = 6, q = 2, reorder = FALSE)
M_wrapper
  6 x 12 sparse Matrix of class "dgCMatrix"
                              
  [1,] 1 . . . 1 . . . 1 . . 1
  [2,] . 1 . 1 . . . 1 . . . 1
  [3,] . . 1 . 1 . 1 . . . 1 .
  [4,] . . 1 . . 1 . 1 . 1 . .
  [5,] . 1 . . . 1 . . 1 . 1 .
  [6,] 1 . . 1 . . 1 . . 1 . .
```

[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md)
is a package wrapper around
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md).

Both
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
and
[`HyperDesign()`](https://myaseen208.github.io/ppgt/reference/HyperDesign.md)
are matrix generators. They do not, by themselves, implement the full
HYPER workflow of stage-1 putative-positive identification followed by
stage-2 individual retesting.

In other words, this vignette is about the construction of \\\mathbf{M}
\in \\0,1\\^{J \times N}\\, not about a complete diagnostic protocol
layered above that matrix.

## 6 Other Supported Cases

The package constructor supports `q = 1`, `q = 2`, and `q = 3`.

``` r
M_hyper_q1 <- hyper_matrix(n = 5, m = 5, q = 1, reorder = FALSE)
M_hyper_q1
  5 x 5 sparse Matrix of class "dgCMatrix"
                
  [1,] 1 . . . .
  [2,] . 1 . . .
  [3,] . . 1 . .
  [4,] . . . 1 .
  [5,] . . . . 1

M_hyper_q3 <- hyper_matrix(n = 24, m = 12, q = 3, reorder = FALSE)
M_hyper_q3
  12 x 24 sparse Matrix of class "dgCMatrix"
                                                       
   [1,] 1 . . . . 1 . . 1 . . . . . 1 . . 1 . . . . . 1
   [2,] . . . 1 1 . . . . 1 . . 1 . . . . . 1 . . 1 . .
   [3,] . 1 . . . . . 1 1 . . . . 1 . . 1 . . . . . 1 .
   [4,] . . 1 . . 1 . . . . . 1 1 . . . . 1 . . 1 . . .
   [5,] . . 1 . . . 1 . . 1 . . . . . 1 1 . . . . 1 . .
   [6,] . . . 1 . . 1 . . . 1 . . 1 . . . . . 1 1 . . .
   [7,] . . . 1 . . . 1 . . 1 . . . 1 . . 1 . . . . . 1
   [8,] . 1 . . . . . 1 . . . 1 . . 1 . . . 1 . . 1 . .
   [9,] . . 1 . . 1 . . . . . 1 . . . 1 . . 1 . . . 1 .
  [10,] 1 . . . . . 1 . . 1 . . . . . 1 . . . 1 . . 1 .
  [11,] . 1 . . 1 . . . . . 1 . . 1 . . . . . 1 . . . 1
  [12,] 1 . . . 1 . . . 1 . . . 1 . . . 1 . . . 1 . . .

dim(M_hyper_q3)
  [1] 12 24
```

## 7 Scope Of The Package Claim

The package’s claim for HYPER is compatibility at the matrix-generation
level. The HYPER paper itself studies a full testing procedure with
stage-1 putative positives and stage-2 retesting. This vignette stays at
the matrix-construction level, which is what
[`hyper_matrix()`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md)
implements.

## 8 References

- Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
  for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
  method. *Nature Communications*, 13, Article 3626.
- Upstream HYPER generator logic as summarized in the package
  documentation.
