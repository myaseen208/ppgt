# Decoding Algorithms for Pooled Test Results

Algorithms for decoding pooled test results: COMP (Combinatorial
Orthogonal Matching Pursuit) for zero-false-negative detection, and GPSR
(Gradient Projection for Sparse Reconstruction) for candidate ranking.

## Details

In the package-wide notation, observed pool outcomes are denoted by
\\\mathbf{z} \in \\0,1\\^J\\. The decoding algorithms in this file use a
decoder-local notation, standard in the sparse-reconstruction
literature, in which \\\mathbf{A}\\ denotes the design matrix,
\\\mathbf{x}\\ denotes the unknown sample-status vector, and
\\\mathbf{y}\\ denotes the response vector. This local notation is
restricted to the decoder layer and does not replace the package-wide
canonical symbols.

Given a binary design matrix \\\mathbf{A} \in \\0,1\\^{m \times n}\\ and
observed pooled outcomes \\\mathbf{y} \in \\0,1\\^m\\, the noiseless
decoder-local model is \$\$\mathbf{y} = \psi(\mathbf{A}\mathbf{x}),\$\$
where \\\mathbf{x} \in \\0,1\\^n\\ is the unknown sample-status vector
and \\\psi(v) = \mathbf{1}\_{v \> 0}\\ is the element-wise OR map.

COMP identifies samples as definitely negative if they appear in at
least one negative pool, and as candidates otherwise. It guarantees no
false negatives but may have false positives when more than
\\k\_{\max}\\ samples are positive.

When COMP produces too many candidates, GPSR scores them by solving the
decoder-local sparse reconstruction problem \$\$\min\_{\mathbf{x} \geq
0} \frac{1}{2}\\\mathbf{y} - \mathbf{A}\mathbf{x}\\\_2^2 +
\tau\\\mathbf{x}\\\_1.\$\$ Higher scores indicate higher likelihood of
being truly positive.

## References

Chan, C. L., et al. (2011). Non-adaptive probabilistic group testing
with noisy measurements. *Allerton Conference*.

Figueiredo, M. A. T., Nowak, R. D., & Wright, S. J. (2007). Gradient
Projection for Sparse Reconstruction. *IEEE Journal of Selected Topics
in Signal Processing*, 1(4), 586-597.

## See also

[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md)
for creating pooling matrices

Other decoding:
[`pp_comp()`](https://myaseen208.github.io/ppgt/reference/pp_comp.md),
[`pp_decode()`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
[`pp_gpsr()`](https://myaseen208.github.io/ppgt/reference/pp_gpsr.md),
[`pp_validate()`](https://myaseen208.github.io/ppgt/reference/pp_validate.md)
