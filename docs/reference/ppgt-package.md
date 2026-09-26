# ppgt: Pooled Group Testing Design Toolkit

The ppgt package provides pooling matrix constructors, decoders,
finite-field utilities, workflow-layer helpers, and comparison helpers
for pooled group testing.

## Details

The package architecture is layered: matrix construction returns
\\\mathbf{M}\\, decoding interprets \\\mathbf{z}\\ relative to a fixed
design, simulation generates \\\widetilde{\mathbf{z}}\\ and
\\\mathbf{z}\\, workflow execution coordinates stage-1 and optional
stage-2 logic, and reporting summarizes stored workflow results.
[`PoolMatrix()`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md)
belongs only to the constructor layer.

The package documentation follows the notation defined in
`Notations.tex`. Use
[`ppgt_notations`](https://myaseen208.github.io/ppgt/reference/ppgt_notations.md)
to inspect the package notation registry programmatically.

The package documentation is organized into one introductory vignette
and a set of topic-specific articles:

- `ppgt_theory`: notation, package overview, and short workflow

- `ppgt_pp`: Polynomial Pools construction and decoding

- `ppgt_pbest`: PP-based P-BEST and clinical context

- `ppgt_hyper`: HYPER matrix generation

- `ppgt_hyperec`: package-specific experimental HyperEC

- `ppgt_galois_field`: finite-field arithmetic

- `ppgt_separable_disjunct`: separable and disjunct matrices

- `ppgt_comparisons`: package comparison helpers

- `ppgt_diagnostics`: incidence-matrix overlap diagnostics based on
  \\\mathbf{M}\mathbf{M}^\top\\ and \\\mathbf{M}^\top\mathbf{M}\\

In package notation, \$\$ \mathbf{M} \in \\0,1\\^{J \times N}, \qquad
\widetilde{\mathbf{y}} = (\widetilde{y}\_1, \ldots,
\widetilde{y}\_N)^\top, \qquad \widetilde{\mathbf{z}} =
(\widetilde{z}\_1, \ldots, \widetilde{z}\_J)^\top, \qquad \mathbf{z} =
(z_1, \ldots, z_J)^\top. \$\$ The latent pool outcomes are induced by
\$\$ \widetilde{z}\_j = \mathbb{I}\left(\sum\_{i = 1}^N
M\_{ji}\widetilde{y}\_i \> 0\right), \qquad j \in \\1, \ldots, J\\. \$\$
The help pages summarize the exported interfaces, while the vignettes
carry the longer derivations, worked examples, and cross-family
comparisons.

## References

Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
information theory perspective. *Foundations and Trends in
Communications and Information Theory*, 15(3-4), 196-392.

Tan, Y. H. I. (2020). Pooling matrix designs for group testing. *SIAM
Undergraduate Research Online*, 13, 1-21.

Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing
for SARS-CoV-2: A flexible, easy-to-implement, and highly efficient
method. *Nature Communications*, 13, Article 3626.

## See also

[`PoolMatrix`](https://myaseen208.github.io/ppgt/reference/PoolMatrix.md),
[`pp_matrix`](https://myaseen208.github.io/ppgt/reference/pp_matrix.md),
[`pp_decode`](https://myaseen208.github.io/ppgt/reference/pp_decode.md),
[`hyper_matrix`](https://myaseen208.github.io/ppgt/reference/hyper_matrix.md),
[`SeparableMatrix`](https://myaseen208.github.io/ppgt/reference/SeparableMatrix.md),
[`DisjunctMatrix`](https://myaseen208.github.io/ppgt/reference/DisjunctMatrix.md),
[`simulate_group_testing`](https://myaseen208.github.io/ppgt/reference/simulate_group_testing.md),
[`run_testing_workflow`](https://myaseen208.github.io/ppgt/reference/run_testing_workflow.md),
[`protocol_summary`](https://myaseen208.github.io/ppgt/reference/protocol_summary.md),
[`ppgt_notations`](https://myaseen208.github.io/ppgt/reference/ppgt_notations.md)

## Author

**Maintainer**: Muhammad Yaseen <myaseen208@gmail.com> \[copyright
holder\]

Authors:

- Muhammad Yaseen <myaseen208@gmail.com> \[copyright holder\]

Other contributors:

- Christopher McMahan ([ORCID](https://orcid.org/0000-0001-5056-9615))
  \[contributor\]

- Christopher Bilder ([ORCID](https://orcid.org/0000-0002-2848-8576))
  \[contributor\]

- Joshua Tebbs ([ORCID](https://orcid.org/0000-0002-6762-7241))
  \[contributor\]

- Pranta Das \[contributor\]
