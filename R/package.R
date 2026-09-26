#' ppgt: Pooled Group Testing Design Toolkit
#'
#' @description
#' The \pkg{ppgt} package provides pooling matrix constructors, decoders,
#' finite-field utilities, workflow-layer helpers, and comparison helpers for
#' pooled group testing.
#'
#' @details
#' The package architecture is layered: matrix construction returns
#' \eqn{\mathbf{M}}, decoding interprets \eqn{\mathbf{z}} relative to a fixed
#' design, simulation generates \eqn{\widetilde{\mathbf{z}}} and
#' \eqn{\mathbf{z}}, workflow execution coordinates stage-1 and optional
#' stage-2 logic, and reporting summarizes stored workflow results.
#' \code{PoolMatrix()} belongs only to the constructor layer.
#'
#' The package documentation follows the notation defined in
#' \code{Notations.tex}. Use \code{\link{ppgt_notations}} to inspect the
#' package notation registry programmatically.
#'
#' The package documentation is organized into one introductory vignette and a
#' set of topic-specific articles:
#' \itemize{
#'   \item \code{ppgt_theory}: notation, package overview, and short workflow
#'   \item \code{ppgt_pp}: Polynomial Pools construction and decoding
#'   \item \code{ppgt_pbest}: PP-based P-BEST and clinical context
#'   \item \code{ppgt_hyper}: HYPER matrix generation
#'   \item \code{ppgt_hyperec}: package-specific experimental HyperEC
#'   \item \code{ppgt_galois_field}: finite-field arithmetic
#'   \item \code{ppgt_separable_disjunct}: separable and disjunct matrices
#'   \item \code{ppgt_comparisons}: package comparison helpers
#'   \item \code{ppgt_diagnostics}: incidence-matrix overlap diagnostics based on
#'   \eqn{\mathbf{M}\mathbf{M}^\top} and \eqn{\mathbf{M}^\top\mathbf{M}}
#' }
#'
#' In package notation,
#' \deqn{
#'   \mathbf{M} \in \{0,1\}^{J \times N}, \qquad
#'   \widetilde{\mathbf{y}} = (\widetilde{y}_1, \ldots, \widetilde{y}_N)^\top,
#'   \qquad
#'   \widetilde{\mathbf{z}} = (\widetilde{z}_1, \ldots, \widetilde{z}_J)^\top,
#'   \qquad
#'   \mathbf{z} = (z_1, \ldots, z_J)^\top.
#' }
#' The latent pool outcomes are induced by
#' \deqn{
#'   \widetilde{z}_j =
#'   \mathbb{I}\left(\sum_{i = 1}^N M_{ji}\widetilde{y}_i > 0\right),
#'   \qquad j \in \{1, \ldots, J\}.
#' }
#' The help pages summarize the exported interfaces, while the vignettes carry
#' the longer derivations, worked examples, and cross-family comparisons.
#'
#' @references
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' Tan, Y. H. I. (2020). Pooling matrix designs for group testing.
#' \emph{SIAM Undergraduate Research Online}, 13, 1-21.
#'
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' @seealso
#' \code{\link{PoolMatrix}}, \code{\link{pp_matrix}}, \code{\link{pp_decode}},
#' \code{\link{hyper_matrix}}, \code{\link{SeparableMatrix}},
#' \code{\link{DisjunctMatrix}}, \code{\link{simulate_group_testing}},
#' \code{\link{run_testing_workflow}}, \code{\link{protocol_summary}},
#' \code{\link{ppgt_notations}}
#'
#' @importFrom stats runif sd
#' @importFrom utils combn
"_PACKAGE"
