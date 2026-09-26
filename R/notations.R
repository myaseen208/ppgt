#' @title Inspect the Package Notation Registry
#'
#' @description
#' Return the package-wide notation registry used by \pkg{ppgt}
#' documentation. The registry records canonical symbols, family-local
#' notation, decoder-local notation, display aliases used in summaries, and
#' deprecated or forbidden aliases that should not be used as package-wide
#' notation.
#'
#' @return A base \code{data.frame} with one row per notation entry and the
#'   columns:
#'   \describe{
#'     \item{\code{symbol}}{Human-readable symbol label.}
#'     \item{\code{render}}{Render-safe LaTeX-ready representation used in the
#'       package documentation.}
#'     \item{\code{meaning}}{Meaning of the symbol in context.}
#'     \item{\code{category}}{One of \code{"global_canonical"},
#'       \code{"family_local"}, \code{"decoder_local"},
#'       \code{"display_alias"}, or \code{"forbidden_alias"}.}
#'     \item{\code{status}}{Notation status, such as \code{"canonical"},
#'       \code{"local"}, \code{"alias"}, or \code{"deprecated"}.}
#'     \item{\code{scope}}{Where the notation should be used.}
#'     \item{\code{family}}{Associated design family, layer, or
#'       package-wide scope.}
#'     \item{\code{alias_of}}{Canonical target symbol for aliases, or
#'       \code{NA_character_} when not applicable.}
#'     \item{\code{source}}{Primary package source for the entry.}
#'     \item{\code{notes}}{Short notes describing usage boundaries or
#'       deprecation status.}
#'   }
#'
#' @details
#' The package follows the notation defined in \code{Notations.tex}. The core
#' package-wide notation is:
#' \deqn{
#'   \mathbf{M} \in \{0,1\}^{J \times N}, \qquad
#'   \widetilde{\mathbf{y}} = (\widetilde{y}_1, \ldots, \widetilde{y}_N)^\top,
#'   \qquad
#'   \widetilde{\mathbf{z}} = (\widetilde{z}_1, \ldots, \widetilde{z}_J)^\top,
#'   \qquad
#'   \mathbf{z} = (z_1, \ldots, z_J)^\top.
#' }
#' The corresponding pool-status map is
#' \deqn{
#'   \widetilde{z}_j =
#'   \mathbb{I}\left(\sum_{i = 1}^N M_{ji}\widetilde{y}_i > 0\right),
#'   \qquad
#'   j \in \{1, \ldots, J\}.
#' }
#'
#' The registry makes three distinctions explicit:
#' \itemize{
#'   \item package-wide canonical notation used throughout the main
#'     documentation;
#'   \item family-local or decoder-local symbols that are valid only inside
#'     those sections;
#'   \item aliases and deprecated forms that should not replace the canonical
#'     package notation.
#' }
#'
#' In particular, design-summary fields such as \code{pool_size},
#' \code{pools_per_sample}, \code{compression}, \code{n_pools}, and
#' \code{n_samples} are recorded as display aliases rather than as primary
#' mathematical symbols. Likewise, deprecated forms such as
#' \eqn{\widetilde{\mathbf{Y}}}, scalar \eqn{M} for the number of pools, or
#' \eqn{n_i} as a package-wide column-weight symbol are tracked explicitly when
#' they would otherwise introduce ambiguity.
#'
#' @examples
#' notations <- ppgt_notations()
#' notations
#'
#' notations[notations$family == "pp/pbest", c("symbol", "meaning", "status")]
#'
#' notations[
#'   notations$category == "display_alias",
#'   c("symbol", "alias_of", "notes")
#' ]
#'
#' @references
#' Tan, Y. H. I. (2020). Pooling matrix designs for group testing.
#' \emph{SIAM Undergraduate Research Online}, 13, 1-21.
#'
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' @seealso
#' \code{\link{ppgt-package}}, \code{\link{PoolMatrix}},
#' \code{\link{pp_matrix}}, \code{\link{pp_decode}},
#' \code{\link{hyper_matrix}}, \code{\link{SeparableMatrix}},
#' \code{\link{DisjunctMatrix}}
#'
#' @export
ppgt_notations <- function() {
  registry <- .ppgt_notation_registry_data()

  registry <- registry[
    order(
      registry$category,
      registry$family,
      registry$symbol
    ),
    ,
    drop = FALSE
  ]

  rownames(registry) <- NULL
  registry
}


.ppgt_notation_registry_data <- function() {
  entry <- function(symbol,
                    meaning,
                    category,
                    status,
                    scope,
                    family,
                    alias_of = NA_character_,
                    source,
                    notes,
                    render = symbol) {
    data.frame(
      symbol = symbol,
      render = render,
      meaning = meaning,
      category = category,
      status = status,
      scope = scope,
      family = family,
      alias_of = alias_of,
      source = source,
      notes = notes,
      stringsAsFactors = FALSE
    )
  }

  rows <- list(
    entry("N", "total number of individuals", "global_canonical", "canonical",
          "package-wide", "all", source = "Notations.tex; package docs",
          notes = "Primary package-wide scalar for sample count."),
    entry("J", "total number of pools", "global_canonical", "canonical",
          "package-wide", "all", source = "Notations.tex; package docs",
          notes = "Primary package-wide scalar for pool count."),
    entry("i", "individual index", "global_canonical", "canonical",
          "package-wide", "all", source = "Notations.tex; package docs",
          notes = "Canonical index for individuals."),
    entry("j", "pool index", "global_canonical", "canonical",
          "package-wide", "all", source = "Notations.tex; package docs",
          notes = "Canonical index for pools."),
    entry("\\mathbf{M}", "pooling design matrix", "global_canonical",
          "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Primary package-wide matrix symbol."),
    entry("M_{ji}", "entry of the pooling matrix", "global_canonical",
          "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Use with \\mathbf{M}, J, and N in package-wide notation."),
    entry("\\mathcal{P}_j", "set of individuals assigned to pool j",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Use to define pool sizes and latent pool-status formulas."),
    entry("\\mathcal{J}_i", "set of pools containing individual i",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Use to define column weights and overlap structure."),
    entry("n_j", "size of pool j", "global_canonical", "canonical",
          "package-wide", "all", source = "Notations.tex; package docs",
          notes = "Canonical package-wide notation for pool size."),
    entry("w_i", "column weight for individual i",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Canonical package-wide notation for the number of pools containing individual i."),
    entry("\\widetilde{\\mathbf{y}}",
          "vector of latent infection indicators", "global_canonical",
          "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Primary package-wide notation for latent infection statuses."),
    entry("\\widetilde{y}_i",
          "latent infection indicator for individual i",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Canonical scalar component of the latent infection vector."),
    entry("\\widetilde{\\mathbf{z}}",
          "vector of latent pool statuses", "global_canonical", "canonical",
          "package-wide", "all", source = "Notations.tex; package docs",
          notes = "Primary package-wide notation for latent pool statuses."),
    entry("\\widetilde{z}_j", "latent status of pool j", "global_canonical",
          "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Canonical scalar component of the latent pool-status vector."),
    entry("\\mathbf{z}", "vector of observed pool results",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Primary package-wide notation for observed pool outcomes."),
    entry("z_j", "observed outcome of pool j", "global_canonical",
          "canonical", "package-wide", "all",
          source = "Notations.tex; package docs",
          notes = "Use lowercase z_j only in plain-text descriptions when math mode is not available."),
    entry("p_i", "infection probability for individual i",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex", notes = "Defined in Notations.tex but used less often than population-level prevalence arguments."),
    entry("S_{e_j}", "pool-level sensitivity for pool j",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex", notes = "Defined in Notations.tex; used mainly in simulation notation."),
    entry("S_{p_j}", "pool-level specificity for pool j",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex", notes = "Defined in Notations.tex; used mainly in simulation notation."),
    entry("\\mathbf{s}_e", "vector of pool-level sensitivities",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex", notes = "Defined in Notations.tex for poolwise simulation parameters."),
    entry("\\mathbf{s}_p", "vector of pool-level specificities",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex", notes = "Defined in Notations.tex for poolwise simulation parameters."),
    entry("\\mathbf{m}_j^\\top", "row j of the pooling matrix",
          "global_canonical", "canonical", "package-wide", "all",
          source = "Notations.tex", notes = "Defined in Notations.tex; used rarely outside the notation reference."),

    entry("q", "family-local field order or pools-per-sample parameter",
          "family_local", "local", "family-local only", "pp/pbest",
          source = "pp_matrix docs; ppgt_pp vignette",
          notes = "Interpret only within the selected family; it is overloaded across families."),
    entry("d", "family-local dimension or target order parameter",
          "family_local", "local", "family-local only", "pp/pbest",
          source = "pp_matrix docs; ppgt_pp vignette",
          notes = "Interpret only within the selected family."),
    entry("n_l", "PP number of layers", "family_local", "local",
          "family-local only", "pp/pbest",
          source = "pp_matrix docs; ppgt_pp vignette",
          notes = "Use only in PP / P-BEST sections; API argument name remains nl."),
    entry("k_{\\max}", "PP guaranteed-positive bound", "family_local",
          "local", "family-local only", "pp/pbest",
          source = "pp_matrix docs; ppgt_pp vignette",
          notes = "Use only in PP / P-BEST sections."),
    entry("m", "HYPER row-count parameter for the stage-1 design",
          "family_local", "local", "family-local only", "hyper",
          source = "hyper docs",
          notes = "Use only in HYPER matrix-generation sections."),
    entry("reorder", "HYPER reorder flag", "family_local", "local",
          "family-local only", "hyper", source = "hyper docs",
          notes = "Use only in HYPER matrix-generation sections."),
    entry("\\mathcal{S}", "active set in separable/disjunct definitions",
          "family_local", "local", "family-local only",
          "separable/disjunct", source = "separable/disjunct docs",
          notes = "Use only in separable/disjunct definitions."),
    entry("\\mathbf{u}_{\\mathcal{S}}",
          "Boolean union pattern induced by an active set",
          "family_local", "local", "family-local only",
          "separable/disjunct", source = "separable/disjunct docs",
          notes = "Use only in separable/disjunct definitions."),
    entry("g", "Dorfman pool size", "family_local", "local",
          "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Use only for Dorfman designs."),
    entry("r", "family-local row count, replication number, or clinical pools-per-sample parameter",
          "family_local", "local", "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Interpret only within the selected family; it is overloaded across families."),
    entry("c", "number of columns in a two-dimensional array design",
          "family_local", "local", "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Valid only in two-dimensional array sections."),
    entry("v", "family-local point count", "family_local", "local",
          "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Use only in classical-design sections such as BIBD or Kirkman."),
    entry("b", "number of blocks in a BIBD", "family_local", "local",
          "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Valid only in BIBD sections."),
    entry("k", "family-local block size or active-count parameter",
          "family_local", "local", "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Interpret only within the selected family."),
    entry("\\lambda", "BIBD pair-incidence parameter", "family_local",
          "local", "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Valid only in BIBD sections."),
    entry("\\alpha", "projective-coordinate parameter or family-local tuning parameter",
          "family_local", "local", "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Valid only in family-local sections such as projective planes."),
    entry("\\beta", "projective-coordinate parameter or family-local tuning parameter",
          "family_local", "local", "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Valid only in family-local sections such as projective planes."),
    entry("\\gamma", "projective-coordinate parameter", "family_local",
          "local", "family-local only", "classical_designs",
          source = "classical-design docs",
          notes = "Valid only in projective-plane sections."),

    entry("\\mathbf{A}", "decoder-local sensing or design matrix",
          "decoder_local", "local", "decoder-local only", "decoder",
          source = "pp_decode docs; tapestry_hyperec docs",
          notes = "Use only in decoder-local least-squares or compressed-sensing notation."),
    entry("\\mathbf{x}", "decoder-local unknown sparse signal",
          "decoder_local", "local", "decoder-local only", "decoder",
          source = "pp_decode docs; tapestry_hyperec docs",
          notes = "Do not replace package-wide latent-status notation with this symbol outside decoder-local sections."),
    entry("\\hat{\\mathbf{x}}", "decoder-local estimate of the sparse signal",
          "decoder_local", "local", "decoder-local only", "decoder",
          source = "pp_decode docs", notes = "Use only inside decoder-local estimation formulas."),
    entry("\\mathbf{y}", "decoder-local response vector",
          "decoder_local", "local", "decoder-local only", "decoder",
          source = "pp_decode docs; tapestry_hyperec docs",
          notes = "Use only inside decoder-local sections; package-wide observed outcomes remain \\mathbf{z}."),
    entry("\\tau", "GPSR regularization parameter", "decoder_local",
          "local", "decoder-local only", "decoder",
          source = "pp_decode docs", notes = "Use only inside GPSR or related decoder derivations."),
    entry("\\mathbf{g}_k", "decoder-local gradient vector at iteration k",
          "decoder_local", "local", "decoder-local only", "decoder",
          source = "pp_decode docs", notes = "Do not confuse with measurement-noise notation used elsewhere."),
    entry("F_\\tau(\\mathbf{z})", "GPSR objective function",
          "decoder_local", "local", "decoder-local only", "decoder",
          source = "pp_decode docs", notes = "Valid only inside GPSR and related decoder-local derivations."),
    entry("\\alpha", "GPSR step-size parameter", "decoder_local",
          "local", "decoder-local only", "decoder",
          source = "pp_decode docs", notes = "Use only inside GPSR or related decoder derivations."),
    entry("\\beta", "GPSR line-search shrinkage factor",
          "decoder_local", "local", "decoder-local only", "decoder",
          source = "pp_decode docs", notes = "Use only inside GPSR or related decoder derivations."),
    entry("\\sigma", "GPSR line-search parameter", "decoder_local",
          "local", "decoder-local only", "decoder",
          source = "pp_decode docs", notes = "Use only inside GPSR or related decoder derivations."),

    entry("pool_size", "display label for pool size", "display_alias",
          "alias", "summary and display outputs", "summaries",
          alias_of = "n_j", source = "pp_matrix docs; compare_* docs",
          notes = "Alias used widely in design summaries for n_j."),
    entry("pools_per_sample", "display label for column weight",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "w_i",
          source = "list_designs docs; additional design attrs",
          notes = "Alias used widely in design summaries for w_i."),
    entry("sample_coverage", "display label for column weight in PP verification output",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "w_i", source = "pp_verify docs",
          notes = "Alias used in PP verification metadata for w_i."),
    entry("n_pools", "display label for the number of pools",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "J",
          source = "pp_verify docs; separable/disjunct attrs",
          notes = "Alias used in return metadata; the canonical mathematical symbol is J."),
    entry("n_samples", "display label for the number of individuals",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "N",
          source = "pp_verify docs; separable/disjunct attrs",
          notes = "Alias used in return metadata; the canonical mathematical symbol is N."),
    entry("compression", "display label for the ratio N / J",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "N / J",
          source = "compare_* docs; pp_design docs",
          notes = "Alias used in design summaries; the canonical mathematical ratio is N / J."),
    entry("Lambda_max", "display label for the maximum column weight",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "\\lambda_{\\max}",
          source = "compare_* docs",
          notes = "Alias used in comparison outputs; the underlying quantity is the maximum column weight."),
    entry("points_per_line", "display label for projective-plane pool size",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "n_j", source = "pg docs",
          notes = "Alias used in projective-plane summaries for pool size."),
    entry("lines_per_point", "display label for projective-plane column weight",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "w_i", source = "pg docs",
          notes = "Alias used in projective-plane summaries for column weight."),
    entry("M", "legacy display alias for the number of pools",
          "display_alias", "alias", "summary and display outputs",
          "summaries", alias_of = "J", source = "notation audit; legacy text",
          notes = "Legacy output label only; do not use plain M as the package-wide symbol for the number of pools."),

    entry("n_i", "deprecated package-wide alias for column weight",
          "forbidden_alias", "deprecated", "legacy or forbidden usage",
          "legacy", alias_of = "w_i", source = "notation audit; legacy text",
          notes = "Use w_i for package-wide column weight notation."),
    entry("\\widetilde{\\mathbf{Y}}",
          "deprecated uppercase-bold latent infection vector",
          "forbidden_alias", "deprecated", "legacy or forbidden usage",
          "legacy", alias_of = "\\widetilde{\\mathbf{y}}",
          source = "notation audit; legacy text",
          notes = "Do not use uppercase-bold latent infection vectors in package docs."),
    entry("\\mathbf{y}",
          "forbidden package-wide alias for observed pool results",
          "forbidden_alias", "deprecated", "legacy or forbidden usage",
          "legacy", alias_of = "\\mathbf{z}",
          source = "notation audit; legacy text",
          notes = "Outside decoder-local derivations, observed pool outcomes should be written as \\mathbf{z}.")
  )

  do.call(rbind, rows)
}
