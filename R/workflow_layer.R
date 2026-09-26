#' @title Simulate Pool Outcomes From a Pooling Design
#'
#' @description
#' Simulate latent and optionally noisy observed pool outcomes from a pooling
#' design matrix \eqn{\mathbf{M}} and a latent individual-status vector
#' \eqn{\widetilde{\mathbf{y}}}. This function belongs to the simulation layer; it
#' does not construct \eqn{\mathbf{M}}, choose a decoder, or execute a testing
#' protocol.
#'
#' @param design A pooling design supplied either as a binary matrix
#'   \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} or as a design object containing
#'   a component named \code{matrix}. When a design object is supplied, the
#'   simulation uses its matrix component as the pooling design matrix.
#' @param Y_tilde Integer or logical vector of length \eqn{N} giving the latent
#'   individual statuses
#'   \eqn{\widetilde{\mathbf{y}} = (\widetilde{y}_1, \ldots, \widetilde{y}_N)^\top}.
#'   Entries are interpreted as binary indicators.
#' @param s_e Optional numeric scalar or numeric vector of length \eqn{J}
#'   specifying the pool-level sensitivities
#'   \eqn{\mathbf{s}_e = (S_{e_1}, \ldots, S_{e_J})^\top}. If \code{NULL}, the
#'   simulation is noiseless on truly positive pools.
#' @param s_p Optional numeric scalar or numeric vector of length \eqn{J}
#'   specifying the pool-level specificities
#'   \eqn{\mathbf{s}_p = (S_{p_1}, \ldots, S_{p_J})^\top}. If \code{NULL}, the
#'   simulation is noiseless on truly negative pools.
#' @param seed Optional integer scalar used to initialize the random-number
#'   generator when noisy outcomes are simulated.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{matrix}}{The pooling design matrix \eqn{\mathbf{M}} used by the
#'     simulation.}
#'   \item{\code{Y_tilde}}{The latent individual-status vector
#'     \eqn{\widetilde{\mathbf{y}}}.}
#'   \item{\code{z_tilde}}{The latent pool-status vector
#'     \eqn{\widetilde{\mathbf{z}}}.}
#'   \item{\code{z}}{The observed pool-outcome vector \eqn{\mathbf{z}}.}
#'   \item{\code{s_e}}{The realized pool-level sensitivity vector used in the
#'     simulation, or \code{NULL}.}
#'   \item{\code{s_p}}{The realized pool-level specificity vector used in the
#'     simulation, or \code{NULL}.}
#'   \item{\code{mode}}{Either \code{"noiseless"} or \code{"noisy"}.}
#' }
#'
#' @details
#' Let \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denote the pooling design
#' matrix and let
#' \eqn{\widetilde{\mathbf{y}} = (\widetilde{y}_1, \ldots, \widetilde{y}_N)^\top}
#' denote the latent individual-status vector. The latent pool statuses are
#' computed as
#' \deqn{
#'   \widetilde{z}_j
#'   =
#'   \mathbb{I}\left(
#'     \sum_{i = 1}^N M_{ji}\widetilde{y}_i > 0
#'   \right),
#'   \qquad
#'   j \in \{1, \ldots, J\}.
#' }
#' Equivalently,
#' \deqn{
#'   \widetilde{\mathbf{z}}
#'   =
#'   \left(
#'     \mathbb{I}\left((\mathbf{M}\widetilde{\mathbf{y}})_1 > 0\right),
#'     \ldots,
#'     \mathbb{I}\left((\mathbf{M}\widetilde{\mathbf{y}})_J > 0\right)
#'   \right)^\top.
#' }
#'
#' If \code{s_e} and \code{s_p} are both \code{NULL}, then the observed outcome
#' vector equals the latent pool-status vector, \eqn{\mathbf{z} =
#' \widetilde{\mathbf{z}}}. Otherwise the function samples
#' \eqn{z_j \mid \widetilde{z}_j} using the supplied pool-level sensitivity and
#' specificity parameters.
#'
#' @examples
#' M <- pp_matrix(q = 4, d = 3, nl = 5)
#' M
#'
#' Y_tilde <- rep(0L, 64L)
#' Y_tilde[c(2L, 13L)] <- 1L
#' Y_tilde
#'
#' sim <- simulate_group_testing(M, Y_tilde)
#' sim
#' sim$z_tilde
#' sim$z
#'
#' sim_noisy <- simulate_group_testing(M, Y_tilde, s_e = 0.95, s_p = 0.99, seed = 42L)
#' sim_noisy
#' sim_noisy$z
#'
#' @references
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' @seealso
#' \code{\link{PoolMatrix}}, \code{\link{run_testing_workflow}},
#' \code{\link{protocol_summary}}
#'
#' @export
simulate_group_testing <- function(design,
                                   Y_tilde,
                                   s_e = NULL,
                                   s_p = NULL,
                                   seed = NULL) {
  M <- .workflow_matrix(design)
  Y_tilde <- as.integer(as.vector(Y_tilde) > 0L)

  if (length(Y_tilde) != ncol(M)) {
    stop(
      "Length of `Y_tilde` (", length(Y_tilde),
      ") must equal ncol(M) (", ncol(M), ")."
    )
  }

  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  z_tilde <- as.integer(as.vector(M %*% Y_tilde) > 0L)

  if (is.null(s_e) && is.null(s_p)) {
    z <- z_tilde
    mode <- "noiseless"
  } else {
    s_e_vec <- .workflow_prob_vec(
      x = if (is.null(s_e)) 1 else s_e,
      J = nrow(M),
      name = "s_e"
    )
    s_p_vec <- .workflow_prob_vec(
      x = if (is.null(s_p)) 1 else s_p,
      J = nrow(M),
      name = "s_p"
    )

    z <- integer(nrow(M))
    pos_idx <- which(z_tilde == 1L)
    neg_idx <- which(z_tilde == 0L)

    if (length(pos_idx) > 0L) {
      z[pos_idx] <- as.integer(stats::runif(length(pos_idx)) < s_e_vec[pos_idx])
    }
    if (length(neg_idx) > 0L) {
      z[neg_idx] <- as.integer(stats::runif(length(neg_idx)) >= s_p_vec[neg_idx])
    }

    s_e <- s_e_vec
    s_p <- s_p_vec
    mode <- "noisy"
  }

  list(
    matrix = M,
    Y_tilde = Y_tilde,
    z_tilde = z_tilde,
    z = z,
    s_e = s_e,
    s_p = s_p,
    mode = mode
  )
}


#' @title Run a Group-Testing Workflow Above the Constructor Layer
#'
#' @description
#' Coordinate stage-1 decoding and optional stage-2 individual retesting above
#' a fixed pooling design. This function belongs to the workflow layer. It does
#' not construct the pooling matrix \eqn{\mathbf{M}}; that remains the role of
#' \code{\link{PoolMatrix}} and the underlying matrix constructors.
#'
#' @param design A pooling design supplied either as a binary matrix
#'   \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} or as a design object containing
#'   a component named \code{matrix}. The function extracts \eqn{\mathbf{M}} for
#'   structural bookkeeping, but preserves the original object for decoder calls
#'   when required.
#' @param z Integer, logical, or numeric vector of length \eqn{J} containing the
#'   observed stage-1 pool outcomes \eqn{\mathbf{z}}.
#' @param workflow Character scalar selecting the workflow layer. The current
#'   implementation supports exactly two labels:
#'   \code{"non_adaptive"} and \code{"adaptive_retest"}.
#' @param decoder Either a function or a character scalar naming the stage-1
#'   decoder. Supported character values in the current implementation are
#'   \code{"pp_decode"} and \code{"hyper_ec_decode"}.
#' @param decoder_args Named list of additional arguments passed to the selected
#'   decoder.
#' @param stage2 Character scalar describing the stage-2 action. The current
#'   implementation supports exactly two values: \code{"none"} and
#'   \code{"individual_retest"}.
#' @param individual_results Optional vector of individual retest outcomes used
#'   only when \code{stage2 = "individual_retest"}. This may be:
#'   \describe{
#'     \item{a length-\eqn{N} vector}{indexed by individual position;}
#'     \item{a named vector}{whose names are retested individual indices;}
#'     \item{a vector of length equal to the planned retest set}{matched in the
#'       planned retest order.}
#'   }
#' @param finalize Logical scalar; if \code{TRUE}, compute final calls from the
#'   available stage-1 and optional stage-2 information. If \code{FALSE}, return
#'   the workflow state without final call consolidation.
#'
#' @return A named list with class \code{"ppgt_workflow"} and components:
#' \describe{
#'   \item{\code{workflow}}{Workflow label.}
#'   \item{\code{matrix}}{The pooling design matrix \eqn{\mathbf{M}}.}
#'   \item{\code{z}}{Observed stage-1 pool outcomes \eqn{\mathbf{z}}.}
#'   \item{\code{stage1}}{A named list describing the stage-1 decoder, its raw
#'     output, and standardized stage-1 calls.}
#'   \item{\code{stage2}}{A named list describing planned retests, supplied
#'     individual results, and stage-2 completion status.}
#'   \item{\code{final_calls}}{A named list with integer vectors
#'     \code{positive}, \code{negative}, and \code{unresolved}.}
#' }
#'
#' @details
#' Let \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denote the stage-1 pooling
#' design matrix and let \eqn{\mathbf{z} = (z_1, \ldots, z_J)^\top} denote the
#' observed pool outcomes. This function organizes workflow execution into two
#' explicit layers:
#'
#' 1. stage 1: decode \eqn{\mathbf{z}} relative to \eqn{\mathbf{M}}, and
#' 2. stage 2: optionally retest selected individuals and consolidate final
#'    calls.
#'
#' The function is intentionally generic and architecture-oriented. It does not
#' claim a single decoder or protocol is correct for every design family.
#' Instead, the decoder is explicit, and the stage-2 policy is explicit.
#'
#' In the current minimal implementation, the supported combinations are:
#' \describe{
#'   \item{\code{workflow = "non_adaptive"}, \code{stage2 = "none"}}{Stage 1
#'   only.}
#'   \item{\code{workflow = "adaptive_retest"},
#'   \code{stage2 = "individual_retest"}}{Stage-1 decoding followed by
#'   individual retesting of the stage-1 non-negative set.}
#' }
#'
#' Standardized stage-1 calls are recorded as integer index sets named
#' \code{positive}, \code{suspected}, \code{negative}, and \code{candidates}.
#' The consolidation rule for the current skeleton is conservative:
#'
#' - stage-1 negatives remain negative;
#' - stage-1 positives remain positive only in the non-adaptive stage-1-only
#'   case;
#' - when stage 2 is requested, all stage-1 non-negative individuals are placed
#'   into the planned retest set;
#' - any planned retest without an observed individual result remains
#'   unresolved.
#'
#' This interface is stable as a layering boundary, but it is intentionally
#' narrow. It should be read as a workflow coordinator for the currently
#' implemented combinations above, not as a claim that every design family or
#' every protocol variant is already implemented here.
#'
#' @examples
#' M <- pp_matrix(q = 4, d = 3, nl = 5)
#' M
#'
#' Y_tilde <- rep(0L, 64L)
#' Y_tilde[c(2L, 13L)] <- 1L
#' Y_tilde
#'
#' sim <- simulate_group_testing(M, Y_tilde)
#' sim
#'
#' wf1 <- run_testing_workflow(
#'   design = M,
#'   z = sim$z,
#'   workflow = "non_adaptive",
#'   decoder = "pp_decode",
#'   decoder_args = list(verbose = FALSE),
#'   stage2 = "none"
#' )
#' wf1
#' wf1$final_calls
#'
#' retest_results <- integer(ncol(M))
#' retest_results[c(2L, 13L)] <- 1L
#' retest_results
#'
#' wf2 <- run_testing_workflow(
#'   design = M,
#'   z = sim$z,
#'   workflow = "adaptive_retest",
#'   decoder = "pp_decode",
#'   decoder_args = list(verbose = FALSE),
#'   stage2 = "individual_retest",
#'   individual_results = retest_results
#' )
#' wf2
#' wf2$stage2
#' wf2$final_calls
#'
#' @references
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' @seealso
#' \code{\link{PoolMatrix}}, \code{\link{simulate_group_testing}},
#' \code{\link{protocol_summary}}, \code{\link{pp_decode}},
#' \code{\link{hyper_ec_decode}}
#'
#' @export
run_testing_workflow <- function(design,
                                 z,
                                 workflow = c("non_adaptive", "adaptive_retest"),
                                 decoder,
                                 decoder_args = list(),
                                 stage2 = c("none", "individual_retest"),
                                 individual_results = NULL,
                                 finalize = TRUE) {
  M <- .workflow_matrix(design)
  workflow <- match.arg(workflow)
  stage2 <- match.arg(stage2)

  z <- as.integer(as.vector(z) > 0L)
  if (length(z) != nrow(M)) {
    stop("Length of `z` must equal nrow(M).")
  }
  if (!is.list(decoder_args)) {
    stop("`decoder_args` must be a named list.")
  }
  if (workflow == "non_adaptive" && stage2 != "none") {
    stop("`workflow = \"non_adaptive\"` requires `stage2 = \"none\"`.")
  }
  if (workflow == "adaptive_retest" && stage2 != "individual_retest") {
    stop(
      "`workflow = \"adaptive_retest\"` requires ",
      "`stage2 = \"individual_retest\"`."
    )
  }
  if (stage2 == "none" && !is.null(individual_results)) {
    stop("`individual_results` may be supplied only when `stage2 = \"individual_retest\"`.")
  }

  decoder_name <- .workflow_decoder_name(decoder)
  stage1_raw <- .workflow_run_decoder(
    decoder = decoder,
    design = design,
    M = M,
    z = z,
    decoder_args = decoder_args
  )
  stage1_std <- .workflow_standardize_stage1(stage1_raw, n = ncol(M))

  planned_stage2 <- integer(0L)
  if (workflow == "adaptive_retest" && stage2 == "individual_retest") {
    planned_stage2 <- sort(unique(c(
      stage1_std$positive,
      stage1_std$suspected,
      stage1_std$candidates
    )))
  }

  stage2_results <- .workflow_standardize_individual_results(
    individual_results = individual_results,
    planned = planned_stage2,
    N = ncol(M)
  )

  final_calls <- .workflow_finalize_calls(
    stage1 = stage1_std,
    planned_stage2 = planned_stage2,
    stage2_results = stage2_results,
    N = ncol(M),
    workflow = workflow,
    stage2 = stage2,
    finalize = isTRUE(finalize)
  )

  out <- list(
    workflow = workflow,
    matrix = M,
    z = z,
    stage1 = list(
      decoder = decoder_name,
      raw = stage1_raw,
      positive = stage1_std$positive,
      suspected = stage1_std$suspected,
      negative = stage1_std$negative,
      candidates = stage1_std$candidates
    ),
    stage2 = list(
      method = stage2,
      planned = planned_stage2,
      results = stage2_results,
      completed = length(planned_stage2) == 0L ||
        all(!is.na(stage2_results[planned_stage2]))
    ),
    final_calls = final_calls
  )

  class(out) <- "ppgt_workflow"
  out
}


#' @title Summarize a Workflow-Layer Result
#'
#' @description
#' Summarize the output of \code{\link{run_testing_workflow}} into a compact
#' reporting object. This function belongs to the reporting layer; it does not
#' construct matrices, simulate outcomes, or run decoders.
#'
#' @param x A workflow result object returned by \code{\link{run_testing_workflow}}.
#'
#' @return A named list containing compact workflow counts, including the number
#'   of pools, the number of individuals, the number of stage-1 non-negative
#'   calls, the number of planned stage-2 retests, and the counts of final
#'   positive, negative, and unresolved individuals.
#'
#' @details
#' Let \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denote the stage-1 pooling
#' design matrix and let \eqn{\mathbf{z}} denote the stage-1 pool outcomes. A
#' workflow result records stage-1 calls, optional stage-2 retests, and final
#' consolidated calls. The summary returned by \code{protocol_summary()} is a
#' reporting view of those stored objects. It does not re-run any protocol
#' logic.
#'
#' @examples
#' M <- pp_matrix(q = 4, d = 3, nl = 5)
#' M
#'
#' Y_tilde <- rep(0L, 64L)
#' Y_tilde[c(2L, 13L)] <- 1L
#' Y_tilde
#'
#' sim <- simulate_group_testing(M, Y_tilde)
#' sim
#'
#' wf <- run_testing_workflow(
#'   design = M,
#'   z = sim$z,
#'   workflow = "non_adaptive",
#'   decoder = "pp_decode",
#'   decoder_args = list(verbose = FALSE)
#' )
#' wf
#'
#' summary_obj <- protocol_summary(wf)
#' summary_obj
#'
#' @seealso
#' \code{\link{run_testing_workflow}}, \code{\link{simulate_group_testing}},
#' \code{\link{PoolMatrix}}
#'
#' @export
protocol_summary <- function(x) {
  if (!inherits(x, "ppgt_workflow")) {
    stop("`x` must inherit from \"ppgt_workflow\".")
  }

  list(
    workflow = x$workflow,
    n_pools = nrow(x$matrix),
    n_individuals = ncol(x$matrix),
    n_stage1_positive = length(x$stage1$positive),
    n_stage1_suspected = length(x$stage1$suspected),
    n_stage1_negative = length(x$stage1$negative),
    n_stage2_planned = length(x$stage2$planned),
    n_final_positive = length(x$final_calls$positive),
    n_final_negative = length(x$final_calls$negative),
    n_final_unresolved = length(x$final_calls$unresolved)
  )
}


# Internal helpers ------------------------------------------------------------

.workflow_matrix <- function(design) {
  if (inherits(design, "Matrix") || is.matrix(design)) {
    return(design)
  }
  if (is.list(design) && "matrix" %in% names(design)) {
    return(design$matrix)
  }
  stop("`design` must be a matrix or a design object with a `matrix` component.")
}


.workflow_prob_vec <- function(x, J, name) {
  if (length(x) == 1L) {
    x <- rep(as.numeric(x), J)
  } else {
    x <- as.numeric(x)
  }
  if (length(x) != J) {
    stop("`", name, "` must have length 1 or length nrow(M).")
  }
  if (any(is.na(x)) || any(x < 0) || any(x > 1)) {
    stop("`", name, "` must lie in [0, 1].")
  }
  x
}


.workflow_decoder_name <- function(decoder) {
  if (is.function(decoder)) {
    return("custom")
  }
  if (!is.character(decoder) || length(decoder) != 1L) {
    stop("`decoder` must be a function or a single character string.")
  }
  decoder
}


.workflow_run_decoder <- function(decoder, design, M, z, decoder_args) {
  if (is.function(decoder)) {
    args <- c(list(design = design, M = M, z = z), decoder_args)
    return(do.call(decoder, args))
  }

  if (identical(decoder, "pp_decode")) {
    args <- c(list(M = M, y = z), decoder_args)
    if (!("verbose" %in% names(args))) {
      args$verbose <- FALSE
    }
    return(do.call(pp_decode, args))
  }

  if (identical(decoder, "hyper_ec_decode")) {
    args <- c(list(y_obs = z, design = design), decoder_args)
    return(do.call(hyper_ec_decode, args))
  }

  stop("Unsupported decoder: ", decoder)
}


.workflow_standardize_stage1 <- function(x, n) {
  if (!is.list(x)) {
    stop("Decoder output must be a list.")
  }

  positive <- integer(0L)
  suspected <- integer(0L)
  negative <- integer(0L)
  candidates <- integer(0L)

  if ("positives" %in% names(x)) {
    positive <- as.integer(x$positives)
  }
  if ("suspected" %in% names(x)) {
    suspected <- as.integer(x$suspected)
  }
  if ("negative" %in% names(x)) {
    negative <- as.integer(x$negative)
  }
  if ("negatives" %in% names(x)) {
    negative <- as.integer(x$negatives)
  }
  if ("candidates" %in% names(x)) {
    candidates <- as.integer(x$candidates)
  }
  if ("x_decoded" %in% names(x)) {
    candidates <- which(as.integer(x$x_decoded) == 1L)
    suspected <- sort(unique(c(suspected, candidates)))
  }

  if (length(candidates) == 0L) {
    candidates <- sort(unique(c(positive, suspected)))
  }
  if (length(suspected) == 0L && length(candidates) > 0L && length(positive) == 0L) {
    suspected <- candidates
  }

  positive <- sort(unique(positive[positive >= 1L & positive <= n]))
  suspected <- sort(unique(suspected[suspected >= 1L & suspected <= n]))
  negative <- sort(unique(negative[negative >= 1L & negative <= n]))
  candidates <- sort(unique(candidates[candidates >= 1L & candidates <= n]))

  listed <- sort(unique(c(positive, suspected, negative, candidates)))
  unresolved <- setdiff(seq_len(n), sort(unique(c(positive, negative))))
  if (length(listed) == 0L) {
    stop("Decoder output could not be standardized into stage-1 calls.")
  }

  list(
    positive = positive,
    suspected = suspected,
    negative = negative,
    candidates = candidates,
    unresolved = unresolved
  )
}


.workflow_standardize_individual_results <- function(individual_results, planned, N) {
  out <- rep(NA_integer_, N)

  if (is.null(individual_results)) {
    return(out)
  }

  vals <- as.integer(as.vector(individual_results) > 0L)

  if (length(individual_results) == N) {
    return(vals)
  }

  nm <- names(individual_results)
  if (!is.null(nm)) {
    idx <- suppressWarnings(as.integer(nm))
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > N)) {
      stop("Named `individual_results` must use valid individual indices.")
    }
    out[idx] <- vals
    return(out)
  }

  if (length(vals) == length(planned)) {
    out[planned] <- vals
    return(out)
  }

  stop(
    "`individual_results` must have length N, use names for individual indices, ",
    "or have length equal to the planned retest set."
  )
}


.workflow_finalize_calls <- function(stage1,
                                     planned_stage2,
                                     stage2_results,
                                     N,
                                     workflow,
                                     stage2,
                                     finalize) {
  if (!isTRUE(finalize)) {
    return(list(
      positive = integer(0L),
      negative = integer(0L),
      unresolved = seq_len(N)
    ))
  }

  if (workflow == "non_adaptive" || stage2 == "none") {
    positive <- sort(unique(stage1$positive))
    negative <- sort(unique(stage1$negative))
    unresolved <- setdiff(seq_len(N), sort(unique(c(positive, negative))))
    return(list(
      positive = positive,
      negative = negative,
      unresolved = unresolved
    ))
  }

  observed_stage2 <- planned_stage2[!is.na(stage2_results[planned_stage2])]
  positive <- observed_stage2[stage2_results[observed_stage2] == 1L]
  negative <- sort(unique(c(
    stage1$negative,
    observed_stage2[stage2_results[observed_stage2] == 0L]
  )))
  unresolved <- setdiff(seq_len(N), sort(unique(c(positive, negative))))

  list(
    positive = sort(unique(positive)),
    negative = negative,
    unresolved = unresolved
  )
}
