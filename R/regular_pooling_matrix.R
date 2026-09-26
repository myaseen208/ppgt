#' Score a Pooling Matrix
#'
#' Summarize marginal and pairwise-overlap properties of a pooling matrix.
#'
#' @param A A binary pooling matrix with pools in rows and specimens in columns.
#' @return A list containing dimensions, marginal tables, overlap summaries,
#'   density, and an indicator for PP metadata.
#' @examples
#' A <- regular_pooling_matrix(12, 18, 6, 4, seed = 4L)
#' score_pooling_matrix(A)
#' @export
score_pooling_matrix <- function(A) {
  if (!is.matrix(A) && !inherits(A, "Matrix")) {
    stop("`A` must be a matrix or a Matrix object.")
  }
  rs <- Matrix::rowSums(A)
  cs <- Matrix::colSums(A)
  C <- Matrix::crossprod(A)
  Matrix::diag(C) <- 0
  col_vals <- as.vector(C)
  P <- Matrix::tcrossprod(A)
  Matrix::diag(P) <- 0
  pool_vals <- as.vector(P)
  col_positive <- col_vals[col_vals > 0]
  pool_positive <- pool_vals[pool_vals > 0]
  list(
    dim = dim(A),
    row_sum_table = table(rs),
    col_sum_table = table(cs),
    max_col_overlap = if (length(col_positive)) max(col_positive) else 0,
    n_col_pairs_overlap_gt_1 = sum(col_vals > 1) / 2,
    col_overlap_table = table(col_positive),
    max_pool_overlap = if (length(pool_positive)) max(pool_positive) else 0,
    pool_overlap_table = table(pool_positive),
    density = Matrix::nnzero(A) / prod(dim(A)),
    has_pp_attr = !is.null(attr(A, "pp"))
  )
}

.try_pp_matrix_for_target <- function(n_pools, n_samples, pool_size,
                                      pools_per_sample, d_values = 2:8) {
  if (pools_per_sample < 1L) return(NULL)
  q_raw <- n_pools / pools_per_sample
  if (!is.finite(q_raw) || q_raw != floor(q_raw)) return(NULL)
  q <- as.integer(q_raw)
  nl <- pools_per_sample
  for (d in d_values) {
    A <- tryCatch(
      pp_matrix(q = q, d = d, nl = nl, N = n_samples),
      error = function(e) NULL
    )
    if (is.null(A)) next
    exact <- identical(dim(A), c(n_pools, n_samples)) &&
      all(Matrix::rowSums(A) == pool_size) &&
      all(Matrix::colSums(A) == pools_per_sample)
    if (exact) return(A)
  }
  NULL
}

.regular_candidate <- function(n_pools, n_samples, pool_size,
                               pools_per_sample, overlap_optimized) {
  total_ones <- n_samples * pools_per_sample
  if (total_ones == 0L) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                dims = c(n_pools, n_samples)))
  }
  remaining <- rep.int(pool_size, n_pools)
  rows_i <- integer(total_ones)
  cols_j <- rep(seq_len(n_samples), each = pools_per_sample)
  pool_members <- vector("list", n_pools)
  position <- 1L

  for (sample_id in seq_len(n_samples)) {
    columns_left <- n_samples - sample_id + 1L
    mandatory <- which(remaining == columns_left)
    if (length(mandatory) > pools_per_sample) return(NULL)
    selected <- mandatory
    seen <- integer(n_samples)
    if (length(selected)) {
      for (pool in selected) {
        members <- pool_members[[pool]]
        if (length(members)) seen[members] <- seen[members] + 1L
      }
    }

    while (length(selected) < pools_per_sample) {
      candidates <- which(remaining > 0L)
      if (length(selected)) candidates <- setdiff(candidates, selected)
      if (!length(candidates)) return(NULL)
      if (!overlap_optimized) {
        best <- candidates[remaining[candidates] == max(remaining[candidates])]
      } else {
        severity <- integer(length(candidates))
        duplicates <- integer(length(candidates))
        for (idx in seq_along(candidates)) {
          members <- pool_members[[candidates[idx]]]
          if (length(members)) {
            new_overlap <- seen[members] + 1L
            severity[idx] <- max(new_overlap)
            duplicates[idx] <- sum(new_overlap > 1L)
          }
        }
        keep <- severity == min(severity)
        keep <- keep & duplicates == min(duplicates[keep])
        best_remaining <- max(remaining[candidates[keep]])
        best <- candidates[keep & remaining[candidates] == best_remaining]
      }
      chosen <- if (length(best) == 1L) best else sample(best, 1L)
      selected <- c(selected, chosen)
      members <- pool_members[[chosen]]
      if (length(members)) seen[members] <- seen[members] + 1L
    }

    idx <- seq.int(position, length.out = pools_per_sample)
    rows_i[idx] <- selected
    remaining[selected] <- remaining[selected] - 1L
    for (pool in selected) {
      pool_members[[pool]] <- c(pool_members[[pool]], sample_id)
    }
    position <- position + pools_per_sample
  }
  if (any(remaining != 0L)) return(NULL)
  Matrix::sparseMatrix(i = rows_i, j = cols_j,
                       x = rep.int(1L, total_ones),
                       dims = c(n_pools, n_samples))
}

.score_key <- function(score) {
  c(score$max_col_overlap, score$n_col_pairs_overlap_gt_1,
    score$max_pool_overlap)
}

.is_better_score <- function(candidate, incumbent) {
  if (is.null(incumbent)) return(TRUE)
  candidate <- .score_key(candidate)
  incumbent <- .score_key(incumbent)
  first_difference <- which(candidate != incumbent)[1L]
  !is.na(first_difference) &&
    candidate[first_difference] < incumbent[first_difference]
}

#' Construct a Smart Regular Pooling Matrix
#'
#' Construct a sparse binary pooling matrix with prescribed constant row and
#' column sums. By default, an exact polynomial-pools construction is returned
#' when available; otherwise an overlap-aware regular construction is used.
#'
#' @param n_pools Integer number of pools (rows).
#' @param n_samples Integer number of samples (columns).
#' @param pool_size Required number of samples in every pool.
#' @param pools_per_sample Required number of pools containing every sample.
#' @param seed Optional integer seed for randomized fallback construction.
#' @param max_tries Maximum randomized candidates for the optimized fallback.
#' @param method One of `"auto"`, `"pp"`, `"regular"`, or
#'   `"overlap_optimized"`. Automatic mode tries PP first. PP mode requires an
#'   exact PP construction. Regular mode uses one generic realization.
#' @param prefer_pp In automatic mode, try an exact PP construction first?
#' @param d_values Polynomial dimensions tried, in order, during PP detection.
#'
#' @return A sparse binary matrix with the requested dimensions and marginals.
#'   An exact PP result is returned directly and retains its `"pp"` attribute;
#'   a generic result has a `"design"` attribute.
#'
#' @details
#' Feasibility requires
#' \deqn{n_{pools} pool\_size = n_{samples} pools\_per\_sample,}
#' `pool_size <= n_samples`, and `pools_per_sample <= n_pools`.
#'
#' PP detection derives `q = n_pools / pools_per_sample` and
#' `nl = pools_per_sample`, calls [pp_matrix()] for each requested `d`, and
#' verifies all dimensions and row and column sums. The generic fallback keeps
#' exact marginals and ranks candidates by maximum specimen overlap, number of
#' specimen pairs with overlap greater than one, and maximum pool overlap.
#'
#' Generic fallback matrices do not automatically claim PP, P-BEST,
#' disjunctness, separability, or protocol-specific decoding guarantees. Use
#' [score_pooling_matrix()] or other diagnostics to evaluate overlaps.
#'
#' @examples
#' M1 <- regular_pooling_matrix(94, 188, 8, 4, seed = 4L)
#' M2 <- regular_pooling_matrix(186, 372, 8, 4, seed = 4L)
#' M3 <- regular_pooling_matrix(188, 282, 6, 4, seed = 4L)
#' attr(M1, "pp")
#' attr(M2, "pp")
#' attr(M3, "pp")
#' M3_pp <- pp_matrix(q = 47, d = 3, nl = 4, N = 282)
#' identical(as.matrix(M3), as.matrix(M3_pp))
#' score_pooling_matrix(M1)
#' score_pooling_matrix(M2)
#' score_pooling_matrix(M3)
#'
#' @seealso [pp_matrix()], [pbest_clinical_matrix()], [PoolMatrix()]
#' @family pooling_designs
#' @export
regular_pooling_matrix <- function(
    n_pools, n_samples, pool_size, pools_per_sample, seed = NULL,
    max_tries = 1000,
    method = c("auto", "pp", "regular", "overlap_optimized"),
    prefer_pp = TRUE, d_values = 2:8) {
  integer_scalar <- function(x, name, minimum) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x != floor(x) ||
        x < minimum) {
      stop("`", name, "` must be a single integer >= ", minimum, ".")
    }
    as.integer(x)
  }
  n_pools <- integer_scalar(n_pools, "n_pools", 1L)
  n_samples <- integer_scalar(n_samples, "n_samples", 1L)
  pool_size <- integer_scalar(pool_size, "pool_size", 0L)
  pools_per_sample <- integer_scalar(pools_per_sample, "pools_per_sample", 0L)
  max_tries <- integer_scalar(max_tries, "max_tries", 1L)
  if (!is.null(seed)) seed <- integer_scalar(seed, "seed", 0L)
  method <- match.arg(method)
  if (length(prefer_pp) != 1L || is.na(prefer_pp) || !is.logical(prefer_pp)) {
    stop("`prefer_pp` must be TRUE or FALSE.")
  }
  if (!length(d_values) || anyNA(d_values) || any(!is.finite(d_values)) ||
      any(d_values != floor(d_values)) || any(d_values < 2L)) {
    stop("`d_values` must be a non-empty vector of integers >= 2.")
  }
  d_values <- as.integer(d_values)

  if (pool_size > n_samples) {
    stop("`pool_size` cannot exceed `n_samples` in a binary matrix.")
  }
  if (pools_per_sample > n_pools) {
    stop("`pools_per_sample` cannot exceed `n_pools` in a binary matrix.")
  }
  ones_from_rows <- as.double(n_pools) * pool_size
  ones_from_columns <- as.double(n_samples) * pools_per_sample
  if (ones_from_rows != ones_from_columns) {
    stop("Incompatible marginals: n_pools * pool_size = ", ones_from_rows,
         ", but n_samples * pools_per_sample = ", ones_from_columns, ".")
  }

  try_pp <- method == "pp" || (method == "auto" && prefer_pp)
  if (try_pp) {
    pp_result <- .try_pp_matrix_for_target(
      n_pools, n_samples, pool_size, pools_per_sample, d_values
    )
    if (!is.null(pp_result)) return(pp_result)
    if (method == "pp") {
      stop("No exact pp_matrix construction exists for the requested target.")
    }
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed)) set.seed(seed)

  optimized <- method != "regular"
  tries <- if (optimized) max_tries else 1L
  best <- NULL
  best_score <- NULL
  attempts_used <- 0L
  for (attempt in seq_len(tries)) {
    candidate <- .regular_candidate(n_pools, n_samples, pool_size,
                                    pools_per_sample, optimized)
    if (is.null(candidate)) next
    attempts_used <- attempt
    candidate_score <- score_pooling_matrix(candidate)
    if (.is_better_score(candidate_score, best_score)) {
      best <- candidate
      best_score <- candidate_score
    }
    if (candidate_score$max_col_overlap <= 1) break
  }
  if (is.null(best)) {
    stop("Could not realize the requested regular matrix after ", tries,
         " randomized attempt", if (tries == 1L) "" else "s", ".")
  }
  if (!identical(dim(best), c(n_pools, n_samples)) ||
      !all(Matrix::rowSums(best) == pool_size) ||
      !all(Matrix::colSums(best) == pools_per_sample)) {
    stop("Internal error: fallback construction did not realize exact marginals.")
  }
  attr(best, "design") <- list(
    type = "regular_bipartite", n_pools = n_pools, n_samples = n_samples,
    pool_size = pool_size, pools_per_sample = pools_per_sample, seed = seed,
    method = if (optimized) "overlap_optimized" else "regular_random_ties",
    attempts = attempts_used, score = .score_key(best_score)
  )
  best
}

