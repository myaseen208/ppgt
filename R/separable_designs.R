# =============================================================================
# Combinatorial verification: d-separable and d-disjunct properties
# Internal functions only — not exported.
# References:
#   Kautz & Singleton (1964) — d-disjunct matrices
#   D'yachkov & Rykov (1982) — d-separable / superimposed codes
# =============================================================================


# -----------------------------------------------------------------------------
# Private helpers
# -----------------------------------------------------------------------------

# Coerce \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} to a plain logical matrix
# with rows indexed by \eqn{j \in \{1,\ldots,J\}} and columns indexed by
# \eqn{i \in \{1,\ldots,N\}}.
.to_dense_logical <- function(M) {
  as.matrix(M) != 0L
}

# For an active set \eqn{\mathcal{S} \subseteq \{1,\ldots,N\}} given by `idx`,
# compute the induced pool-status vector
# \eqn{\mathbf{u}_{\mathcal{S}} = \bigvee_{i \in \mathcal{S}} \mathbf{M}_{\cdot i}}
# in \eqn{\{0,1\}^{J}}.
.row_union <- function(M_dense, idx) {
  if (length(idx) == 1L) {
    as.integer(M_dense[, idx])
  } else {
    as.integer(rowSums(M_dense[, idx, drop = FALSE]) > 0L)
  }
}

# Collapse an integer (0/1) vector to a single string key for hash lookup.
.fp_key <- function(v) {
  paste(v, collapse = "")
}

# Attach standard attributes to a logical result, with `m = J` and `n = N`.
.make_result <- function(value, d, m, n, verified) {
  attr(value, "d")         <- d
  attr(value, "n_pools")   <- m
  attr(value, "n_samples") <- n
  attr(value, "verified")  <- verified
  value
}

# Constant: exact-enumeration threshold (number of size-d subsets).
.EXACT_THRESHOLD <- 1e5


# -----------------------------------------------------------------------------
# .is_d_separable
# -----------------------------------------------------------------------------

#' @title Check Whether a Pooling Matrix Is d-Separable
#'
#' @description
#' A binary pooling design matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}
#' (\eqn{J} pools, \eqn{N} individuals) is \eqn{d}-**separable** if no two
#' distinct active sets \eqn{\mathcal{S}_1, \mathcal{S}_2 \subseteq
#' \{1,\ldots,N\}} with \eqn{|\mathcal{S}_1|, |\mathcal{S}_2| \leq d} produce
#' the same union activation vector
#' \eqn{\mathbf{u}_\mathcal{S} = \bigvee_{i \in \mathcal{S}} \mathbf{m}_i}.
#' Equivalently, the observed pool-outcome vector \eqn{\mathbf{z}} uniquely
#' identifies the active set whenever at most \eqn{d} individuals are positive.
#'
#' @details
#' # Algorithm
#'
#' **Exact** (when \eqn{\binom{N}{d} \leq 10^5}): enumerate every non-empty
#' active set \eqn{\mathcal{S} \subseteq \{1,\ldots,N\}} with
#' \eqn{|\mathcal{S}| \in \{1,\ldots,d\}}, compute the induced pool-status
#' vector
#' \deqn{
#'   \mathbf{u}_{\mathcal{S}}
#'   =
#'   \bigvee_{i \in \mathcal{S}} \mathbf{M}_{\cdot i}
#'   \in \{0,1\}^{J},
#' }
#' and store a fingerprint of \eqn{\mathbf{u}_{\mathcal{S}}} in a hash table.
#' A collision means that two distinct active sets induce the same pool outcome
#' pattern and are therefore not distinguishable from \eqn{\mathbf{z}}.
#'
#' **Approximate** (when \eqn{\binom{N}{d} > 10^5}): draw \eqn{10^5} random
#' active-set pairs
#' \eqn{(\mathcal{S}_1, \mathcal{S}_2)} with
#' \eqn{\mathcal{S}_1, \mathcal{S}_2 \subseteq \{1,\ldots,N\}},
#' \eqn{|\mathcal{S}_1|, |\mathcal{S}_2| \leq d}, and
#' \eqn{\mathcal{S}_1 \neq \mathcal{S}_2}, and test whether
#' \eqn{\mathbf{u}_{\mathcal{S}_1} = \mathbf{u}_{\mathcal{S}_2}}. Returns
#' \code{FALSE} as soon as a collision is found; returns \code{TRUE} if none are
#' (conservative: \code{attr(., "verified") = "approximate"}).
#'
#' # Attributes on the Return Value
#'
#' | Attribute | Content |
#' |-----------|---------|
#' | \code{d} | The \eqn{d} tested |
#' | \code{n_pools} | \eqn{J} = number of pools (rows of \eqn{\mathbf{M}}) |
#' | \code{n_samples} | \eqn{N} = number of individuals (columns of \eqn{\mathbf{M}}) |
#' | \code{verified} | \code{"exact"} or \code{"approximate"} |
#'
#' @param M Pooling design matrix
#'   \eqn{\mathbf{M} = (M_{ji}) \in \{0,1\}^{J \times N}}, where row
#'   \eqn{j \in \{1,\ldots,J\}} indexes pool \eqn{\mathcal{P}_j} and column
#'   \eqn{i \in \{1,\ldots,N\}} indexes individual \eqn{i}. Accepted classes:
#'   \code{matrix}, \code{dgCMatrix}, or any \code{Matrix} subclass.
#' @param d Positive integer giving the largest active-set cardinality
#'   \eqn{|\mathcal{S}| \leq d} to verify.
#'
#' @return A single \code{logical} scalar with the attributes listed above.
#'
#' @examples
#' # Identity matrix: each individual i belongs to exactly one unique pool j.
#' # Any two distinct active sets have disjoint pool memberships, so it is always separable.
#' I5 <- diag(5)
#' res <- ppgt:::.is_d_separable(I5, d = 2)
#' stopifnot(isTRUE(res), attr(res, "verified") == "exact")
#'
#' # Two identical columns (individuals with identical pool membership) are
#' # indistinguishable, so it is NOT 1-separable.
#' M_dup <- cbind(c(1, 0, 1), c(1, 0, 1), c(0, 1, 0))
#' stopifnot(isFALSE(ppgt:::.is_d_separable(M_dup, d = 1)))
#'
#' # P-BEST: N = 384 individuals, J = 48 pools (approximate for d = 3)
#' \dontrun{
#'   M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#'   ppgt:::.is_d_separable(M, d = 2)   # exact
#'   ppgt:::.is_d_separable(M, d = 3)   # approximate (choose(384,3) >> 1e5)
#' }
#'
#' @references
#' D'yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of disjunctive
#' codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.
#'
#' @keywords internal
.is_d_separable <- function(M, d) {
  d       <- as.integer(d)
  M_dense <- .to_dense_logical(M)
  m       <- nrow(M_dense)
  n       <- ncol(M_dense)

  if (d < 1L) stop("d must be a positive integer")
  if (n < 1L) stop("M must have at least one column")

  exact <- choose(n, d) <= .EXACT_THRESHOLD

  if (!exact) {
    warning(
      "choose(", n, ", ", d, ") > 1e5: using Monte Carlo approximation ",
      "(1e5 random subset pairs). ",
      "attr(result, \"verified\") == \"approximate\"."
    )

    n_mc <- 1e5L
    for (i in seq_len(n_mc)) {
      k1 <- sample.int(d, 1L)
      k2 <- sample.int(d, 1L)
      S1 <- sort(sample.int(n, k1))
      S2 <- sort(sample.int(n, k2))
      if (identical(S1, S2)) next

      fp1 <- .row_union(M_dense, S1)
      fp2 <- .row_union(M_dense, S2)
      if (identical(fp1, fp2)) {
        return(.make_result(FALSE, d, m, n, "approximate"))
      }
    }
    return(.make_result(TRUE, d, m, n, "approximate"))
  }

  # Exact: enumerate all active sets \mathcal{S} with 1 <= |\mathcal{S}| <= d,
  # then check whether any two induce the same J-vector of pool outcomes.
  seen <- new.env(hash = TRUE, parent = emptyenv())

  for (k in seq_len(min(d, n))) {
    combs <- combn(n, k)
    for (j in seq_len(ncol(combs))) {
      fp  <- .row_union(M_dense, combs[, j])
      key <- .fp_key(fp)
      if (exists(key, envir = seen, inherits = FALSE)) {
        return(.make_result(FALSE, d, m, n, "exact"))
      }
      assign(key, TRUE, envir = seen)
    }
  }

  .make_result(TRUE, d, m, n, "exact")
}


# -----------------------------------------------------------------------------
# .is_d_disjunct
# -----------------------------------------------------------------------------

#' @title Check Whether a Pooling Matrix Is d-Disjunct
#'
#' @description
#' A binary matrix \eqn{\mathbf{M}} (pools \eqn{\times} samples) is
#' \eqn{d}-**disjunct** if for every set \eqn{S} of \eqn{d} columns and
#' every column \eqn{c \notin S}, the row support of \eqn{c} is **not** a
#' subset of the union of row supports of \eqn{S}.
#'
#' Equivalently: every sample can be uniquely identified even if any \eqn{d}
#' other samples are simultaneously positive.  A \eqn{d}-disjunct matrix
#' guarantees zero false negatives for the COMP decoder when at most \eqn{d}
#' samples are truly positive.
#'
#' @details
#' # Algorithm
#'
#' **Exact** (when \eqn{\binom{N}{d} \leq 10^5}): for every size-\eqn{d}
#' active set \eqn{\mathcal{S} \subseteq \{1,\ldots,N\}} with
#' \eqn{|\mathcal{S}| = d}, compute the pooled union
#' \deqn{
#'   \mathbf{u}_{\mathcal{S}}
#'   =
#'   \bigvee_{i \in \mathcal{S}} \mathbf{M}_{\cdot i},
#' }
#' then for every individual \eqn{i^\star \notin \mathcal{S}} test whether
#' \eqn{\mathbf{M}_{\cdot i^\star} \leq \mathbf{u}_{\mathcal{S}}} entry-wise.
#' Equivalently, this checks whether
#' \eqn{\mathcal{J}_{i^\star} \subseteq \bigcup_{i \in \mathcal{S}}
#' \mathcal{J}_{i}}. Returns \code{FALSE} immediately upon finding a violating
#' pair \eqn{(\mathcal{S}, i^\star)}.
#'
#' **Approximate** (when \eqn{\binom{N}{d} > 10^5}): draw \eqn{10^5} random
#' pairs \eqn{(\mathcal{S}, i^\star)} with \eqn{|\mathcal{S}| = d} and
#' \eqn{i^\star \notin \mathcal{S}}, and test the same containment condition.
#' Returns \code{FALSE} as soon as a violation is found; otherwise returns
#' \code{TRUE} (conservative, \code{attr(., "verified") = "approximate"}).
#'
#' # Attributes on the Return Value
#'
#' | Attribute | Content |
#' |-----------|---------|
#' | \code{d} | The \eqn{d} tested |
#' | \code{n_pools} | \eqn{J} = number of pools (rows of \eqn{\mathbf{M}}) |
#' | \code{n_samples} | \eqn{N} = number of individuals (columns of \eqn{\mathbf{M}}) |
#' | \code{verified} | \code{"exact"} or \code{"approximate"} |
#'
#' @param M Pooling design matrix
#'   \eqn{\mathbf{M} = (M_{ji}) \in \{0,1\}^{J \times N}}, where row
#'   \eqn{j \in \{1,\ldots,J\}} indexes pool \eqn{\mathcal{P}_j} and column
#'   \eqn{i \in \{1,\ldots,N\}} indexes individual \eqn{i}. Accepted classes:
#'   \code{matrix}, \code{dgCMatrix}, or any \code{Matrix} subclass.
#' @param d Positive integer giving the disjunctness level: every active set
#'   \eqn{\mathcal{S}} with \eqn{|\mathcal{S}| = d} must fail to cover any
#'   outside individual \eqn{i^\star \notin \mathcal{S}}.
#'
#' @return A single \code{logical} scalar with the attributes listed above.
#'
#' @examples
#' # Identity matrix is 1-disjunct: each sample occupies a unique pool,
#' # so no other column's support can cover it.
#' I5 <- diag(5)
#' res <- ppgt:::.is_d_disjunct(I5, d = 1)
#' stopifnot(isTRUE(res), attr(res, "verified") == "exact")
#'
#' # A matrix of all ones is NOT 1-disjunct:
#' # for any column c and set S = {s}, support(c) = all rows is a subset of support(s).
#' M_ones <- matrix(1L, nrow = 3, ncol = 4)
#' stopifnot(isFALSE(ppgt:::.is_d_disjunct(M_ones, d = 1)))
#'
#' # 2-disjunct check on a small PP matrix
#' \dontrun{
#'   M <- pp_matrix(q = 4, d = 3, nl = 5)   # 20 pools, 64 samples
#'   ppgt:::.is_d_disjunct(M, d = 2)        # exact (choose(64,2) = 2016)
#' }
#'
#' @references
#' Kautz, W. H., & Singleton, R. C. (1964). Nonrandom binary superimposed
#' codes. *IEEE Transactions on Information Theory*, 10(4), 363-377.
#'
#' @keywords internal
.is_d_disjunct <- function(M, d) {
  d       <- as.integer(d)
  M_dense <- .to_dense_logical(M)
  m       <- nrow(M_dense)
  n       <- ncol(M_dense)

  if (d < 1L) stop("d must be a positive integer")
  if (n <= d)  stop("M must have more columns than d")

  exact <- choose(n, d) <= .EXACT_THRESHOLD

  if (!exact) {
    warning(
      "choose(", n, ", ", d, ") > 1e5: using Monte Carlo approximation ",
      "(1e5 random (S, c) pairs). ",
      "attr(result, \"verified\") == \"approximate\"."
    )

    n_mc    <- 1e5L
    all_idx <- seq_len(n)
    for (i in seq_len(n_mc)) {
      S       <- sample.int(n, d)
      others  <- setdiff(all_idx, S)
      c_idx   <- sample(others, 1L)
      union_S <- rowSums(M_dense[, S, drop = FALSE]) > 0L
      col_c   <- M_dense[, c_idx]
      if (all(!col_c | union_S)) {
        return(.make_result(FALSE, d, m, n, "approximate"))
      }
    }
    return(.make_result(TRUE, d, m, n, "approximate"))
  }

  # Exact: enumerate every size-d active set \mathcal{S}, then test whether any
  # outside individual i^* has \mathcal{J}_{i^*} covered by the union over
  # i in \mathcal{S}.
  all_idx <- seq_len(n)
  combs   <- combn(n, d)

  for (j in seq_len(ncol(combs))) {
    S       <- combs[, j]
    union_S <- if (d == 1L) {
      M_dense[, S]
    } else {
      rowSums(M_dense[, S, drop = FALSE]) > 0L
    }

    others <- setdiff(all_idx, S)
    for (c_idx in others) {
      col_c <- M_dense[, c_idx]
      # \mathcal{J}_{i^*} \subseteq \bigcup_{i \in \mathcal{S}} \mathcal{J}_i
      # iff every pool containing i^* is active in the pooled union.
      if (all(!col_c | union_S)) {
        return(.make_result(FALSE, d, m, n, "exact"))
      }
    }
  }

  .make_result(TRUE, d, m, n, "exact")
}


# =============================================================================
# SeparableMatrix — exported constructor
# =============================================================================

# -----------------------------------------------------------------------------
# Private construction helpers
# -----------------------------------------------------------------------------

# Check if n is prime (trial division).
.sm_is_prime <- function(n) {
  n <- as.integer(n)
  if (n < 2L) return(FALSE)
  if (n == 2L) return(TRUE)
  if (n %% 2L == 0L) return(FALSE)
  lim <- as.integer(sqrt(n)) + 1L
  if (lim < 3L) return(TRUE)
  for (f in seq.int(3L, lim, by = 2L)) {
    if (n %% f == 0L) return(FALSE)
  }
  TRUE
}

# Check if n is a prime power (p^e, e >= 1).
.sm_is_prime_power <- function(n) {
  n <- as.integer(n)
  if (n < 2L) return(FALSE)
  if (.sm_is_prime(n)) return(TRUE)
  lim <- as.integer(sqrt(n))
  for (p in 2L:lim) {
    if (!.sm_is_prime(p)) next
    pk <- p * p
    while (pk <= n) {
      if (pk == n) return(TRUE)
      pk <- pk * p
    }
  }
  FALSE
}

# Smallest prime power \eqn{q} satisfying \eqn{q \geq x}.
.sm_next_prime_power <- function(x) {
  n <- as.integer(ceiling(x))
  if (n < 2L) n <- 2L
  while (!.sm_is_prime_power(n)) n <- n + 1L
  n
}

# Convert the integer label `x - 1` to the coefficient vector
# \eqn{(c_0, \ldots, c_{k-1})} in \eqn{GF(q)^k}.
.sm_base_q_coeffs <- function(x, q, k) {
  powers <- as.integer(q)^(seq_len(k) - 1L)
  ((as.integer(x) - 1L) %/% powers) %% as.integer(q)
}

# Evaluate the polynomial
# \eqn{f_i(x) = c_0 + c_1 x + \cdots + c_{k-1} x^{k-1}} in \eqn{GF(q)} at a
# finite field point `point`; `NA_integer_` denotes the point at infinity, whose
# evaluation is the leading coefficient \eqn{c_{k-1}}.
.sm_poly_eval <- function(coeffs, point, gf_obj) {
  if (is.na(point)) {
    return(as.integer(coeffs[length(coeffs)]))
  }

  value <- 0L
  point_pow <- 1L
  for (k in seq_along(coeffs)) {
    term <- gf_mult(as.integer(coeffs[k]), point_pow, gf_obj)
    value <- gf_add(value, term, gf_obj)
    point_pow <- gf_mult(point_pow, as.integer(point), gf_obj)
  }
  as.integer(value)
}

# Smallest Reed-Solomon / Kautz-Singleton parameters for a d-separable design.
# Let k be the smallest integer with q^k >= N. A Reed-Solomon code of length L
# has agreement at most k - 1, so the Kautz-Singleton condition
# L > d (k - 1) yields d-disjunctness and therefore d-separability.
.sm_rs_params <- function(N, d) {
  q <- 2L
  repeat {
    if (.sm_is_prime_power(q)) {
      k <- 1L
      capacity <- q
      while (capacity < N) {
        k <- k + 1L
        capacity <- capacity * q
      }
      L <- d * (k - 1L) + 1L
      if (L <= q + 1L) {
        return(list(q = q, k = k, L = L, J = L * q))
      }
    }
    q <- q + 1L
  }
}

# Reed-Solomon / Kautz-Singleton construction. Individual \eqn{i} is encoded by
# a polynomial \eqn{f_i} over \eqn{GF(q)} of degree at most \eqn{k-1}, where
# \eqn{k} is the smallest integer satisfying \eqn{q^k \geq N}. Evaluate the
# codeword at \eqn{L = d(k-1)+1} projective points, then expand each q-ary
# symbol into a one-hot block of length q. Returns
# \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} with \eqn{J = L q}.
.build_rs_sep_matrix <- function(N, d) {
  N <- as.integer(N)
  d <- as.integer(d)
  params <- .sm_rs_params(N, d)
  q <- params$q
  k <- params$k
  L <- params$L

  gf_obj <- gf(q)
  eval_points <- c(NA_integer_, seq.int(0L, q - 1L))[seq_len(L)]
  J <- params$J

  rows_i <- integer(0L)
  cols_j <- integer(0L)

  for (i in seq_len(N)) {
    coeffs <- .sm_base_q_coeffs(i, q, k)
    for (ell in seq_len(L)) {
      value <- .sm_poly_eval(coeffs, eval_points[ell], gf_obj)
      rows_i <- c(rows_i, (ell - 1L) * q + value + 1L)
      cols_j <- c(cols_j, i)
    }
  }

  Matrix::sparseMatrix(
    i    = rows_i,
    j    = cols_j,
    x    = rep(1L, length(rows_i)),
    dims = c(J, N)
  )
}

# Run `expr` under a fixed RNG seed, saving and restoring the caller's
# global RNG state so the enclosing constructor has no visible side effect
# on the user's random-number stream.
.sd_with_fixed_seed <- function(seed, expr) {
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
  })
  set.seed(seed)
  force(expr)
}

# Bernoulli construction with independent entries
# \eqn{M_{ji} \sim \mathrm{Bernoulli}(p)} and
# \eqn{p = (d \log N)^{1/d} / N}.
.build_random_sep_matrix <- function(N, d, J, seed) {
  N    <- as.integer(N)
  d    <- as.integer(d)
  J    <- as.integer(J)
  seed <- as.integer(seed)

  .sd_with_fixed_seed(seed, {
    p <- ((d * log(N))^(1.0 / d)) / N
    p <- min(max(p, 0.0), 1.0)
    draw <- matrix(stats::rbinom(J * N, 1L, p), nrow = J, ncol = N)
    nz <- which(draw != 0L, arr.ind = TRUE)

    Matrix::sparseMatrix(
      i    = nz[, 1L],
      j    = nz[, 2L],
      x    = rep(1L, nrow(nz)),
      dims = c(J, N)
    )
  })
}

# Bernoulli construction for d-disjunct search with independent entries
# \eqn{M_{ji} \sim \mathrm{Bernoulli}(1/d)}.
.build_random_disjunct_matrix <- function(N, d, J, seed) {
  N    <- as.integer(N)
  d    <- as.integer(d)
  J    <- as.integer(J)
  seed <- as.integer(seed)

  .sd_with_fixed_seed(seed, {
    p <- min(1.0, 1.0 / d)
    draw <- matrix(stats::rbinom(J * N, 1L, p), nrow = J, ncol = N)
    nz <- which(draw != 0L, arr.ind = TRUE)

    Matrix::sparseMatrix(
      i    = nz[, 1L],
      j    = nz[, 2L],
      x    = rep(1L, nrow(nz)),
      dims = c(J, N)
    )
  })
}


# -----------------------------------------------------------------------------
# SeparableMatrix
# -----------------------------------------------------------------------------

#' @title Construct a \eqn{d}-Separable Pooling Design Matrix
#'
#' @description
#' Returns a binary pooling design matrix
#' \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} that is (or is designed to be)
#' \eqn{d}-**separable**: for any two distinct active sets
#' \eqn{\mathcal{S}_1, \mathcal{S}_2 \subseteq \{1,\ldots,N\}} with
#' \eqn{|\mathcal{S}_1|, |\mathcal{S}_2| \leq d}, the union pool-activation
#' vectors satisfy \eqn{\mathbf{u}_{\mathcal{S}_1} \neq \mathbf{u}_{\mathcal{S}_2}}.
#' The constructor supports automatic method selection, a random Bernoulli
#' design, and a Reed-Solomon-style construction over \eqn{GF(q)}.
#'
#' @details
#' Let \eqn{N} denote the number of individuals, let \eqn{J} denote the number
#' of pools, and let \eqn{\mathbf{M} = (M_{ji}) \in \{0,1\}^{J \times N}} be the
#' pooling design matrix, where
#' \eqn{M_{ji} = 1} if and only if individual \eqn{i \in \{1,\ldots,N\}}
#' belongs to pool \eqn{\mathcal{P}_j}. For an active set
#' \eqn{\mathcal{S} \subseteq \{1,\ldots,N\}} with \eqn{|\mathcal{S}| \leq d},
#' define the induced pool-status vector
#' \deqn{
#'   \mathbf{u}_{\mathcal{S}}
#'   =
#'   \bigvee_{i \in \mathcal{S}} \mathbf{M}_{\cdot i}
#'   \in \{0,1\}^{J}.
#' }
#' The matrix \eqn{\mathbf{M}} is \eqn{d}-separable if
#' \deqn{
#'   \mathbf{u}_{\mathcal{S}_1} = \mathbf{u}_{\mathcal{S}_2}
#'   \;\Longrightarrow\;
#'   \mathcal{S}_1 = \mathcal{S}_2
#'   \quad \text{for all }
#'   \mathcal{S}_1, \mathcal{S}_2 \subseteq \{1,\ldots,N\}
#'   \text{ with }
#'   |\mathcal{S}_1|, |\mathcal{S}_2| \leq d.
#' }
#'
#' # Reed-Solomon Construction (\code{method = "reed-solomon"})
#'
#' Let \eqn{q} be the smallest prime power for which there exists an integer
#' \eqn{k \geq 1} such that \eqn{q^k \geq N} and
#' \eqn{L = d(k-1)+1 \leq q+1}. Associate individual \eqn{i} with the
#' polynomial \eqn{f_i(x)} of degree at most \eqn{k-1} whose coefficient vector
#' is the base-\eqn{q} expansion of \eqn{i-1}:
#' \deqn{
#'   f_i(x)
#'   =
#'   c_0(i) + c_1(i)x + \cdots + c_{k-1}(i)x^{k-1},
#'   \quad
#'   i - 1 = \sum_{\ell=0}^{k-1} c_\ell(i) q^\ell.
#' }
#' Choose \eqn{L = d(k-1)+1} distinct evaluation points
#' \eqn{x_1, \ldots, x_L} in the projective line over \eqn{GF(q)} and index the
#' pools by \eqn{(\ell, v)} with \eqn{\ell \in \{1,\ldots,L\}} and
#' \eqn{v \in GF(q)}. Then \eqn{J = L q} and
#' \deqn{
#'   M_{ji} = 1
#'   \;\Longleftrightarrow\;
#'   f_i(x_\ell) = v,
#'   \quad j = (\ell-1)q + v + 1.
#' }
#' Because any two Reed-Solomon codewords agree in at most \eqn{k-1}
#' evaluation positions, the choice \eqn{L > d(k-1)} satisfies the
#' Kautz-Singleton disjunctness criterion, hence the resulting binary design is
#' \eqn{d}-disjunct and therefore \eqn{d}-separable. The argument `M` is
#' ignored for this method because \eqn{J} is fixed by \eqn{q}, \eqn{k}, and
#' \eqn{d}.
#'
#' # Random Bernoulli Construction (\code{method = "random"})
#'
#' If the requested number of pools is not supplied, use the theoretical lower
#' bound
#' \deqn{
#'   J \geq 2 \binom{d}{1} \log_2(N) = 2 d \log_2(N).
#' }
#' The random construction uses independent Bernoulli entries
#' \deqn{
#'   M_{ji} \sim \mathrm{Bernoulli}(p),
#'   \quad
#'   p = \frac{(d \log N)^{1/d}}{N}.
#' }
#' Up to 10 seeds are tried until \code{.is_d_separable} returns \code{TRUE},
#' or the retry budget is exhausted.
#'
#' # Automatic Selection (\code{method = "auto"})
#'
#' Uses the Reed-Solomon construction when \eqn{d \leq 3}; otherwise it falls
#' back to the random Bernoulli construction.
#'
#' @param N Positive integer; total number of individuals \eqn{N}.
#' @param d Positive integer giving the separability order \eqn{d}.
#' @param M Positive integer or \code{NULL}. This argument requests the number
#'   of pools, corresponding to \eqn{J} in the notation above. If \code{NULL},
#'   the random method uses \eqn{\lceil 2 d \log_2(N) \rceil}; the
#'   Reed-Solomon method ignores this argument and uses the Kautz-Singleton
#'   value \eqn{J = q[d(k-1)+1]}.
#' @param method Character string specifying the construction rule:
#'   \code{"auto"}, \code{"random"}, or \code{"reed-solomon"}.
#' @param seed Integer seed used for the random Bernoulli construction and its
#'   retry sequence. The Reed-Solomon method records the supplied seed without
#'   using it in construction.
#' @param verify Logical; if \code{TRUE}, run \code{.is_d_separable} on the
#'   returned matrix and attach the result as \code{attr(result, "separable_verified")}.
#'
#' @return A binary \code{dgCMatrix} with \eqn{J = nrow(\mathbf{M})}
#'   rows and \eqn{N = ncol(\mathbf{M})} columns. The result carries the
#'   attributes \code{d}, \code{method_used}, \code{seed}, and, when
#'   \code{verify = TRUE}, \code{separable_verified}. The attribute
#'   \code{separable_verified} is a single logical value returned by
#'   \code{.is_d_separable()} and may itself carry the metadata attributes
#'   \code{d}, \code{n_pools}, \code{n_samples}, and \code{verified} indicating
#'   whether the check was exact or approximate.
#'
#' @examples
#' M <- SeparableMatrix(N = 50L, d = 2L)
#' M
#' attr(M, "d")
#' attr(M, "method_used")
#'
#' M_rs <- SeparableMatrix(N = 16L, d = 2L, method = "reed-solomon")
#' M_rs
#' attr(M_rs, "separable_verified")
#'
#' @seealso
#'   \code{\link{.is_d_separable}} for standalone verification,
#'   \code{\link{.is_d_disjunct}} for the stronger disjunctness property,
#'   \code{\link{DisjunctMatrix}} for the stronger exported constructor, and
#'   \code{\link{pp_matrix}} for the P-BEST polynomial-pool design.
#'
#' @references
#' D'yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of disjunctive
#' codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.
#'
#' Kautz, W. H., & Singleton, R. C. (1964). Nonrandom binary superimposed
#' codes. *IEEE Transactions on Information Theory*, 10(4), 363-377.
#'
#' @export
SeparableMatrix <- function(N,
                            d,
                            M      = NULL,
                            method = c("auto", "random", "reed-solomon"),
                            seed   = 42L,
                            verify = TRUE) {
  N      <- as.integer(N)
  d      <- as.integer(d)
  seed   <- as.integer(seed)
  method <- match.arg(method)

  if (N < 1L) stop("N must be a positive integer")
  if (d < 1L) stop("d must be a positive integer")
  if (d >= N) stop("d must be less than N")
  if (!is.null(M)) {
    M <- as.integer(M)
    if (M < 1L) stop("M must be a positive integer")
  }

  # Choose the construction rule. For d <= 3, the requested Reed-Solomon path
  # is available because a prime power q >= N^(1/d) always exists.
  if (method == "auto") {
    method <- if (d <= 3L) "reed-solomon" else "random"
  }

  if (method == "random" && is.null(M)) {
    M <- as.integer(ceiling(2.0 * d * log2(N)))
  }

  mat <- NULL
  seed_used <- seed
  verified <- NULL

  if (method == "reed-solomon") {
    mat <- .build_rs_sep_matrix(N, d)
    if (verify) {
      verified <- .is_d_separable(mat, d)
    }
  } else {
    for (attempt in 0:9) {
      seed_used <- seed + attempt
      mat <- .build_random_sep_matrix(N, d, M, seed_used)
      if (!verify) {
        break
      }
      verified <- .is_d_separable(mat, d)
      if (isTRUE(as.logical(verified))) {
        break
      }
    }

    if (verify && !isTRUE(as.logical(verified))) {
      warning(
        "SeparableMatrix: random construction did not verify as ",
        d, "-separable after 10 attempts; returning the last draw."
      )
    }
  }

  attr(mat, "d") <- d
  attr(mat, "method_used") <- method
  attr(mat, "seed") <- seed_used
  if (verify) {
    attr(mat, "separable_verified") <- verified
  }

  mat
}


# -----------------------------------------------------------------------------
# DisjunctMatrix
# -----------------------------------------------------------------------------

#' @title Construct a \eqn{d}-Disjunct Pooling Design Matrix
#'
#' @description
#' Returns a binary pooling design matrix
#' \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} that is (or is designed to be)
#' \eqn{d}-disjunct: for every active set
#' \eqn{\mathcal{S} \subseteq \{1,\ldots,N\}} with \eqn{|\mathcal{S}| = d} and
#' every outside individual \eqn{i^\star \notin \mathcal{S}}, the pool set
#' \eqn{\mathcal{J}_{i^\star}} is not contained in
#' \eqn{\bigcup_{i \in \mathcal{S}} \mathcal{J}_i}. The constructor supports
#' automatic method selection, a random Bernoulli design, and a Reed-Solomon /
#' Kautz-Singleton construction over \eqn{GF(q)}.
#'
#' @details
#' Let \eqn{N} denote the number of individuals, let \eqn{J} denote the number
#' of pools, and let \eqn{\mathbf{M} = (M_{ji}) \in \{0,1\}^{J \times N}} be the
#' pooling design matrix. For each individual \eqn{i}, let
#' \eqn{\mathcal{J}_i = \{j : M_{ji} = 1\}} denote the set of pools containing
#' individual \eqn{i}. Then \eqn{\mathbf{M}} is \eqn{d}-disjunct if
#' \deqn{
#'   \mathcal{J}_{i^\star}
#'   \nsubseteq
#'   \bigcup_{i \in \mathcal{S}} \mathcal{J}_i
#'   \quad \text{for all }
#'   \mathcal{S} \subseteq \{1,\ldots,N\},\ |\mathcal{S}| = d,\ i^\star \notin \mathcal{S}.
#' }
#' Equivalently, for every \eqn{i^\star \notin \mathcal{S}} there exists a pool
#' \eqn{j \in \mathcal{J}_{i^\star}} with \eqn{j \notin \bigcup_{i \in \mathcal{S}}
#' \mathcal{J}_i}.
#'
#' The Reed-Solomon branch uses the same Kautz-Singleton expansion as
#' \code{SeparableMatrix()}. Let \eqn{q} be the smallest prime power for which
#' there exists an integer \eqn{k \geq 1} satisfying \eqn{q^k \geq N} and
#' \eqn{L = d(k-1)+1 \leq q+1}. Associate individual \eqn{i} with the
#' polynomial \eqn{f_i(x)} over \eqn{GF(q)} of degree at most \eqn{k-1}, then
#' evaluate at \eqn{L} projective points and expand each q-ary symbol into a
#' one-hot block of length \eqn{q}, giving \eqn{J = L q} pools.
#'
#' When \eqn{q > d}, the Kautz-Singleton construction yields a
#' \eqn{d}-disjunct matrix (Kautz and Singleton, 1964). The argument `M` is
#' ignored for this method because \eqn{J} is fixed by \eqn{q}, \eqn{k}, and
#' \eqn{d}.
#'
#' If the requested number of pools is not supplied, use the lower bound
#' \deqn{
#'   J \geq (d+1)^2 \log_2(N).
#' }
#' The random branch uses independent Bernoulli entries
#' \deqn{
#'   M_{ji} \sim \mathrm{Bernoulli}(1/d).
#' }
#' Up to 10 seeds are tried until \code{.is_d_disjunct} returns \code{TRUE}, or
#' the retry budget is exhausted.
#'
#' Uses the Reed-Solomon construction when the corresponding prime-power
#' parameter satisfies \eqn{q > d}; otherwise falls back to the random
#' Bernoulli construction.
#'
#' @section Construction Methods:
#' \describe{
#'   \item{\code{method = "reed-solomon"}}{Uses the Reed-Solomon /
#'   Kautz-Singleton path and ignores the user-supplied \code{M} argument
#'   because \eqn{J} is determined by the resulting field parameters.}
#'   \item{\code{method = "random"}}{Uses an independent Bernoulli design with
#'   default size \eqn{\lceil (d + 1)^2 \log_2(N) \rceil} when \code{M} is not
#'   supplied.}
#'   \item{\code{method = "auto"}}{Uses the Reed-Solomon path when its
#'   prime-power parameter satisfies \eqn{q > d}; otherwise falls back to the
#'   random Bernoulli path.}
#' }
#'
#' @param N Positive integer; total number of individuals \eqn{N}.
#' @param d Positive integer giving the disjunctness order \eqn{d}.
#' @param M Positive integer or \code{NULL}. This argument requests the number
#'   of pools, corresponding to \eqn{J} in the notation above. If \code{NULL},
#'   the random method uses \eqn{\lceil (d+1)^2 \log_2(N) \rceil}; the
#'   Reed-Solomon method ignores this argument and uses the Kautz-Singleton
#'   value \eqn{J = q[d(k-1)+1]}.
#' @param method Character string specifying the construction rule:
#'   \code{"auto"}, \code{"random"}, or \code{"reed-solomon"}.
#' @param seed Integer seed used for the random Bernoulli construction and its
#'   retry sequence. The Reed-Solomon method records the supplied seed without
#'   using it in construction.
#' @param verify Logical; if \code{TRUE}, run \code{.is_d_disjunct} on the
#'   returned matrix and attach the result as \code{attr(result, "disjunct_verified")}.
#'
#' @return A binary \code{dgCMatrix} with \eqn{J = nrow(\mathbf{M})}
#'   rows and \eqn{N = ncol(\mathbf{M})} columns. The result carries the
#'   attributes \code{d}, \code{method_used}, \code{seed}, and, when
#'   \code{verify = TRUE}, \code{disjunct_verified}. The attribute
#'   \code{disjunct_verified} is a single logical value returned by
#'   \code{.is_d_disjunct()} and may itself carry the metadata attributes
#'   \code{d}, \code{n_pools}, \code{n_samples}, and \code{verified} indicating
#'   whether the check was exact or approximate.
#'
#' @note Every \eqn{d}-disjunct matrix is also \eqn{d}-separable, but the
#'   converse need not hold. Longer derivations and comparisons belong in the
#'   separable/disjunct vignette rather than in this constructor help page.
#'
#' @examples
#' M <- DisjunctMatrix(N = 50L, d = 2L)
#' M
#' attr(M, "d")
#' attr(M, "method_used")
#'
#' M_rs <- DisjunctMatrix(N = 16L, d = 2L, method = "reed-solomon")
#' M_rs
#' attr(M_rs, "disjunct_verified")
#'
#' @seealso
#'   \code{\link{SeparableMatrix}} for the weaker separability guarantee,
#'   \code{\link{.is_d_disjunct}} for standalone verification,
#'   \code{\link{.is_d_separable}} for the related separability property, and
#'   \code{\link{pp_matrix}} for the polynomial-pool construction.
#'
#' @references
#' Kautz, W. H., & Singleton, R. C. (1964). Nonrandom binary superimposed
#' codes. *IEEE Transactions on Information Theory*, 10(4), 363-377.
#'
#' D'yachkov, A. G., & Rykov, V. V. (1982). Bounds on the length of disjunctive
#' codes. *Problemy Peredachi Informatsii*, 18(3), 7-13.
#'
#' @export
DisjunctMatrix <- function(N,
                           d,
                           M      = NULL,
                           method = c("auto", "random", "reed-solomon"),
                           seed   = 42L,
                           verify = TRUE) {
  N      <- as.integer(N)
  d      <- as.integer(d)
  seed   <- as.integer(seed)
  method <- match.arg(method)

  if (N < 1L) stop("N must be a positive integer")
  if (d < 1L) stop("d must be a positive integer")
  if (d >= N) stop("d must be less than N")
  if (!is.null(M)) {
    M <- as.integer(M)
    if (M < 1L) stop("M must be a positive integer")
  }

  rs_params <- .sm_rs_params(N, d)

  if (method == "auto") {
    method <- if (rs_params$q > d) "reed-solomon" else "random"
  }

  if (method == "reed-solomon" && rs_params$q <= d) {
    stop("reed-solomon construction requires q > d for the d-disjunct guarantee; got q = ",
         rs_params$q, " and d = ", d)
  }

  if (method == "random" && is.null(M)) {
    M <- as.integer(ceiling((d + 1L)^2 * log2(N)))
  }

  mat <- NULL
  seed_used <- seed
  verified <- NULL

  if (method == "reed-solomon") {
    mat <- .build_rs_sep_matrix(N, d)
    if (verify) {
      verified <- .is_d_disjunct(mat, d)
    }
  } else {
    for (attempt in 0:9) {
      seed_used <- seed + attempt
      mat <- .build_random_disjunct_matrix(N, d, M, seed_used)
      if (!verify) {
        break
      }
      verified <- .is_d_disjunct(mat, d)
      if (isTRUE(as.logical(verified))) {
        break
      }
    }

    if (verify && !isTRUE(as.logical(verified))) {
      warning(
        "DisjunctMatrix: random construction did not verify as ",
        d, "-disjunct after 10 attempts; returning the last draw."
      )
    }
  }

  attr(mat, "d") <- d
  attr(mat, "method_used") <- method
  attr(mat, "seed") <- seed_used
  if (verify) {
    attr(mat, "disjunct_verified") <- verified
  }

  mat
}
