# =============================================================================
# Hypergraph-based Designs (HyperDesign equivalent)
# Replicates functionality of HyperDesign package
# =============================================================================

#' @title Generate A HYPER-Style Pooling Matrix
#'
#' @description
#' Construct a binary HYPER-style pooling matrix by cyclically developing
#' combinatorial starter blocks for \eqn{q = 1}, \eqn{q = 2}, or \eqn{q = 3}.
#' The returned design matches the package's upstream-compatible HYPER generator
#' and produces a matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} with exactly
#' \eqn{q} ones in each column.
#'
#' @param n Integer scalar giving the total number of individuals. This is the
#'   column dimension \eqn{N} of the returned design matrix and must satisfy
#'   \eqn{N \ge 1}.
#' @param m Integer scalar giving the total number of pools. This is the row
#'   dimension \eqn{J} of the returned design matrix and must satisfy
#'   \eqn{J \ge 1}. Additional arithmetic constraints depend on \code{q}.
#' @param q Integer scalar giving the number of pools containing each individual.
#'   Allowed values are \code{1L}, \code{2L}, and \code{3L}; the default is
#'   \code{1L}. In the returned matrix,
#'   \eqn{w_i = |\mathcal{J}_i| = q} for every individual \eqn{i}.
#' @param reorder Logical scalar controlling whether the generated column blocks
#'   are reordered using the upstream HYPER grouping pattern. The default is
#'   \code{TRUE}. This changes column order only; it does not change the row set
#'   or the column weights.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where \eqn{J = m} and \eqn{N = n}. Every column sum equals
#'   \eqn{q}. The matrix carries \code{attr(., "design")} with named fields
#'   \code{type}, \code{n}, \code{m}, \code{q}, and \code{pools_per_sample}.
#'
#' @details
#' Let \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denote the returned pooling
#' matrix, with \eqn{J = m} pools and \eqn{N = n} individuals. The construction
#' assigns each individual \eqn{i} to a set \eqn{\mathcal{J}_i} of exactly
#' \eqn{q} pools, so that
#' \deqn{w_i = |\mathcal{J}_i| = q \qquad \text{for all } i \in \{1,\ldots,N\}.}
#' The design is generated from a canonical list of \eqn{q}-subsets of the pool
#' index set, and the columns cycle through this factor list when
#' \eqn{N} exceeds the number of available canonical factors.
#'
#' For \eqn{q = 1}, the factors are singleton pools and the design cycles through
#' the \eqn{J} pools. For \eqn{q = 2}, the construction uses cyclic development
#' of starter pairs over a projective line of size \eqn{m-1} augmented by an
#' infinite point, which requires
#' \deqn{m \equiv 0 \pmod 2.}
#' For \eqn{q = 3}, the construction uses developed projective orbits and
#' requires
#' \deqn{m \equiv 0 \pmod 6}
#' together with the primality condition
#' \deqn{m - 1 \text{ is prime}.}
#' Under all three constructions, the row sets
#' \eqn{\mathcal{P}_j = \{i : M_{ji} = 1\}} are intended to be balanced as evenly
#' as permitted by the cyclic factor development and by the finite truncation to
#' \eqn{N} columns. When \code{reorder = TRUE}, the columns are permuted by
#' factor blocks to match the upstream HYPER ordering convention.
#'
#' @examples
#' M1 <- hyper_matrix(n = 10, m = 5, q = 1)
#' M1
#' Matrix::colSums(M1)
#'
#' M2 <- hyper_matrix(n = 12, m = 6, q = 2)
#' M2
#' Matrix::colSums(M2)
#'
#' M3 <- hyper_matrix(n = 24, m = 12, q = 3)
#' M3
#' Matrix::colSums(M3)
#'
#' @references
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' McMahan, C. S., Tebbs, J. M., & Bilder, C. R. (2012). Informative Dorfman
#' screening. \emph{Biometrics}, 68(1), 287-296.
#'
#' @seealso
#' \code{\link{HyperDesign}}, \code{\link{hyper_ec_matrix}},
#' \code{\link{compare_all_designs}}
#'
#' @family pooling_designs
#' @export
hyper_matrix <- function(n, m, q = 1L, reorder = TRUE) {
  n <- as.integer(n)
  m <- as.integer(m)
  q <- as.integer(q)

  if (n < 1L) stop("n must be a positive integer")
  if (m < 1L) stop("m must be a positive integer")
  if (!q %in% 1:3) stop("q must be 1, 2, or 3")

  # Generate the canonical factorization used by HyperGen::hyperdesign().
  if (q == 1L) {
    factors <- .gen1_factors(m)
    num_factors <- length(factors)
  } else if (q == 2L) {
    if (m %% 2L != 0L) stop("For q = 2, m must be even. Got m = ", m)
    factors <- .gen2_factors(m)
    num_factors <- length(factors)
  } else {  # q == 3
    if (m %% 6L != 0L) stop("For q = 3, m must be divisible by 6. Got m = ", m)
    if (!is_prime(m - 1L)) stop("For q = 3, m-1 must be prime. Got m-1 = ", m - 1L)
    factors <- .gen3_factors(m)
    num_factors <- length(factors)
  }

  if (n > num_factors) {
    factor_index <- ((seq_len(n) - 1L) %% num_factors) + 1L
  } else {
    factor_index <- seq_len(n)
  }

  pools <- factors[factor_index]

  if (reorder && num_factors > 1L) {
    ordered <- vector("list", length(pools))
    pos <- 1L
    for (i in seq_len(num_factors)) {
      if (i <= length(pools)) {
        idx <- seq.int(from = i, to = length(pools), by = num_factors)
        for (j in idx) {
          ordered[[pos]] <- pools[[j]]
          pos <- pos + 1L
        }
      }
    }
    pools <- ordered
  }

  rows_i <- integer(n * q)
  cols_j <- integer(n * q)
  pos <- 1L

  for (j in seq_len(n)) {
    factor <- pools[[j]]
    pool_ids <- vapply(
      factor,
      FUN = function(elem) {
        if (is.infinite(elem)) {
          m
        } else {
          as.integer(elem) + 1L
        }
      },
      FUN.VALUE = integer(1L)
    )
    idx <- seq.int(from = pos, length.out = q)
    rows_i[idx] <- pool_ids
    cols_j[idx] <- j
    pos <- pos + q
  }

  result <- Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(m, n)
  )

  attr(result, "design") <- list(
    type = "hyper",
    n = n,
    m = m,
    q = q,
    pools_per_sample = q
  )

  result
}


# Internal: Generate canonical factors for q = 1
.gen1_factors <- function(m) {
  factors <- vector("list", m)
  for (j in seq_len(m - 1L)) {
    factors[[j]] <- list(j - 1L)
  }
  factors[[m]] <- list(Inf)
  factors
}


# Internal: Generate factors for q = 2
# Uses starter blocks: {0, ∞}, {1, r-1}, {2, r-2}, ...
.gen2_factors <- function(m) {
  if (m %% 2L != 0L) stop("m must be even for q = 2")

  r <- m - 1L  # Field size

  # Generate starter blocks (difference pairs)
  starter <- vector("list", (r + 1L) %/% 2L)

  # First starter: {0, ∞}
  starter[[1L]] <- c(0L, Inf)

  # Remaining starters: {i, -i} = {i, r-i} for i = 1, ..., (r-1)/2
  for (i in 1:((r - 1L) %/% 2L)) {
    starter[[i + 1L]] <- c(i, r - i)
  }

  factors <- vector("list", m * length(starter) - length(starter))
  idx <- 1L

  for (g in 0:(r - 1L)) {
    for (s in starter) {
      new_factor <- sapply(s, function(x) {
        if (is.infinite(x)) {
          Inf
        } else {
          (x + g) %% r
        }
      })
      factors[[idx]] <- as.list(unname(new_factor))
      idx <- idx + 1L
    }
  }

  factors
}


# Internal: Generate factors for q = 3
# Uses difference triples over GF(r) where r = m-1 is prime
.gen3_factors <- function(m) {
  if (m %% 6L != 0L) stop("m must be divisible by 6 for q = 3")

  r <- m - 1L

  if (!is_prime(r)) stop("m-1 must be prime for q = 3")

  omega <- .find_primitive_root(r)
  projective_points <- c(as.list(0:(r - 1L)), list(Inf))

  orbits <- list()
  keys_seen <- character(0)

  for (point in projective_points) {
    orbit <- .projective_orbit(point, r)
    orbit_key <- .factor_key(orbit)
    if (!orbit_key %in% keys_seen) {
      orbits[[length(orbits) + 1L]] <- orbit
      keys_seen <- c(keys_seen, orbit_key)
    }
  }

  lambdas <- vapply(
    seq_len((r - 1L) %/% 2L),
    FUN = function(i) .mod_pow(omega, i, r),
    FUN.VALUE = integer(1L)
  )

  factors <- vector("list", length(orbits) * length(lambdas) * r)
  idx <- 1L

  for (lambda in lambdas) {
    for (g in 0:(r - 1L)) {
      for (orbit in orbits) {
        transformed <- lapply(
          orbit,
          FUN = function(x) .projective_add(.projective_mult(lambda, x, r), g, r)
        )
        factors[[idx]] <- transformed
        idx <- idx + 1L
      }
    }
  }

  factors
}

.projective_orbit <- function(point, p) {
  orbit <- vector("list", 0L)
  current <- point
  max_iter <- p + 1L

  for (i in seq_len(max_iter)) {
    orbit[[length(orbit) + 1L]] <- current
    current <- .projective_transform(current, p)
    if (.projective_equal(current, orbit[[1L]])) {
      break
    }
  }

  orbit
}

.projective_transform <- function(x, p) {
  neg_prod <- .projective_neg(
    .projective_mult(.projective_inv(x, p), .projective_add(1L, x, p), p),
    p
  )
  neg_prod
}

.projective_add <- function(x, y, p) {
  if (is.infinite(x) && is.infinite(y)) {
    return(0)
  }
  if (is.infinite(x)) {
    return(Inf)
  }
  if (is.infinite(y)) {
    return(Inf)
  }
  (as.integer(x) + as.integer(y)) %% p
}

.projective_mult <- function(x, y, p) {
  if (is.infinite(x) && is.infinite(y)) {
    return(Inf)
  }
  if (is.infinite(x)) {
    return(if (as.integer(y) == 0L) 1L else Inf)
  }
  if (is.infinite(y)) {
    return(if (as.integer(x) == 0L) 1L else Inf)
  }
  (as.integer(x) * as.integer(y)) %% p
}

.projective_inv <- function(x, p) {
  if (is.infinite(x)) {
    return(0L)
  }
  if (as.integer(x) == 0L) {
    return(Inf)
  }
  .mod_pow(as.integer(x), p - 2L, p)
}

.projective_neg <- function(x, p) {
  if (is.infinite(x)) {
    return(Inf)
  }
  ((-as.integer(x)) %% p + p) %% p
}

.projective_equal <- function(x, y) {
  if (is.infinite(x) && is.infinite(y)) {
    return(TRUE)
  }
  if (is.infinite(x) || is.infinite(y)) {
    return(FALSE)
  }
  as.integer(x) == as.integer(y)
}

.factor_key <- function(factor) {
  values <- vapply(
    factor,
    FUN = function(x) if (is.infinite(x)) "Inf" else as.character(as.integer(x)),
    FUN.VALUE = character(1L)
  )
  paste(sort(values), collapse = ",")
}


# Internal: Find primitive root modulo prime p
.find_primitive_root <- function(p) {
  if (p == 2L) return(1L)

  # Factor p-1
  phi <- p - 1L
  factors <- .prime_factors(phi)

  for (g in 2:(p - 1L)) {
    is_primitive <- TRUE
    for (f in factors) {
      if (.mod_pow(g, phi %/% f, p) == 1L) {
        is_primitive <- FALSE
        break
      }
    }
    if (is_primitive) return(g)
  }

  stop("No primitive root found for p = ", p)
}


# Internal: Modular exponentiation
.mod_pow <- function(base, exp, mod) {
  result <- 1L
  base <- base %% mod

  while (exp > 0L) {
    if (exp %% 2L == 1L) {
      result <- (result * base) %% mod
    }
    exp <- exp %/% 2L
    base <- (base * base) %% mod
  }

  as.integer(result)
}


# Internal: Prime factorization
.prime_factors <- function(n) {
  factors <- integer(0)
  d <- 2L

  while (d * d <= n) {
    while (n %% d == 0L) {
      if (!d %in% factors) factors <- c(factors, d)
      n <- n %/% d
    }
    d <- d + 1L
  }

  if (n > 1L && !n %in% factors) {
    factors <- c(factors, n)
  }

  factors
}


#' @title HyperDesign Compatibility Wrapper
#'
#' @description
#' Provide a compatibility wrapper with the same argument interface as the former
#' \code{HyperDesign::HyperDesign()} constructor. This function forwards directly
#' to \code{\link{hyper_matrix}} and returns the same matrix object.
#'
#' @param n Integer scalar giving the total number of individuals. This becomes
#'   the column dimension \eqn{N} of the returned matrix.
#' @param m Integer scalar giving the total number of pools. This becomes the row
#'   dimension \eqn{J} of the returned matrix.
#' @param q Integer scalar giving the number of pools containing each individual.
#'   Allowed values are \code{1L}, \code{2L}, and \code{3L}.
#' @param reorder Logical scalar controlling whether the columns are reordered by
#'   the upstream HYPER grouping convention. The default is \code{TRUE}.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}. The object is exactly the output of
#'   \code{\link{hyper_matrix}(n = n, m = m, q = q, reorder = reorder)} and
#'   carries the same \code{design} attribute.
#'
#' @details
#' This wrapper exists for interface compatibility only. If
#' \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denotes the returned matrix, then
#' \eqn{J = m}, \eqn{N = n}, and every individual \eqn{i} satisfies
#' \deqn{w_i = |\mathcal{J}_i| = q.}
#' All arithmetic constraints on \code{m} and \code{q} are inherited unchanged
#' from \code{\link{hyper_matrix}}:
#' \deqn{q \in \{1,2,3\},}
#' with \eqn{m} even for \eqn{q = 2}, and \eqn{m \equiv 0 \pmod 6} together with
#' primality of \eqn{m-1} for \eqn{q = 3}. No additional computation is carried
#' out by this wrapper beyond forwarding the call.
#'
#' @examples
#' design1 <- HyperDesign(n = 10, m = 5, q = 1)
#' design1
#'
#' design2 <- HyperDesign(n = 12, m = 6, q = 2)
#' design2
#' Matrix::colSums(design2)
#'
#' @references
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' @seealso
#' \code{\link{hyper_matrix}}, \code{\link{compare_all_designs}},
#' \code{\link{compare_all_designs_honest}}
#'
#' @family pooling_designs
#' @export
HyperDesign <- function(n, m, q, reorder = TRUE) {
  hyper_matrix(n = n, m = m, q = q, reorder = reorder)
}
