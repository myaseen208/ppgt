# =============================================================================
# Additional Pooling Matrix Designs for Group Testing
# Covers: Dorfman, Array, IBD, Hypercube, Kirkman Triple, Projective Geometry
# =============================================================================

#' @title Additional Pooling Matrix Designs
#'
#' @description
#' This module implements additional pooling matrix designs beyond Polynomial
#' Pools, covering the complete taxonomy of group testing designs:
#'
#' 1. **Dorfman Designs** - Classic two-stage pooling
#' 2. **Array Designs** - Row/column grid pooling
#' 3. **IBD Designs** - Balanced incomplete block designs
#' 4. **Hypercube Designs** - Multi-dimensional grid pooling
#' 5. **Kirkman Triple Designs** - Resolvable triple systems
#' 6. **Projective Geometry Designs** - PG(d,q) based designs
#'
#' @name additional_designs
#' @family pooling_designs
NULL


# =============================================================================
# 1. DORFMAN DESIGNS
# =============================================================================

#' @title Generate A Dorfman Pooling Matrix
#'
#' @description
#' Construct the binary pooling design matrix for classical two-stage Dorfman
#' testing with non-overlapping pools of common target size \code{g}.
#'
#' @param N Integer scalar giving the total number of individuals. This is the
#'   column dimension \eqn{N} of the returned design matrix
#'   \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} and must satisfy \eqn{N \ge 1}.
#' @param g Integer scalar giving the target pool size. It must satisfy
#'   \eqn{2 \le g \le N}. In the returned design, each pool \eqn{\mathcal{P}_j}
#'   has size at most \eqn{g}, and all but possibly the last pool have size
#'   exactly \eqn{g}.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where \eqn{J = \lceil N / g \rceil}. Entry
#'   \eqn{M_{ji} = 1} if and only if individual \eqn{i} belongs to pool
#'   \eqn{j}. Every column sum is one, so \eqn{w_i = |\mathcal{J}_i| = 1} for
#'   every individual \eqn{i}. The matrix carries \code{attr(., "design")} with
#'   named fields \code{type}, \code{N}, \code{M}, \code{pool_size}, and
#'   \code{pools_per_sample}, where \code{M} is a legacy display label for the
#'   number of pools \eqn{J}.
#'
#' @details
#' Let \eqn{N} denote the total number of individuals. Dorfman testing partitions
#' the individuals into
#' \deqn{J = \left\lceil \frac{N}{g} \right\rceil}
#' disjoint pools. The returned matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}
#' satisfies
#' \deqn{\mathcal{P}_j = \{i : M_{ji} = 1\}, \qquad \mathcal{J}_i = \{j :
#' M_{ji} = 1\},}
#' with
#' \deqn{|\mathcal{J}_i| = 1 \quad \text{for all } i \in \{1,\ldots,N\}.}
#' Hence the stage-1 pools form a partition of the individual set. If
#' \eqn{\widetilde{y}_i \in \{0,1\}} denotes the latent status of individual
#' \eqn{i}, then the latent status of pool \eqn{j} is
#' \deqn{\widetilde{z}_j = \mathbb{I}\left(\sum_{i \in \mathcal{P}_j}
#' \widetilde{y}_i > 0\right).}
#' Under a homogeneous prevalence model
#' \eqn{\Pr(\widetilde{y}_i = 1) = p} and perfect testing, the expected number of
#' tests per individual in the two-stage procedure is
#' \deqn{\frac{\mathbb{E}[T]}{N} = \frac{1}{g} + 1 - (1-p)^g.}
#' The heuristic choice
#' \deqn{g_{\mathrm{opt}} \approx p^{-1/2}}
#' motivates \code{\link{dorfman_optimal_g}}.
#'
#' @examples
#' M <- dorfman_matrix(N = 12, g = 3)
#' M
#' Matrix::rowSums(M)
#' Matrix::colSums(M)
#'
#' g_opt <- dorfman_optimal_g(0.01)
#' g_opt
#' M_opt <- dorfman_matrix(N = 33, g = g_opt)
#' M_opt
#'
#' @references
#' Dorfman, R. (1943). The detection of defective members of large populations.
#' \emph{The Annals of Mathematical Statistics}, 14(4), 436-440.
#'
#' Sobel, M., & Groll, P. A. (1959). Group testing to eliminate efficiently all
#' defectives in a binomial sample. \emph{Bell System Technical Journal}, 38(5),
#' 1179-1252.
#'
#' @seealso
#' \code{\link{dorfman_optimal_g}}, \code{\link{array_matrix}},
#' \code{\link{compare_designs}}
#'
#' @family pooling_designs
#' @export
dorfman_matrix <- function(N, g) {
  N <- as.integer(N)
  g <- as.integer(g)

  if (N < 1L) stop("N must be a positive integer")
  if (g < 2L) stop("Pool size g must be at least 2")
  if (g > N) stop("Pool size g cannot exceed N")

  M <- ceiling(N / g)

  # Assign each sample to exactly one pool
  pool_assignments <- rep(1:M, each = g)[1:N]

  result <- Matrix::sparseMatrix(
    i = pool_assignments,
    j = 1:N,
    x = rep(1L, N),
    dims = c(M, N)
  )

  attr(result, "design") <- list(
    type = "dorfman",
    N = N,
    M = M,
    pool_size = g,
    pools_per_sample = 1L
  )

  result
}


#' @title Compute An Approximate Optimal Dorfman Pool Size
#'
#' @description
#' Choose the pool size \eqn{g} that minimizes the expected number of tests per
#' individual under the homogeneous Dorfman model.
#'
#' @param p Numeric scalar giving the common individual prevalence
#'   \eqn{p_i = p \in (0,1)}. It must satisfy \eqn{0 < p < 1}.
#'
#' @return An integer scalar giving the pool size \eqn{g} that minimizes the
#'   objective evaluated by this function over a local search range around
#'   \eqn{1 / \sqrt{p}}.
#'
#' @details
#' In two-stage Dorfman testing, a pool of size \eqn{g} requires one stage-1
#' pooled test and then \eqn{g} follow-up individual tests whenever the pooled
#' result is positive. Under a homogeneous prevalence model
#' \eqn{\Pr(\widetilde{y}_i = 1) = p} with independent individual statuses and
#' perfect testing, the expected number of tests per individual is
#' \deqn{\psi(g; p) = \frac{1}{g} + 1 - (1-p)^g.}
#' This function numerically minimizes \eqn{\psi(g; p)} over integer values in a
#' finite window centered near the classical approximation
#' \deqn{g_{\mathrm{opt}} \approx \frac{1}{\sqrt{p}}.}
#' The returned value can be passed directly to \code{\link{dorfman_matrix}}.
#'
#' @examples
#' g1 <- dorfman_optimal_g(0.01)
#' g1
#' g2 <- dorfman_optimal_g(0.05)
#' g2
#' g3 <- dorfman_optimal_g(0.001)
#' g3
#'
#' @references
#' Dorfman, R. (1943). The detection of defective members of large populations.
#' \emph{The Annals of Mathematical Statistics}, 14(4), 436-440.
#'
#' Sobel, M., & Groll, P. A. (1959). Group testing to eliminate efficiently all
#' defectives in a binomial sample. \emph{Bell System Technical Journal}, 38(5),
#' 1179-1252.
#'
#' @seealso
#' \code{\link{dorfman_matrix}}, \code{\link{compare_designs}}
#'
#' @family pooling_designs
#' @export
dorfman_optimal_g <- function(p) {
  if (p <= 0 || p >= 1) stop("Prevalence must be in (0, 1)")

  # Expected tests per sample
  expected_cost <- function(g, p) {
    1/g + 1 - (1 - p)^g
  }

  # Search around 1/sqrt(p)
  g_approx <- 1 / sqrt(p)
  g_range <- max(2, floor(g_approx) - 3):min(ceiling(g_approx) + 3, 200)

  costs <- sapply(g_range, expected_cost, p = p)
  g_opt <- g_range[which.min(costs)]

  as.integer(g_opt)
}


# =============================================================================
# 2. ARRAY DESIGNS
# =============================================================================

#' @title Generate A Two-Dimensional Array Pooling Matrix
#'
#' @description
#' Construct a row-column array design in which the \eqn{N = rc} individuals are
#' placed on an \eqn{r \times c} grid and pooled by grid rows and columns.
#'
#' @param r Integer scalar giving the number of grid rows. It must satisfy
#'   \eqn{r \ge 2}. This parameter determines the first family of pools and
#'   contributes \eqn{r} rows to \eqn{\mathbf{M}}.
#' @param c Integer scalar giving the number of grid columns. It must satisfy
#'   \eqn{c \ge 2}. This parameter determines the second family of pools and
#'   contributes \eqn{c} rows to \eqn{\mathbf{M}}.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where \eqn{J = r + c} and \eqn{N = rc}. Every individual
#'   belongs to exactly one row pool and one column pool, so
#'   \eqn{w_i = |\mathcal{J}_i| = 2}. The attached \code{design} attribute is a
#'   named list with fields \code{type}, \code{N}, \code{M}, \code{rows},
#'   \code{cols}, \code{pool_size_row}, \code{pool_size_col}, and
#'   \code{pools_per_sample}, where \code{M} is a legacy display label for the
#'   number of pools \eqn{J}.
#'
#' @details
#' Index the individuals by ordered pairs \eqn{(u,v)} with
#' \eqn{u \in \{1,\ldots,r\}} and \eqn{v \in \{1,\ldots,c\}}. After linearizing
#' the grid into column indices \eqn{i \in \{1,\ldots,N\}}, the returned design
#' matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} has two pool families:
#' \deqn{\mathcal{P}^{\mathrm{row}}_u = \{i : \text{individual } i \text{ lies in
#' row } u\}, \qquad u = 1,\ldots,r,}
#' \deqn{\mathcal{P}^{\mathrm{col}}_v = \{i : \text{individual } i \text{ lies in
#' column } v\}, \qquad v = 1,\ldots,c.}
#' Hence
#' \deqn{J = r + c, \qquad n_j \in \{r, c\}, \qquad w_i = 2.}
#' For a single positive individual, the unique positive row and unique positive
#' column identify that individual by intersection. With multiple positives, the
#' stage-1 positive pattern may correspond to several candidate intersections, so
#' the design is primarily useful as a low-weight screening matrix rather than an
#' exact combinatorial code.
#'
#' @examples
#' M <- array_matrix(r = 4, c = 5)
#' M
#' Matrix::rowSums(M)
#' Matrix::colSums(M)
#'
#' M96 <- array_matrix(r = 8, c = 12)
#' M96
#'
#' @references
#' Phatarfod, R. M., & Sudbury, A. (1994). The use of a square array scheme in
#' blood testing. \emph{Statistics in Medicine}, 13(22), 2337-2343.
#'
#' Kim, H. Y., Hudgens, M. G., Dreyfuss, J. M., Westreich, D. J., & Pilcher,
#' C. D. (2007). Comparison of group testing algorithms for case identification
#' in the presence of test error. \emph{Biometrics}, 63(4), 1152-1163.
#'
#' @seealso
#' \code{\link{array_matrix_3d}}, \code{\link{hypercube_matrix}},
#' \code{\link{compare_designs}}
#'
#' @family pooling_designs
#' @export
array_matrix <- function(r, c) {
  r <- as.integer(r)
  c <- as.integer(c)

  if (r < 2 || c < 2) stop("Both r and c must be at least 2")

  N <- r * c
  M <- r + c

  rows_i <- integer(0)
  cols_j <- integer(0)

  # Row pools (first r pools)
  for (i in 1:r) {
    sample_ids <- ((i - 1) * c + 1):(i * c)
    rows_i <- c(rows_i, rep(i, c))
    cols_j <- c(cols_j, sample_ids)
  }

  # Column pools (next c pools)
  for (j in 1:c) {
    sample_ids <- seq(j, N, by = c)
    rows_i <- c(rows_i, rep(r + j, r))
    cols_j <- c(cols_j, sample_ids)
  }

  result <- Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(M, N)
  )

  attr(result, "design") <- list(
    type = "array_2d",
    N = N,
    M = M,
    rows = r,
    cols = c,
    pool_size_row = c,
    pool_size_col = r,
    pools_per_sample = 2L
  )

  result
}


#' @title Generate A Three-Dimensional Array Pooling Matrix
#'
#' @description
#' Construct a three-dimensional array design in which the \eqn{N = d_1 d_2 d_3}
#' individuals are arranged on a rectangular lattice and pooled by the three
#' coordinate directions.
#'
#' @param d1 Integer scalar giving the size of the first coordinate dimension.
#'   It must satisfy \eqn{d_1 \ge 2}. This contributes \eqn{d_1} pools to the
#'   returned matrix.
#' @param d2 Integer scalar giving the size of the second coordinate dimension.
#'   It must satisfy \eqn{d_2 \ge 2}. This contributes \eqn{d_2} pools to the
#'   returned matrix.
#' @param d3 Integer scalar giving the size of the third coordinate dimension.
#'   It must satisfy \eqn{d_3 \ge 2}. This contributes \eqn{d_3} pools to the
#'   returned matrix.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where \eqn{J = d_1 + d_2 + d_3} and
#'   \eqn{N = d_1 d_2 d_3}. Each individual lies in exactly one pool from each
#'   coordinate family, so \eqn{w_i = |\mathcal{J}_i| = 3}. The attached
#'   \code{design} attribute is a named list with fields \code{type},
#'   \code{N}, \code{M}, \code{dimensions}, and \code{pools_per_sample}, where
#'   \code{M} is a legacy display label for the number of pools \eqn{J}.
#'
#' @details
#' Label the individuals by triples \eqn{(u,v,w)} with
#' \eqn{u \in \{1,\ldots,d_1\}}, \eqn{v \in \{1,\ldots,d_2\}}, and
#' \eqn{w \in \{1,\ldots,d_3\}}. The design matrix
#' \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} is formed from three parallel pool
#' families:
#' \deqn{\mathcal{P}^{(1)}_u = \{i : \text{the first coordinate of } i \text{ is }
#' u\},}
#' \deqn{\mathcal{P}^{(2)}_v = \{i : \text{the second coordinate of } i \text{ is }
#' v\},}
#' \deqn{\mathcal{P}^{(3)}_w = \{i : \text{the third coordinate of } i \text{ is }
#' w\}.}
#' Therefore
#' \deqn{J = d_1 + d_2 + d_3}
#' and every individual participates in exactly three pools. The family sizes are
#' \deqn{|\mathcal{P}^{(1)}_u| = d_2 d_3, \qquad
#' |\mathcal{P}^{(2)}_v| = d_1 d_3, \qquad
#' |\mathcal{P}^{(3)}_w| = d_1 d_2.}
#' When there is a single positive individual, the unique positive pool in each
#' family identifies its three coordinates.
#'
#' @examples
#' M <- array_matrix_3d(3, 3, 2)
#' M
#' Matrix::rowSums(M)
#' Matrix::colSums(M)
#'
#' @references
#' Li, X., & Ying, K. (2013). Two- and three-dimensional pooling algorithms for
#' high-throughput screening assays. \emph{Technometrics}, 55(3), 285-296.
#'
#' Kim, H. Y., Hudgens, M. G., Dreyfuss, J. M., Westreich, D. J., & Pilcher,
#' C. D. (2007). Comparison of group testing algorithms for case identification
#' in the presence of test error. \emph{Biometrics}, 63(4), 1152-1163.
#'
#' @seealso
#' \code{\link{array_matrix}}, \code{\link{hypercube_matrix}},
#' \code{\link{compare_designs}}
#'
#' @family pooling_designs
#' @export
array_matrix_3d <- function(d1, d2, d3) {
  d1 <- as.integer(d1)
  d2 <- as.integer(d2)
  d3 <- as.integer(d3)

  if (d1 < 2 || d2 < 2 || d3 < 2) stop("All dimensions must be at least 2")

  N <- d1 * d2 * d3
  M <- d1 + d2 + d3

  # Sample index function: (i,j,k) -> linear index
  sample_idx <- function(i, j, k) (i - 1) * d2 * d3 + (j - 1) * d3 + k

  rows_i <- integer(0)
  cols_j <- integer(0)

  # Dimension 1 pools (fix i, vary j,k)
  for (i in 1:d1) {
    for (j in 1:d2) {
      for (k in 1:d3) {
        rows_i <- c(rows_i, i)
        cols_j <- c(cols_j, sample_idx(i, j, k))
      }
    }
  }

  # Dimension 2 pools (fix j, vary i,k)
  for (j in 1:d2) {
    for (i in 1:d1) {
      for (k in 1:d3) {
        rows_i <- c(rows_i, d1 + j)
        cols_j <- c(cols_j, sample_idx(i, j, k))
      }
    }
  }

  # Dimension 3 pools (fix k, vary i,j)
  for (k in 1:d3) {
    for (i in 1:d1) {
      for (j in 1:d2) {
        rows_i <- c(rows_i, d1 + d2 + k)
        cols_j <- c(cols_j, sample_idx(i, j, k))
      }
    }
  }

  result <- Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(M, N)
  )

  attr(result, "design") <- list(
    type = "array_3d",
    N = N,
    M = M,
    dimensions = c(d1, d2, d3),
    pools_per_sample = 3L
  )

  result
}


# =============================================================================
# 3. IBD DESIGNS (Incomplete Block Designs)
# =============================================================================

#' @title Generate A BIBD Pooling Matrix
#'
#' @description
#' Construct the incidence matrix of a balanced incomplete block design (BIBD)
#' when one of the built-in block constructions is available.
#'
#' @param v Integer scalar giving the number of treatments in the underlying
#'   BIBD. In the pooling interpretation this is the total number of individuals,
#'   so the returned matrix has \eqn{N = v} columns.
#' @param k Integer scalar giving the block size of the BIBD. In pooling terms,
#'   each pool \eqn{\mathcal{P}_j} has size \eqn{k}. It must satisfy
#'   \eqn{2 \le k < v}.
#' @param lambda Integer scalar giving the pairwise concurrence parameter
#'   \eqn{\lambda \ge 1}. Every unordered pair of distinct individuals appears
#'   together in exactly \eqn{\lambda} pools when the design exists. The default
#'   is \code{1L}.
#'
#' @return Either a sparse binary matrix of class \code{dgCMatrix} with
#'   dimensions \eqn{J \times N}, where \eqn{N = v} and \eqn{J = b}, or
#'   \code{NULL} if the required arithmetic conditions fail or no built-in
#'   construction is available. When a matrix is returned, the \code{design}
#'   attribute stores \code{type}, \code{v}, \code{b}, \code{r}, \code{k}, and
#'   \code{lambda}.
#'
#' @details
#' A \eqn{(v, b, r, k, \lambda)} balanced incomplete block design is represented
#' by an incidence matrix \eqn{\mathbf{M} \in \{0,1\}^{b \times v}}. In the
#' pooling interpretation, \eqn{N = v} and \eqn{J = b}. The defining relations
#' are
#' \deqn{n_j = |\mathcal{P}_j| = k \quad \text{for all } j,}
#' \deqn{w_i = |\mathcal{J}_i| = r \quad \text{for all } i,}
#' and every unordered pair of distinct individuals occurs together in exactly
#' \eqn{\lambda} pools. The standard parameter identities are
#' \deqn{bk = vr}
#' and
#' \deqn{\lambda(v-1) = r(k-1).}
#' This function first computes \eqn{r} and \eqn{b} from these identities and
#' then dispatches to a small collection of built-in constructions such as the
#' Fano plane. If no matching construction is coded, the function returns
#' \code{NULL}.
#'
#' @examples
#' M <- bibd_matrix(v = 7, k = 3)
#' M
#' Matrix::rowSums(M)
#' Matrix::colSums(M)
#'
#' M_fail <- bibd_matrix(v = 8, k = 3)
#' M_fail
#'
#' @references
#' Beth, T., Jungnickel, D., & Lenz, H. (1999). \emph{Design Theory}. Cambridge
#' University Press.
#'
#' Colbourn, C. J., & Dinitz, J. H. (Eds.). (2007). \emph{Handbook of
#' Combinatorial Designs} (2nd ed.). Chapman & Hall/CRC.
#'
#' @seealso
#' \code{\link{kirkman_matrix}}, \code{\link{pg_matrix}},
#' \code{\link{compare_designs}}
#'
#' @family pooling_designs
#' @export
bibd_matrix <- function(v, k, lambda = 1L) {
  v <- as.integer(v)
  k <- as.integer(k)
  lambda <- as.integer(lambda)

  if (k < 2) stop("k must be at least 2")
  if (k >= v) stop("k must be less than v")

  # Calculate r from lambda(v-1) = r(k-1)
  if ((lambda * (v - 1)) %% (k - 1) != 0) {
    warning("No BIBD exists: lambda(v-1) not divisible by (k-1)")
    return(NULL)
  }
  r <- lambda * (v - 1) %/% (k - 1)

  # Calculate b from bk = vr
  if ((v * r) %% k != 0) {
    warning("No BIBD exists: vr not divisible by k")
    return(NULL)
  }
  b <- (v * r) %/% k

  # Known BIBD constructions
  blocks <- .bibd_blocks(v, k, b, lambda)

  if (is.null(blocks)) {
    warning("BIBD construction not available for v=", v, ", k=", k)
    return(NULL)
  }

  # Build matrix
  rows_i <- integer(0)
  cols_j <- integer(0)

  for (i in seq_along(blocks)) {
    rows_i <- c(rows_i, rep(i, length(blocks[[i]])))
    cols_j <- c(cols_j, blocks[[i]])
  }

  result <- Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(b, v)
  )

  attr(result, "design") <- list(
    type = "bibd",
    v = v, b = b, r = r, k = k, lambda = lambda
  )

  result
}


# Internal: Known BIBD block constructions
.bibd_blocks <- function(v, k, b, lambda) {
  if (lambda != 1) return(NULL)

  # (7,3,1)-BIBD: Fano plane
  if (v == 7 && k == 3) {
    return(list(
      c(1,2,4), c(2,3,5), c(3,4,6), c(4,5,7),
      c(5,6,1), c(6,7,2), c(7,1,3)
    ))
  }

  # (9,3,1)-BIBD: Affine plane AG(2,3)
  if (v == 9 && k == 3) {
    return(list(
      c(1,2,3), c(4,5,6), c(7,8,9),
      c(1,4,7), c(2,5,8), c(3,6,9),
      c(1,5,9), c(2,6,7), c(3,4,8),
      c(1,6,8), c(2,4,9), c(3,5,7)
    ))
  }

  # (13,4,1)-BIBD
  if (v == 13 && k == 4) {
    return(list(
      c(1,2,4,10), c(2,3,5,11), c(3,4,6,12), c(4,5,7,13),
      c(5,6,8,1), c(6,7,9,2), c(7,8,10,3), c(8,9,11,4),
      c(9,10,12,5), c(10,11,13,6), c(11,12,1,7), c(12,13,2,8),
      c(13,1,3,9)
    ))
  }

  # (21,5,1)-BIBD
  if (v == 21 && k == 5) {
    # Projective plane PG(2,4) complement structure
    return(.pg_complement_blocks(4))
  }

  # (31,6,1)-BIBD
  if (v == 31 && k == 6) {
    return(.difference_set_bibd(31, 6))
  }

  NULL
}


# Internal: Generate BIBD from difference set
.difference_set_bibd <- function(v, k) {
  # For prime v, use Singer difference set if possible
  if (!is_prime(v)) return(NULL)

  # Find a (v,k,1) difference set
  # D is a set of k elements from Z_v such that every non-zero

  # element of Z_v can be expressed uniquely as d_i - d_j

  # For (31,6,1): D = {1, 5, 11, 24, 25, 27} (mod 31)
  if (v == 31 && k == 6) {
    D <- c(1, 5, 11, 24, 25, 27)
  } else {
    return(NULL)
  }

  # Generate all translates
  blocks <- lapply(0:(v-1), function(t) {
    ((D + t - 1) %% v) + 1
  })

  blocks
}


# Internal: PG complement blocks
.pg_complement_blocks <- function(q) {
  # For PG(2,q), the complement gives a BIBD
  n <- q^2 + q + 1  # points and lines
  k <- q^2  # complement block size

  # This is complex; return NULL for now
  NULL
}


# =============================================================================
# 4. HYPERCUBE DESIGNS
# =============================================================================

#' @title Generate A Hypercube Pooling Matrix
#'
#' @description
#' Construct the incidence matrix of a \eqn{d}-dimensional hypercube design by
#' arranging the individuals on a \eqn{q \times \cdots \times q} lattice and
#' pooling along coordinate hyperplanes.
#'
#' @param q Integer scalar giving the side length of each coordinate dimension.
#'   It must satisfy \eqn{q \ge 2}. The returned design contains \eqn{q} pools in
#'   each coordinate family.
#' @param d Integer scalar giving the number of dimensions. It must satisfy
#'   \eqn{d \ge 2}. Each individual belongs to exactly \eqn{d} pools, one from
#'   each coordinate family.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where \eqn{N = q^d} and \eqn{J = dq}. Each row has
#'   constant size \eqn{q^{d-1}} and each column has constant sum \eqn{d}. The
#'   attached \code{design} attribute is a named list with fields \code{type},
#'   \code{N}, \code{M}, \code{q}, \code{d}, \code{pool_size}, and
#'   \code{pools_per_sample}, where \code{M} is a legacy display label for the
#'   number of pools \eqn{J}.
#'
#' @details
#' Label each individual by a coordinate vector
#' \eqn{(x_1,\ldots,x_d) \in \{0,\ldots,q-1\}^d}. For each coordinate index
#' \eqn{\ell \in \{1,\ldots,d\}} and each value
#' \eqn{a \in \{0,\ldots,q-1\}}, define the pool
#' \deqn{\mathcal{P}_{\ell,a} = \{i : x_{\ell}(i) = a\}.}
#' The resulting matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} has
#' \deqn{J = dq}
#' rows and
#' \deqn{N = q^d}
#' columns. Every pool fixes one coordinate and lets the remaining
#' \eqn{d-1} coordinates vary freely, so
#' \deqn{n_j = |\mathcal{P}_j| = q^{d-1}.}
#' Every individual has one coordinate value in each dimension, hence
#' \deqn{w_i = |\mathcal{J}_i| = d.}
#' For \eqn{d = 2}, this recovers the square array design with equally sized row
#' and column families.
#'
#' @examples
#' M2 <- hypercube_matrix(q = 4, d = 2)
#' M2
#' Matrix::rowSums(M2)
#' Matrix::colSums(M2)
#'
#' M3 <- hypercube_matrix(q = 3, d = 3)
#' M3
#'
#' @references
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' Du, D.-Z., & Hwang, F. K. (2000). \emph{Combinatorial Group Testing and Its
#' Applications} (2nd ed.). World Scientific.
#'
#' @seealso
#' \code{\link{array_matrix}}, \code{\link{array_matrix_3d}},
#' \code{\link{compare_designs}}
#'
#' @family pooling_designs
#' @export
hypercube_matrix <- function(q, d) {
  q <- as.integer(q)
  d <- as.integer(d)

  if (q < 2) stop("q must be at least 2")
  if (d < 2) stop("d must be at least 2")

  N <- q^d
  M <- d * q

  # Generate all d-tuples (coordinates)
  coords <- as.matrix(expand.grid(rep(list(0:(q-1)), d)))

  rows_i <- integer(0)
  cols_j <- integer(0)

  # For each dimension, create q pools
  for (dim in 1:d) {
    for (val in 0:(q-1)) {
      pool_id <- (dim - 1) * q + val + 1

      # Samples where coordinate 'dim' equals 'val'
      sample_ids <- which(coords[, dim] == val)

      rows_i <- c(rows_i, rep(pool_id, length(sample_ids)))
      cols_j <- c(cols_j, sample_ids)
    }
  }

  result <- Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(M, N)
  )

  attr(result, "design") <- list(
    type = "hypercube",
    N = N,
    M = M,
    q = q,
    d = d,
    pool_size = q^(d - 1),
    pools_per_sample = d
  )

  result
}


# =============================================================================
# 5. KIRKMAN TRIPLE SYSTEM DESIGNS
# =============================================================================

#' @title Generate A Kirkman Triple System Pooling Matrix
#'
#' @description
#' Construct the incidence matrix of a Kirkman triple system (KTS): a
#' resolvable Steiner triple system for pooled testing, in which the triples
#' (blocks) can be partitioned into parallel classes ("days") that each cover
#' every individual exactly once.
#'
#' @param v Integer scalar giving the number of individuals. A Kirkman triple
#'   system exists if and only if \eqn{v \equiv 3 \pmod 6} and \eqn{v \geq 3}
#'   (Ray-Chaudhuri & Wilson, 1971); any other value is an error.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where \eqn{N = v} and \eqn{J = v(v-1)/6}. Each row has
#'   size three and each column sum is \eqn{(v-1)/2}. The attached
#'   \code{design} attribute is a named list with fields \code{type},
#'   \code{v}, \code{b}, \code{r}, \code{k}, and \code{parallel_class} (an
#'   integer vector of length \eqn{J} giving each row's parallel-class label,
#'   \eqn{1,\ldots,(v-1)/2}). Use \code{\link{is_kts}} to independently verify
#'   that a returned matrix satisfies all four defining KTS properties.
#'
#' @details
#' A Kirkman triple system is a Steiner triple system whose blocks can be
#' partitioned into parallel classes. In the returned incidence matrix
#' \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}, the columns index the \eqn{N=v}
#' individuals and the rows index the triples. The defining combinatorial
#' relations are
#' \deqn{n_j = |\mathcal{P}_j| = 3 \quad \text{for all } j,}
#' \deqn{w_i = |\mathcal{J}_i| = \frac{v-1}{2} \quad \text{for all } i,}
#' every unordered pair of distinct individuals appears together in exactly
#' one pool (so \eqn{J = v(v-1)/6}), and the rows can be partitioned into
#' \eqn{(v-1)/2} parallel classes of \eqn{v/3} rows each, every class
#' covering all \eqn{v} individuals exactly once. The necessary and
#' sufficient existence condition is
#' \deqn{v \equiv 3 \pmod 6, \quad v \geq 3.}
#'
#' # Construction
#'
#' Three dedicated, deterministic (non-search) constructions are used
#' whenever \eqn{v} matches their pattern:
#' \itemize{
#'   \item \eqn{v = 3^n} (\eqn{n = 1, 2, 3, 4, \ldots}, i.e. \eqn{v \in
#'     \{3, 9, 27, 81, 243, \ldots\}}): the affine geometry
#'     \eqn{\mathrm{AG}(n,3)} construction. Points are the \eqn{3^n} vectors
#'     of \eqn{\mathrm{GF}(3)^n}. For each of the \eqn{(3^n-1)/2} canonical
#'     line directions \eqn{d} (nonzero vectors of \eqn{\mathrm{GF}(3)^n}
#'     taken up to the equivalence \eqn{d \sim -d}), the lines \eqn{\{p, p+d,
#'     p+2d\}} for \eqn{p} ranging over \eqn{\mathrm{GF}(3)^n} partition all
#'     \eqn{3^n} points into \eqn{3^{n-1}} parallel lines, giving one
#'     parallel class per direction (Colbourn & Rosa, 1999, Section 3.3).
#'     This is the standard direct construction of a resolvable Steiner
#'     triple system from an affine space over \eqn{\mathrm{GF}(3)}; it
#'     covers \eqn{v = 3} (the trivial system) and \eqn{v = 9} as special
#'     cases, in addition to \eqn{v = 27} and \eqn{v = 81}.
#'   \item \eqn{v = 15}: the construction packs the 35 lines of
#'     \eqn{\mathrm{PG}(3,2)} into 7 pairwise line-disjoint spreads (Colbourn
#'     & Dinitz, 2007, "Kirkman systems"), the classical solution to
#'     Kirkman's 1850 schoolgirl problem.
#' }
#' For every other admissible \eqn{v} (e.g. \eqn{v = 21, 33, 39, \ldots}, and
#' also \eqn{v = 63}, since no deterministic construction is currently
#' implemented for it despite \eqn{63 = 2^6 - 1} suggesting a
#' \eqn{\mathrm{PG}(5,2)}-spread approach analogous to \eqn{v = 15} --
#' packing all 651 lines of \eqn{\mathrm{PG}(5,2)} into 31 disjoint spreads
#' is a substantially harder combinatorial search than the \eqn{v = 15} case
#' and was not solved within the current implementation), the constructor
#' falls back to a general randomized search: it generates a pool of
#' candidate parallel classes (uniformly random partitions of the \eqn{v}
#' points into triples) and performs an exact-cover backtracking search for
#' \eqn{(v-1)/2} of them that are pairwise pair-disjoint and jointly cover
#' every pair exactly once. Every constructed system, from every code path,
#' is independently re-verified with \code{\link{is_kts}} before being
#' returned; the function errors instead of ever returning an object that
#' fails that check.
#'
#' This search-based general fallback is empirically unreliable once
#' \eqn{v} moves beyond the dedicated cases above: it is not guaranteed to
#' find a system within its computational budget even for \eqn{v = 21}, and
#' is not tuned for larger \eqn{v} such as 33 or 63. This is a limitation of
#' the current search heuristic, not evidence that a KTS fails to exist --
#' by Ray-Chaudhuri & Wilson (1971), KTS(v) exists for every \eqn{v \equiv 3
#' \pmod 6}. If the search exhausts its attempt budget, \code{kirkman_matrix()}
#' stops with an informative error naming the values of \eqn{v} that are
#' guaranteed fast: any power of 3 (\eqn{3, 9, 27, 81, \ldots}), and 15.
#'
#' @examples
#' M9 <- kirkman_matrix(9)
#' M9
#' Matrix::rowSums(M9)
#' Matrix::colSums(M9)
#' is_kts(M9)
#'
#' # v = 15: the classic Kirkman schoolgirl problem
#' M15 <- kirkman_matrix(15)
#' is_kts(M15)
#' attr(M15, "design")$parallel_class
#'
#' # v = 27 and v = 81: AG(n,3), also fast and deterministic
#' is_kts(kirkman_matrix(27))
#' is_kts(kirkman_matrix(81))
#'
#' # Invalid v (not == 3 mod 6) is an error, not a silent NULL/warning
#' tryCatch(kirkman_matrix(10), error = function(e) conditionMessage(e))
#'
#' @references
#' Ray-Chaudhuri, D. K., & Wilson, R. M. (1971). Solution of Kirkman's
#' schoolgirl problem. \emph{Proceedings of Symposia in Pure Mathematics},
#' 19, 187-203.
#'
#' Colbourn, C. J., & Rosa, A. (1999). \emph{Triple Systems}. Oxford University
#' Press.
#'
#' Colbourn, C. J., & Dinitz, J. H. (Eds.). (2007). \emph{Handbook of
#' Combinatorial Designs} (2nd ed.). Chapman & Hall/CRC.
#'
#' @seealso
#' \code{\link{is_kts}} to verify the defining properties of a returned (or
#' any other) putative KTS incidence matrix; \code{\link{bibd_matrix}},
#' \code{\link{pg_matrix}}, \code{\link{compare_designs}}
#'
#' @family pooling_designs
#' @export
kirkman_matrix <- function(v) {
  v <- as.integer(v)

  if (is.na(v) || v < 3L) {
    stop("v must be an integer >= 3. Got v = ", v, ".", call. = FALSE)
  }
  if (v %% 6L != 3L) {
    stop("A Kirkman triple system does not exist for v = ", v, ": ",
         "v must satisfy v == 3 (mod 6) (e.g. 3, 9, 15, 21, 27, ...).",
         call. = FALSE)
  }

  b <- (v * (v - 1L)) %/% 6L
  classes <- .kts_construct(v)

  blocks <- unlist(classes, recursive = FALSE)
  rows_i <- integer(3L * length(blocks))
  cols_j <- integer(3L * length(blocks))
  parallel_class <- integer(length(blocks))
  row_id <- 0L
  for (cl in seq_along(classes)) {
    for (blk in classes[[cl]]) {
      row_id <- row_id + 1L
      idx <- seq.int(3L * (row_id - 1L) + 1L, length.out = 3L)
      rows_i[idx] <- row_id
      cols_j[idx] <- blk
      parallel_class[row_id] <- cl
    }
  }

  result <- Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(b, v)
  )

  attr(result, "design") <- list(
    type = "kirkman",
    v = v,
    b = b,
    r = (v - 1L) %/% 2L,
    k = 3L,
    parallel_class = parallel_class
  )

  if (!isTRUE(is_kts(result))) {
    stop("Internal error: kirkman_matrix(v = ", v, ") constructed a system ",
         "that failed is_kts() verification. This should not happen; ",
         "please report it.", call. = FALSE)
  }

  result
}


#' @title Verify the Defining Properties of a Kirkman Triple System
#'
#' @description
#' Independently checks whether a putative Kirkman triple system (KTS)
#' incidence matrix, such as one returned by \code{\link{kirkman_matrix}},
#' actually satisfies all four properties that define a KTS: block size
#' three, existence condition on \eqn{v}, every pair of points covered
#' exactly once, and resolvability into parallel classes that each partition
#' the point set.
#'
#' @param M A candidate KTS incidence matrix: a \eqn{J \times N} binary
#'   matrix (\code{matrix}, \code{dgCMatrix}, or other \code{Matrix} class)
#'   whose rows are blocks and whose columns are the \eqn{N} points. To check
#'   resolvability, \code{M} should carry a \code{design} attribute with a
#'   \code{parallel_class} component (a length-\eqn{J} vector of class
#'   labels), exactly as attached by \code{\link{kirkman_matrix}}; without
#'   it, resolvability is reported as a failure rather than skipped, since an
#'   unresolved block set is not a Kirkman triple system.
#'
#' @return A single \code{logical} scalar, \code{TRUE} only if every one of
#'   the four properties holds. On \code{FALSE}, the result carries a
#'   character-vector attribute \code{"reasons"} describing every violation
#'   found; pass \code{verbose = TRUE} to also print them.
#' @param verbose Logical; if \code{TRUE} and the check fails, print the
#'   failure reasons via \code{message()}. Default \code{FALSE}.
#'
#' @details
#' Writing \eqn{v = N} for the number of points and \eqn{b = J} for the
#' number of blocks, this function checks, directly against the incidence
#' matrix (not against any internal bookkeeping the constructor may have
#' used):
#' \enumerate{
#'   \item \eqn{v \equiv 3 \pmod 6} and \eqn{b = v(v-1)/6};
#'   \item every row has exactly 3 nonzero entries;
#'   \item every unordered pair of columns co-occurs (both 1) in exactly one
#'     row;
#'   \item the rows, grouped by \code{attr(M, "design")$parallel_class},
#'     form exactly \eqn{(v-1)/2} classes of \eqn{v/3} rows each, and within
#'     each class the rows partition \eqn{\{1,\ldots,v\}} (every point
#'     covered exactly once, none missing, none repeated).
#' }
#'
#' @examples
#' M15 <- kirkman_matrix(15)
#' is_kts(M15)
#'
#' # A deliberately broken KTS(15) (the pre-fix package bug): a hand-built
#' # block set whose last three parallel classes repeat some points and
#' # omit others, so several pairs are covered zero or more than one time.
#' broken_blocks <- list(
#'   c(1,2,3), c(4,5,6), c(7,8,9), c(10,11,12), c(13,14,15),
#'   c(1,4,7), c(2,5,8), c(3,6,9), c(10,13,14), c(11,12,15),
#'   c(1,5,9), c(2,6,7), c(3,4,8), c(10,12,14), c(11,13,15),
#'   c(1,6,8), c(2,4,9), c(3,5,7), c(10,11,14), c(12,13,15),
#'   c(1,10,15), c(2,11,13), c(3,12,14), c(4,8,13), c(5,7,12),
#'   c(1,11,14), c(2,10,12), c(3,13,15), c(4,7,15), c(5,8,10),
#'   c(1,12,13), c(2,14,15), c(3,10,11), c(4,9,14), c(5,6,13)
#' )
#' rows_i <- unlist(lapply(seq_along(broken_blocks), function(i) rep(i, 3)))
#' cols_j <- unlist(broken_blocks)
#' M_broken <- Matrix::sparseMatrix(
#'   i = rows_i, j = cols_j, x = rep(1L, length(rows_i)), dims = c(35, 15)
#' )
#' attr(M_broken, "design") <- list(parallel_class = rep(1:7, each = 5))
#' is_kts(M_broken)
#' attr(is_kts(M_broken), "reasons")
#'
#' @references
#' Ray-Chaudhuri, D. K., & Wilson, R. M. (1971). Solution of Kirkman's
#' schoolgirl problem. \emph{Proceedings of Symposia in Pure Mathematics},
#' 19, 187-203.
#'
#' @seealso \code{\link{kirkman_matrix}}
#' @family pooling_designs
#' @export
is_kts <- function(M, verbose = FALSE) {
  reasons <- character(0)
  fail <- function(msg) reasons <<- c(reasons, msg)

  if (is.null(M) || (!is.matrix(M) && !inherits(M, "Matrix"))) {
    fail("M is not a matrix")
    return(.kts_verdict(FALSE, reasons, verbose))
  }

  v <- ncol(M)
  b <- nrow(M)

  if (v < 3L || v %% 6L != 3L) {
    fail(sprintf(
      "ncol(M) = %d is not a valid KTS order (must be >= 3 and == 3 mod 6)",
      v
    ))
    return(.kts_verdict(FALSE, reasons, verbose))
  }

  b_expected <- (v * (v - 1L)) %/% 6L
  if (b != b_expected) {
    fail(sprintf("nrow(M) = %d, expected v(v-1)/6 = %d", b, b_expected))
    return(.kts_verdict(FALSE, reasons, verbose))
  }

  Mi <- as.matrix(M) != 0
  row_sizes <- rowSums(Mi)
  if (!all(row_sizes == 3L)) {
    fail(sprintf(
      "%d block(s) do not have exactly 3 points",
      sum(row_sizes != 3L)
    ))
  }

  pair_count <- matrix(0L, v, v)
  for (i in seq_len(b)) {
    pts <- which(Mi[i, ])
    if (length(pts) == 3L) {
      pr <- utils::combn(pts, 2)
      for (k in seq_len(ncol(pr))) {
        a <- pr[1L, k]; c2 <- pr[2L, k]
        pair_count[a, c2] <- pair_count[a, c2] + 1L
        pair_count[c2, a] <- pair_count[c2, a] + 1L
      }
    }
  }
  upper <- pair_count[upper.tri(pair_count)]
  if (!all(upper == 1L)) {
    fail(sprintf(
      "%d pair(s) of points are covered 0 or >1 times (each pair must occur in exactly one block); %d pair(s) never covered, %d pair(s) covered more than once",
      sum(upper != 1L), sum(upper == 0L), sum(upper > 1L)
    ))
  }

  design <- attr(M, "design")
  pc <- design[["parallel_class"]]
  if (is.null(pc)) {
    fail("M has no design$parallel_class attribute; resolvability cannot be verified, so M is not confirmed to be a Kirkman (resolvable) triple system")
  } else if (length(pc) != b) {
    fail(sprintf(
      "design$parallel_class has length %d, expected nrow(M) = %d",
      length(pc), b
    ))
  } else {
    expected_classes <- (v - 1L) %/% 2L
    expected_block_per_class <- v %/% 3L
    classes <- split(seq_len(b), pc)
    if (length(classes) != expected_classes) {
      fail(sprintf(
        "found %d parallel class(es), expected (v-1)/2 = %d",
        length(classes), expected_classes
      ))
    }
    for (cl_name in names(classes)) {
      idx <- classes[[cl_name]]
      if (length(idx) != expected_block_per_class) {
        fail(sprintf(
          "parallel class '%s' has %d block(s), expected v/3 = %d",
          cl_name, length(idx), expected_block_per_class
        ))
        next
      }
      pts <- unlist(lapply(idx, function(i) which(Mi[i, ])))
      if (length(pts) != v || anyDuplicated(pts) || !setequal(pts, seq_len(v))) {
        fail(sprintf(
          "parallel class '%s' does not partition the %d points exactly once",
          cl_name, v
        ))
      }
    }
  }

  .kts_verdict(length(reasons) == 0L, reasons, verbose)
}


# Internal: attach the "reasons" attribute and optionally message() them.
.kts_verdict <- function(ok, reasons, verbose) {
  if (verbose && !ok) {
    message(paste(reasons, collapse = "\n"))
  }
  structure(ok, reasons = reasons)
}


# Internal: dispatch to a fast dedicated construction when one is known for
# v (v = 3^n via AG(n,3), or v = 15 via the PG(3,2) spread packing), or to
# the general search-based constructor otherwise. Always returns a list of
# parallel classes, each a list of length-3 integer vectors; errors instead
# of returning an unverifiable or invalid result.
.kts_construct <- function(v) {
  n3 <- .kts_pow3_exponent(v)
  if (!is.na(n3)) {
    return(.kts_ag3(n3))
  }
  if (v == 15L) {
    return(.kts15_pg32())
  }

  classes <- .kts_general_search(v)
  if (is.null(classes)) {
    stop(
      "kirkman_matrix(): the general search-based constructor could not ",
      "find a Kirkman triple system for v = ", v, " within its ",
      "computational budget. This is a limitation of the current search, ",
      "not evidence that KTS(", v, ") does not exist (it does, for every ",
      "v == 3 mod 6, by Ray-Chaudhuri & Wilson 1971). Values of v currently ",
      "guaranteed to construct quickly: 3, 9, 15, 27, 81 (via AG(n,3)), ",
      "and any other power of 3.",
      call. = FALSE
    )
  }
  classes
}


# Internal: if v is an exact power of 3 (v = 3^n, n >= 1), return n;
# otherwise return NA_integer_. Used to dispatch v = 3, 9, 27, 81, ... to
# the AG(n,3) construction below.
.kts_pow3_exponent <- function(v) {
  if (v < 3L) return(NA_integer_)
  n <- 0L
  x <- v
  while (x %% 3L == 0L) {
    x <- x %/% 3L
    n <- n + 1L
  }
  if (x == 1L) n else NA_integer_
}


# Internal: deterministic construction of a Kirkman triple system for
# v = 3^n points via the affine geometry AG(n,3) over GF(3). Points are the
# 3^n vectors of GF(3)^n (encoded here as base-3 integers 0..3^n-1, then
# shifted to the package's 1-indexed point labels). For each of the
# (3^n-1)/2 canonical line "directions" d (nonzero vectors of GF(3)^n, taken
# up to the equivalence d ~ -d = 2d since GF(3)^* = {1,2}, with the
# representative normalized so its first nonzero coordinate is 1), the
# lines {p, p+d, p+2d} (addition mod 3, componentwise) for p ranging over
# GF(3)^n partition all 3^n points into 3^(n-1) parallel lines -- one
# parallel class per direction. This is the standard direct (non-search)
# construction of a resolvable Steiner triple system from an affine space
# over GF(3) (Colbourn & Rosa, 1999, Triple Systems, Section 3.3), used here
# for v = 3, 9, 27, 81 (n = 1, 2, 3, 4); kirkman_matrix() re-verifies the
# result with is_kts() regardless. Returns a list of (3^n-1)/2 parallel
# classes, each a list of 3^(n-1) length-3 integer vectors of 1-indexed
# point labels.
.kts_ag3 <- function(n) {
  v <- 3L^n
  pow3 <- 3L^(0:(n - 1L))

  # pts[id + 1, ] is the base-3 digit vector (GF(3)^n coordinates) of the
  # 0-indexed point `id`.
  pts <- matrix(0L, v, n)
  rem <- 0:(v - 1L)
  for (k in seq_len(n)) {
    pts[, k] <- rem %% 3L
    rem <- rem %/% 3L
  }

  # Canonical direction representatives: nonzero vectors whose first
  # nonzero coordinate is 1, i.e. one representative per {d, 2d} pair.
  dir_ids <- integer(0)
  for (id in seq_len(v - 1L)) {
    vec <- pts[id + 1L, ]
    first_nz <- vec[vec != 0L][1L]
    if (first_nz == 1L) dir_ids <- c(dir_ids, id)
  }
  stopifnot(length(dir_ids) == (v - 1L) %/% 2L)

  classes <- vector("list", length(dir_ids))
  for (ci in seq_along(dir_ids)) {
    d <- pts[dir_ids[ci] + 1L, ]
    assigned <- logical(v)
    cls <- vector("list", v %/% 3L)
    cnt <- 0L
    for (id in 0:(v - 1L)) {
      if (!assigned[id + 1L]) {
        p <- pts[id + 1L, ]
        p1 <- (p + d) %% 3L
        p2 <- (p + 2L * d) %% 3L
        id1 <- sum(p1 * pow3)
        id2 <- sum(p2 * pow3)
        assigned[c(id + 1L, id1 + 1L, id2 + 1L)] <- TRUE
        cnt <- cnt + 1L
        cls[[cnt]] <- sort(c(id + 1L, id1 + 1L, id2 + 1L))
      }
    }
    classes[[ci]] <- cls
  }
  classes
}


# Internal: pairwise key for an unordered pair of points, vectorized.
.kts_pair_key <- function(i, j) pmin(i, j) * 100000L + pmax(i, j)

# Internal: a uniformly random partition of 1:v into v/3 triples.
.kts_random_class <- function(v) {
  pts <- sample.int(v)
  m <- v %/% 3L
  lapply(seq_len(m), function(k) {
    sort(pts[(3L * (k - 1L) + 1L):(3L * k)])
  })
}

# Internal: sorted vector of pair-keys covered by a parallel class.
.kts_class_pairs <- function(cls) {
  unlist(lapply(cls, function(b) {
    pr <- utils::combn(b, 2)
    .kts_pair_key(pr[1L, ], pr[2L, ])
  }))
}

# Internal: general KTS constructor for arbitrary v == 3 (mod 6), v != 15.
# Two-phase randomized construction: (1) generate a pool of candidate
# parallel classes (random partitions of the v points into triples -- always
# trivially valid on their own, no pairing constraints involved), then (2)
# an MRV-ordered exact-cover backtracking search picks (v-1)/2 of them that
# are pairwise pair-disjoint and jointly cover every pair exactly once. Sound
# by construction (any result returned satisfies the KTS properties, which
# kirkman_matrix() re-verifies with is_kts() regardless); not guaranteed to
# be complete within the bounded attempt/node budget, in which case it
# returns NULL rather than searching indefinitely.
.kts_general_search <- function(v, n_candidates = 2500L, max_nodes = 150000L,
                                 max_attempts = 5L, seed = 1L) {
  n_classes <- (v - 1L) %/% 2L

  attempt_once <- function() {
    seen <- new.env(parent = emptyenv())
    cand_classes <- vector("list", 0L)
    cand_pairs <- vector("list", 0L)
    stale <- 0L
    while (length(cand_classes) < n_candidates && stale < 500L) {
      cls <- .kts_random_class(v)
      pk <- sort(.kts_class_pairs(cls))
      key <- paste(pk, collapse = ",")
      if (is.null(seen[[key]])) {
        seen[[key]] <- TRUE
        cand_classes[[length(cand_classes) + 1L]] <- cls
        cand_pairs[[length(cand_pairs) + 1L]] <- pk
        stale <- 0L
      } else {
        stale <- stale + 1L
      }
    }
    n_cand <- length(cand_classes)

    all_pairs <- integer(0)
    if (v >= 2L) {
      for (i in seq_len(v - 1L)) {
        all_pairs <- c(all_pairs, .kts_pair_key(i, (i + 1L):v))
      }
    }
    all_pairs <- sort(unique(all_pairs))
    n_pairs_total <- length(all_pairs)
    pair_index <- stats::setNames(seq_along(all_pairs), all_pairs)
    cand_pair_idx <- lapply(cand_pairs, function(pk) {
      unname(pair_index[as.character(pk)])
    })

    cand_by_pair <- vector("list", n_pairs_total)
    for (ci in seq_len(n_cand)) {
      for (pidx in cand_pair_idx[[ci]]) {
        cand_by_pair[[pidx]] <- c(cand_by_pair[[pidx]], ci)
      }
    }

    covered <- logical(n_pairs_total)
    used_cand <- logical(n_cand)
    chosen <- integer(n_classes)
    nodes <- 0L

    rec <- function(depth) {
      nodes <<- nodes + 1L
      if (nodes > max_nodes) return(NA)
      if (depth > n_classes) return(TRUE)
      unc <- which(!covered)
      best_opts <- NULL
      best_n <- Inf
      for (pidx in unc) {
        opts <- cand_by_pair[[pidx]]
        opts <- opts[!used_cand[opts]]
        if (length(opts) < best_n) {
          best_n <- length(opts)
          best_opts <- opts
          if (best_n <= 1L) break
        }
      }
      if (identical(best_n, 0)) return(FALSE)
      for (ci in best_opts) {
        pidxs <- cand_pair_idx[[ci]]
        if (any(covered[pidxs])) next
        covered[pidxs] <<- TRUE
        used_cand[ci] <<- TRUE
        chosen[depth] <<- ci
        res <- rec(depth + 1L)
        if (isTRUE(res)) return(TRUE)
        covered[pidxs] <<- FALSE
        used_cand[ci] <<- FALSE
        if (is.na(res)) return(NA)
      }
      FALSE
    }

    res <- rec(1L)
    if (!isTRUE(res)) return(NULL)
    lapply(chosen, function(ci) cand_classes[[ci]])
  }

  .kts_with_fixed_seed(seed, {
    for (attempt in seq_len(max_attempts)) {
      set.seed(seed * 1000003L + attempt)
      classes <- attempt_once()
      if (!is.null(classes)) return(classes)
    }
    NULL
  })
}

# Internal: run `expr` under a fixed RNG seed, saving and restoring the
# caller's global RNG state so this function has no visible side effect on
# the user's random-number stream.
.kts_with_fixed_seed <- function(seed, expr) {
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

# Internal: dedicated construction for v = 15, the classic Kirkman schoolgirl
# problem. Packs the 35 lines of PG(3,2) (points 1..15 as nonzero vectors of
# GF(2)^4, lines {a, b, a xor b}) into 7 pairwise line-disjoint spreads (each
# a set of 5 pairwise point-disjoint lines covering all 15 points). This is a
# standard classical construction (Colbourn & Dinitz, 2007) and is fully
# deterministic. Returns a list of 7 parallel classes, each a list of 5
# length-3 integer vectors.
.kts15_pg32 <- function() {
  seen <- character(0)
  lines <- list()
  for (a in 1:14) {
    for (b in (a + 1L):15) {
      cc <- bitwXor(a, b)
      if (cc == a || cc == b) next
      trip <- sort(c(a, b, cc))
      key <- paste(trip, collapse = ",")
      if (!(key %in% seen)) {
        seen <- c(seen, key)
        lines[[length(lines) + 1L]] <- trip
      }
    }
  }
  ord <- order(
    vapply(lines, `[`, integer(1L), 1L),
    vapply(lines, `[`, integer(1L), 2L),
    vapply(lines, `[`, integer(1L), 3L)
  )
  lines <- lines[ord]
  stopifnot(length(lines) == 35L)

  spreads <- list()
  rec_spread <- function(cover, chosen) {
    if (length(cover) == 15L) {
      spreads[[length(spreads) + 1L]] <<- chosen
      return(invisible())
    }
    p <- min(setdiff(1:15, cover))
    for (i in seq_along(lines)) {
      l <- lines[[i]]
      if (p %in% l && length(intersect(l, cover)) == 0L) {
        rec_spread(union(cover, l), c(chosen, i))
      }
    }
  }
  rec_spread(integer(0), integer(0))

  sol <- NULL
  pack <- function(used, chosen) {
    if (!is.null(sol)) return(invisible())
    if (length(used) == 35L) { sol <<- chosen; return(invisible()) }
    first <- min(setdiff(1:35, used))
    for (sp in spreads) {
      if (!is.null(sol)) return(invisible())
      if (first %in% sp && length(intersect(sp, used)) == 0L) {
        pack(union(used, sp), c(chosen, list(sp)))
      }
    }
  }
  pack(integer(0), list())
  stopifnot(!is.null(sol), length(sol) == 7L)

  lapply(sol, function(sp) lapply(sp, function(i) lines[[i]]))
}


# =============================================================================
# 6. PROJECTIVE GEOMETRY DESIGNS
# =============================================================================

#' @title Generate A Projective Plane Pooling Matrix
#'
#' @description
#' Construct the incidence matrix of the projective plane \eqn{\mathrm{PG}(2,q)}
#' over the finite field \eqn{\mathrm{GF}(q)}, where points act as individuals
#' and lines act as pools.
#'
#' @param q Integer scalar giving the field order. It must be a prime power so
#'   that \eqn{\mathrm{GF}(q)} exists. Supported values are those handled by
#'   \code{\link{gf}} and \code{\link{is_prime_power}}.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where
#'   \eqn{N = J = q^2 + q + 1}. Every row sum and every column sum equals
#'   \eqn{q+1}. The attached \code{design} attribute is a named list with fields
#'   \code{type}, \code{q}, \code{N}, \code{M}, \code{points_per_line}, and
#'   \code{lines_per_point}.
#'
#' @details
#' The projective plane \eqn{\mathrm{PG}(2,q)} consists of the one-dimensional
#' subspaces of \eqn{\mathrm{GF}(q)^3}. This implementation uses canonical
#' homogeneous representatives of the form
#' \deqn{(1, a, b), \qquad (0, 1, c), \qquad (0, 0, 1),}
#' where \eqn{a,b,c \in \mathrm{GF}(q)}. These representatives index both the
#' \eqn{N = q^2 + q + 1} points and, by duality, the \eqn{J = q^2 + q + 1}
#' lines. A point \eqn{P = [x:y:z]} is incident with a line
#' \eqn{L = [\alpha:\beta:\gamma]} precisely when
#' \deqn{\alpha x + \beta y + \gamma z = 0 \quad \text{in } \mathrm{GF}(q).}
#' The returned incidence matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}
#' therefore satisfies
#' \deqn{M_{ji} = 1 \iff P_i \text{ lies on line } j.}
#' The standard projective-plane identities become
#' \deqn{n_j = |\mathcal{P}_j| = q + 1, \qquad
#' w_i = |\mathcal{J}_i| = q + 1,}
#' and every pair of distinct individuals belongs to exactly one common pool.
#'
#' @examples
#' M2 <- pg_matrix(2)
#' M2
#' Matrix::rowSums(M2)
#' Matrix::colSums(M2)
#'
#' M3 <- pg_matrix(3)
#' M3
#'
#' @references
#' Hirschfeld, J. W. P. (1998). \emph{Projective Geometries over Finite Fields}
#' (2nd ed.). Oxford University Press.
#'
#' Colbourn, C. J., & Dinitz, J. H. (Eds.). (2007). \emph{Handbook of
#' Combinatorial Designs} (2nd ed.). Chapman & Hall/CRC.
#'
#' @seealso
#' \code{\link{gf}}, \code{\link{is_prime_power}}, \code{\link{bibd_matrix}}
#'
#' @family pooling_designs
#' @export
pg_matrix <- function(q) {
  pp <- is_prime_power(q)
  if (!pp$is_prime_power) {
    stop("q must be a prime power. Got q =", q)
  }

  q <- as.integer(q)
  n <- q^2 + q + 1  # Number of points = number of lines

  gf_obj <- gf(q)

  # Generate all points in homogeneous coordinates
  # Points: equivalence classes [x:y:z] with (x,y,z) != (0,0,0)
  # Canonical forms: (1,a,b), (0,1,c), (0,0,1)

  points <- list()
  idx <- 1

  # Type (1, a, b)
  for (a in 0:(q-1)) {
    for (b in 0:(q-1)) {
      points[[idx]] <- c(1L, a, b)
      idx <- idx + 1
    }
  }

  # Type (0, 1, c)
  for (c in 0:(q-1)) {
    points[[idx]] <- c(0L, 1L, c)
    idx <- idx + 1
  }

  # Type (0, 0, 1)
  points[[idx]] <- c(0L, 0L, 1L)

  # Lines have the same representation (duality)
  lines <- points

  # Incidence: point P is on line L iff P * L = 0 in GF(q)
  rows_i <- integer(0)
  cols_j <- integer(0)

  for (i in 1:n) {
    L <- lines[[i]]
    for (j in 1:n) {
      P <- points[[j]]

      # Compute dot product in GF(q)
      dot <- 0L
      for (k in 1:3) {
        prod <- gf_mult(P[k], L[k], gf_obj)
        dot <- gf_add(dot, prod, gf_obj)
      }

      if (dot == 0L) {
        rows_i <- c(rows_i, i)
        cols_j <- c(cols_j, j)
      }
    }
  }

  result <- Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(n, n)
  )

  attr(result, "design") <- list(
    type = "projective_plane",
    q = q,
    N = n,
    M = n,
    points_per_line = q + 1L,
    lines_per_point = q + 1L
  )

  result
}


# =============================================================================
# 7. P-BEST CLINICAL DESIGN (Zismanov et al. 2024)
# =============================================================================

#' @title Generate The Clinical P-BEST 2024 Pooling Matrix
#'
#' @description
#' Construct the pseudo-random balanced \eqn{J \times N} pooling matrix used for
#' the 384-to-94 clinical P-BEST design reported by Zismanov et al. (2024).
#'
#' @param N Integer scalar giving the total number of individuals. The default
#'   is \code{384L}. In the returned design this is the column dimension
#'   \eqn{N} of \eqn{\mathbf{M}}.
#' @param J Integer scalar giving the total number of pools. The default is
#'   \code{94L}. In the returned design this is the row dimension \eqn{J} of
#'   \eqn{\mathbf{M}}.
#' @param seed Integer scalar used to initialize the pseudo-random tie-breaking
#'   step. The default is \code{42L}. Different seeds generate different balanced
#'   realizations of the same construction heuristic.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}. Under the default \eqn{(N,J) = (384,94)} setting, each
#'   column sum is seven and each row sum lies in \eqn{\{28,29\}}. The attached
#'   \code{design} attribute is a named list with fields \code{type},
#'   \code{N}, \code{J}, \code{seed}, \code{pools_per_sample}, and
#'   \code{pool_size_range}.
#'
#' @details
#' This function targets the clinical design with
#' \deqn{N = 384, \qquad J = 94.}
#' The implementation fixes the number of pools per individual at
#' \deqn{w_i = |\mathcal{J}_i| = r = 7}
#' and therefore creates
#' \deqn{T = Nr}
#' total incidences. For the default parameter values,
#' \deqn{T = 384 \times 7 = 2688.}
#' Distributing these incidences over \eqn{J=94} pools implies an average pool
#' size of \eqn{T/J \approx 28.6}. The construction therefore aims for a nearly
#' balanced row-load profile with
#' \deqn{56 \text{ pools of size } 29 \quad \text{and} \quad 38 \text{ pools of
#' size } 28.}
#' The algorithm processes the individuals sequentially and assigns each
#' individual to the currently least-loaded pools, using pseudo-random jitter to
#' break ties. The output is thus a balanced heuristic realization rather than a
#' closed-form algebraic design. It is documented as a package-specific
#' heuristic constructor motivated by the reported clinical matrix dimensions,
#' not as a claim that the exact clinical deployment matrix has been recovered.
#'
#' @examples
#' M <- pbest_clinical_matrix()
#' M
#' range(Matrix::rowSums(M))
#' range(Matrix::colSums(M))
#'
#' M2 <- pbest_clinical_matrix(seed = 42L)
#' M2
#' identical(as.matrix(M), as.matrix(M2))
#'
#' @references
#' Zismanov, S., Yelin, I., Klochendler, A., et al. (2024). High capacity
#' clinical SARS-CoV-2 molecular testing using combinatorial pooling.
#' \emph{Communications Medicine}, 4, Article 121.
#' \doi{10.1038/s43856-024-00531-w}.
#'
#' @seealso
#' \code{\link{pp_matrix}}, \code{\link{compare_designs}}, \code{\link{list_designs}}
#'
#' @family pooling_designs
#' @export
pbest_clinical_matrix <- function(N = 384L, J = 94L, seed = 42L) {
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
  set.seed(seed)
  N <- as.integer(N)
  J <- as.integer(J)

  if (N < 1L) stop("N must be a positive integer")
  if (J < 1L) stop("J must be a positive integer")

  # Target pools per sample: r=7 gives T=2688 -> 28-29 samples/pool
  # (Zismanov et al. 2024, Table 1: 384->94 design)
  r <- 7L
  if (r > J) stop("pools per sample (7) cannot exceed J")

  pool_load   <- integer(J)
  assignments <- vector("list", N)

  # Greedy balanced assignment: each sample picks the r most under-loaded pools.
  # Uniform random jitter breaks ties and produces pseudo-random variation.
  for (s in seq_len(N)) {
    score  <- -pool_load + runif(J) * 0.5
    chosen <- order(score, decreasing = TRUE)[seq_len(r)]
    assignments[[s]] <- sort(chosen)
    pool_load[chosen] <- pool_load[chosen] + 1L
  }

  rows_i <- unlist(assignments)
  cols_j <- rep(seq_len(N), each = r)

  result <- Matrix::sparseMatrix(
    i    = rows_i,
    j    = cols_j,
    x    = rep(1L, length(rows_i)),
    dims = c(J, N)
  )

  attr(result, "design") <- list(
    type             = "pbest_clinical_2024",
    N                = N,
    J                = J,
    seed             = seed,
    pools_per_sample = r,
    pool_size_range  = range(Matrix::rowSums(result))
  )

  result
}


# =============================================================================
# DESIGN COMPARISON AND SELECTION
# =============================================================================

#' @title Compare Candidate Pooling Designs
#'
#' @description
#' Summarize several built-in pooling design families for a target number of
#' individuals by reporting the achieved matrix dimensions and simple structural
#' metrics.
#'
#' @param N Integer scalar giving the target number of individuals. This is the
#'   desired column dimension \eqn{N} for the candidate design matrices.
#' @param designs Character vector selecting which design families to evaluate.
#'   Allowed values are drawn from \code{"dorfman"}, \code{"array"},
#'   \code{"hypercube"}, and \code{"pp"}. The default evaluates all four.
#'
#' @return A data frame with one row per retained design candidate. The columns
#'   are \code{type}, \code{N}, \code{M}, \code{pool_size},
#'   \code{pools_per_sample}, and \code{compression}. Here \code{N} is the
#'   achieved number of individuals, \code{M} is a legacy display label for the
#'   number of pools \eqn{J}, and
#'   \code{compression} is the ratio \eqn{N / J}. The function returns
#'   \code{NULL} if no candidate design is available.
#'
#' @details
#' For each selected design family, this function computes a candidate pooling
#' matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} or a parameter proxy and then
#' reports the basic design dimensions. The reported compression metric is
#' \deqn{\mathrm{compression} = \frac{N}{J},}
#' where \eqn{N} is the achieved number of individuals and \eqn{J} is the number
#' of pools in the candidate design. Some families, such as array and hypercube
#' designs, may adjust the achieved \eqn{N} upward from the target value in order
#' to satisfy structural constraints like \eqn{N = rc} or \eqn{N = q^d}. The
#' summary is therefore comparative rather than an exact design optimization
#' routine.
#'
#' @examples
#' comparison_100 <- compare_designs(100)
#' comparison_100
#' comparison_384 <- compare_designs(384, designs = c("dorfman", "pp"))
#' comparison_384
#'
#' @references
#' Dorfman, R. (1943). The detection of defective members of large populations.
#' \emph{The Annals of Mathematical Statistics}, 14(4), 436-440.
#'
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' @seealso
#' \code{\link{list_designs}}, \code{\link{dorfman_matrix}},
#' \code{\link{hypercube_matrix}}, \code{\link{pp_design}}
#'
#' @family pooling_designs
#' @export
compare_designs <- function(N, designs = c("dorfman", "array", "hypercube", "pp")) {
  N <- as.integer(N)
  results <- list()

  if ("dorfman" %in% designs) {
    g_opt <- dorfman_optimal_g(0.01)
    M_dorf <- ceiling(N / g_opt)
    results$dorfman <- data.frame(
      type = "dorfman",
      N = N,
      M = M_dorf,
      pool_size = g_opt,
      pools_per_sample = 1,
      compression = round(N / M_dorf, 2)
    )
  }

  if ("array" %in% designs) {
    r <- ceiling(sqrt(N))
    c <- ceiling(N / r)
    N_actual <- r * c
    M_arr <- r + c
    results$array <- data.frame(
      type = "array_2d",
      N = N_actual,
      M = M_arr,
      pool_size = paste0(c, "/", r),
      pools_per_sample = 2,
      compression = round(N_actual / M_arr, 2)
    )
  }

  if ("hypercube" %in% designs && N >= 8) {
    # Find best (q, d)
    best <- NULL
    for (d in 2:4) {
      q <- ceiling(N^(1/d))
      N_actual <- q^d
      if (N_actual >= N && N_actual <= 2*N) {
        M_hc <- d * q
        comp <- N_actual / M_hc
        if (is.null(best) || comp > best$compression) {
          best <- data.frame(
            type = "hypercube",
            N = N_actual,
            M = M_hc,
            pool_size = q^(d-1),
            pools_per_sample = d,
            compression = round(comp, 2)
          )
        }
      }
    }
    if (!is.null(best)) results$hypercube <- best
  }

  if ("pp" %in% designs) {
    # Find best PP configuration
    configs <- pp_design(N_min = N * 0.8, N_max = N * 1.5)
    if (!is.null(configs) && nrow(configs) > 0) {
      # Pick closest to N with k_max >= 2
      configs <- configs[configs$k_max >= 2, ]
      if (nrow(configs) > 0) {
        idx <- which.min(abs(configs$N - N))
        cfg <- configs[idx, ]
        results$pp <- data.frame(
          type = paste0("pp(", cfg$q, ",", cfg$d, ",", cfg$nl, ")"),
          N = cfg$N,
          M = cfg$M,
          pool_size = cfg$pool_size,
          pools_per_sample = cfg$nl,
          compression = round(cfg$compression, 2)
        )
      }
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}


#' @title List Available Pooling Design Families
#'
#' @description
#' Return a compact catalog of the pooling design families exposed by this
#' package layer, together with their structural interpretations.
#'
#' @return A data frame with columns \code{type}, \code{function_name},
#'   \code{description}, \code{pools_per_sample}, and
#'   \code{detection_guarantee}. Each row describes one design family and its
#'   associated constructor or parameterization.
#'
#' @details
#' The returned table is a descriptive index rather than a design matrix. Its
#' entries summarize how each family maps individuals to pools in a binary
#' matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}. In particular, the
#' \code{pools_per_sample} column records the intended column weight
#' \eqn{w_i = |\mathcal{J}_i|}, while \code{detection_guarantee} records a brief
#' family-specific statement about the positive-pattern resolution guaranteed by
#' the corresponding construction.
#'
#' @examples
#' design_catalog <- list_designs()
#' design_catalog
#'
#' @references
#' Du, D.-Z., & Hwang, F. K. (2000). \emph{Combinatorial Group Testing and Its
#' Applications} (2nd ed.). World Scientific.
#'
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' @seealso
#' \code{\link{compare_designs}}, \code{\link{dorfman_matrix}},
#' \code{\link{pp_matrix}}
#'
#' @family pooling_designs
#' @export
list_designs <- function() {
  data.frame(
    type = c("dorfman", "array_2d", "array_3d", "bibd", "hypercube",
             "kirkman", "pg", "pp", "pbest"),
    function_name = c("dorfman_matrix", "array_matrix", "array_matrix_3d",
                      "bibd_matrix", "hypercube_matrix", "kirkman_matrix",
                      "pg_matrix", "pp_matrix", "pp_matrix(q=8,d=3,nl=6,N=384)"),
    description = c(
      "Dorfman two-stage pooling (non-overlapping)",
      "2D array (row/column) design",
      "3D array design",
      "Balanced Incomplete Block Design",
      "d-dimensional hypercube",
      "Kirkman Triple System",
      "Projective Geometry PG(2,q)",
      "Polynomial Pools algorithm",
      "P-BEST: 384 samples in 48 pools"
    ),
    pools_per_sample = c(1, 2, 3, "r", "d", "(v-1)/2", "q+1", "nl", 6),
    detection_guarantee = c("1", "1", "1", "depends", "1", "1",
                            "depends", "floor((nl-1)/(d-1))", "2"),
    stringsAsFactors = FALSE
  )
}
