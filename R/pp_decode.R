#' @title Decoding Algorithms for Pooled Test Results
#'
#' @description
#' Algorithms for decoding pooled test results: COMP (Combinatorial Orthogonal
#' Matching Pursuit) for zero-false-negative detection, and GPSR (Gradient
#' Projection for Sparse Reconstruction) for candidate ranking.
#'
#' @details
#' In the package-wide notation, observed pool outcomes are denoted by
#' \eqn{\mathbf{z} \in \{0,1\}^J}. The decoding algorithms in this file use a
#' decoder-local notation, standard in the sparse-reconstruction literature, in
#' which \eqn{\mathbf{A}} denotes the design matrix, \eqn{\mathbf{x}} denotes
#' the unknown sample-status vector, and \eqn{\mathbf{y}} denotes the response
#' vector. This local notation is restricted to the decoder layer and does not
#' replace the package-wide canonical symbols.
#'
#' Given a binary design matrix \eqn{\mathbf{A} \in \{0,1\}^{m \times n}} and
#' observed pooled outcomes \eqn{\mathbf{y} \in \{0,1\}^m}, the noiseless
#' decoder-local model is
#' \deqn{\mathbf{y} = \psi(\mathbf{A}\mathbf{x}),}
#' where \eqn{\mathbf{x} \in \{0,1\}^n} is the unknown sample-status vector and
#' \eqn{\psi(v) = \mathbf{1}_{v > 0}} is the element-wise OR map.
#'
#' COMP identifies samples as definitely negative if they appear in at least one
#' negative pool, and as candidates otherwise. It guarantees no false negatives
#' but may have false positives when more than \eqn{k_{\max}} samples are
#' positive.
#'
#' When COMP produces too many candidates, GPSR scores them by solving the
#' decoder-local sparse reconstruction problem
#' \deqn{\min_{\mathbf{x} \geq 0} \frac{1}{2}\|\mathbf{y} - \mathbf{A}\mathbf{x}\|_2^2 + \tau\|\mathbf{x}\|_1.}
#' Higher scores indicate higher likelihood of being truly positive.
#'
#' @name decoding_algorithms
#' @family decoding
#'
#' @references
#' Chan, C. L., et al. (2011). Non-adaptive probabilistic group testing with
#' noisy measurements. *Allerton Conference*.
#'
#' Figueiredo, M. A. T., Nowak, R. D., & Wright, S. J. (2007). Gradient
#' Projection for Sparse Reconstruction. *IEEE Journal of Selected Topics
#' in Signal Processing*, 1(4), 586-597.
#'
#' @seealso \code{\link{pp_matrix}} for creating pooling matrices
NULL


#' @title COMP Decoding Algorithm
#'
#' @description
#' Decodes pooled test results using the Combinatorial Orthogonal Matching
#' Pursuit (COMP) algorithm, which guarantees zero false negatives.
#'
#' @details
#' # Algorithm
#'
#' For each sample \eqn{j \in \{1, \ldots, n\}}:
#'
#' \deqn{\hat{x}_j = \begin{cases}
#'   0 \text{ (negative)} & \text{if } \exists i : M_{ij} = 1 \land y_i = 0 \\
#'   1 \text{ (candidate)} & \text{if } \forall i : M_{ij} = 1 \Rightarrow y_i = 1
#' \end{cases}}
#'
#' **In words**: A sample is marked negative if it appears in ANY negative pool;
#' it's a candidate (potentially positive) if ALL its pools are positive.
#'
#' # Theoretical Properties
#'
#' 1. **No false negatives**: COMP never misses a true positive
#' 2. **May have false positives**: Negative samples in only positive pools
#'    are incorrectly marked as candidates
#' 3. **Exact for PP designs**: When \eqn{k \leq k_{max}}, COMP is exact
#'    (no false positives) due to the intersection bound
#'
#' # Exactness Theorem for PP Matrices
#'
#' For a PP matrix with parameters \eqn{(q, d, n_l)}:
#'
#' If there are at most \eqn{k_{max} = \lfloor(n_l-1)/(d-1)\rfloor} true
#' positives, COMP identifies them exactly.
#'
#' **Proof**: Each negative sample shares at most \eqn{d-1} pools with each
#' positive (intersection bound). With \eqn{k} positives, a negative shares
#' at most \eqn{k(d-1)} pools. If \eqn{n_l > k(d-1)}, at least one pool
#' containing the negative is "clean" -> COMP marks it negative.
#'
#' # Time Complexity
#'
#' \eqn{O(m \times n)} where \eqn{m} = pools, \eqn{n} = samples.
#'
#' @param M Binary pooling matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}. It
#'   may be dense or sparse. Row \eqn{j} corresponds to pool \eqn{\mathcal{P}_j}
#'   and column \eqn{i} corresponds to individual \eqn{i}.
#' @param y Numeric vector of length \eqn{J} containing observed pool outcomes.
#'   Allowed inputs are either binary indicators \eqn{Z_j \in \{0,1\}} or
#'   quantitative Ct-like measurements encoded so that \eqn{0} means undetected
#'   and \eqn{0 < y_j < \texttt{ct\_threshold}} is treated as a positive pool.
#' @param ct_threshold Numeric scalar giving the positivity threshold applied when
#'   \code{y} contains Ct-like values. The default is \code{40}.
#'
#' @return A list with components:
#' \describe{
#'   \item{candidates}{Integer vector of potentially positive sample indices (1-indexed)}
#'   \item{negatives}{Integer vector of definitely negative sample indices}
#'   \item{n_positive_pools}{Integer; number of positive pools}
#' }
#'
#' @examples
#' # P-BEST example
#' M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#' M
#'
#' # Simulate 2 true positives
#' x_true <- rep(0L, 384)
#' x_true[c(72, 142)] <- 1L
#' x_true[c(72, 142)]
#'
#' # Compute pool outcomes
#' y <- as.integer(as.vector(M %*% x_true) > 0)
#' y
#' sum(y)
#'
#' # Decode with COMP
#' result <- pp_comp(M, y)
#' result
#' result$candidates
#' length(result$negatives)
#'
#' # Julia example_1.jl replication
#' M_small <- pp_matrix(q = 4, d = 3, nl = 5)
#' M_small
#' x_small <- rep(0L, 64)
#' x_small[c(2, 13)] <- 1L
#' x_small[c(2, 13)]
#' y_small <- as.integer(as.vector(M_small %*% x_small) > 0)
#' y_small
#' pp_comp(M_small, y_small)
#'
#' # With Ct values (simulated)
#' ct_values <- ifelse(y == 1, 28 + rnorm(length(y), 0, 3), 0)
#' result_ct <- pp_comp(M, ct_values, ct_threshold = 40)
#' result_ct
#'
#' # Verify no false negatives property
#' all(c(72, 142) %in% result$candidates)
#'
#' @references
#' Chan, C. L., Che, P. H., Jaggi, S., & Saligrama, V. (2011). Non-adaptive
#' probabilistic group testing with noisy measurements: Near-optimal bounds with
#' efficient algorithms. \emph{49th Annual Allerton Conference on Communication,
#' Control, and Computing}, 1832-1839.
#'
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' @seealso
#' \code{\link{pp_decode}}, \code{\link{pp_gpsr}}, \code{\link{pp_matrix}}
#'
#' @family decoding
#' @export
pp_comp <- function(M, y, ct_threshold = 40) {

  m <- nrow(M)
  n <- ncol(M)

  # Convert to binary
  if (any(y > 1, na.rm = TRUE)) {
    y_bin <- as.integer(y > 0 & y < ct_threshold)
  } else {
    y_bin <- as.integer(y > 0)
  }

  if (length(y_bin) != m) {
    stop("Length of y (", length(y_bin), ") must equal number of pools (", m, ")")
  }

  # COMP: candidate iff in >=1 pool AND not in any negative pool.
  # t(M) %*% (1 - y_bin) counts negative pools per sample; must be 0.
  x <- as.integer(
    (Matrix::colSums(M) > 0L) &
    (as.vector(Matrix::crossprod(M, 1L - y_bin)) == 0L)
  )

  list(
    candidates = which(x == 1L),
    negatives  = which(x == 0L),
    n_positive_pools = sum(y_bin)
  )
}


#' @title GPSR Compressed Sensing Algorithm
#'
#' @description
#' Solves the \eqn{\ell_1}-regularized least squares problem for scoring
#' candidate samples when COMP produces too many candidates.
#'
#' @details
#' GPSR (Gradient Projection for Sparse Reconstruction) solves the convex
#' optimization problem
#' \deqn{
#'   \min_{\mathbf{x} \geq 0}
#'   \frac{1}{2}\|\mathbf{y} - \mathbf{A}\mathbf{x}\|_2^2 + \tau\|\mathbf{x}\|_1.
#' }
#' The \eqn{\ell_1} penalty promotes sparse solutions while the non-negativity
#' constraint keeps the candidate scores interpretable as non-negative weights.
#'
#' @section Algorithm:
#' The implementation follows a GPSR-BB style projected-gradient scheme with a
#' short continuation warm start. Let
#' \eqn{F_\tau(\mathbf{z}) = \frac{1}{2}\|\mathbf{A}\mathbf{z} - \mathbf{y}\|_2^2 +
#' \tau \mathbf{1}^T\mathbf{z}} with gradient
#' \eqn{\nabla F_\tau(\mathbf{z}) = \mathbf{A}^T\mathbf{A}\mathbf{z} -
#' \mathbf{A}^T\mathbf{y} + \tau\mathbf{1}}. Before the main loop, the algorithm
#' applies \code{continuation_steps} projected-gradient updates at geometrically
#' decreasing regularization values from
#' \deqn{
#'   \tau_{\text{start}} =
#'   \texttt{first\_tau\_factor} \cdot \max|\mathbf{A}^T\mathbf{y}|
#' }
#' down to the target \eqn{\tau}. The main loop then alternates gradient
#' evaluation, projection onto the non-negative orthant, backtracking line
#' search, and a Barzilai-Borwein step-size update. Convergence is declared when
#' \deqn{
#'   \frac{|F_\tau(\mathbf{z}_{k+1}) - F_\tau(\mathbf{z}_k)|}{
#'     \max(1, |F_\tau(\mathbf{z}_k)|)
#'   } < \texttt{tol}.
#' }
#'
#' @section Role in Decoding:
#' When COMP produces more candidates than \eqn{k_{\max}}, GPSR assigns
#' continuous scores to those candidates. Larger scores indicate stronger
#' support for inclusion in a sparse positive set, but the scores themselves are
#' not final binary calls.
#'
#' @param A Numeric design matrix with dimensions \eqn{J \times N}. In the
#'   decoding pipeline it is typically the submatrix of positive pools by
#'   candidate individuals, but any numeric matrix is allowed.
#' @param y Numeric vector of length \eqn{J} containing the observed response
#'   associated with \code{A}. In the decoding pipeline this is usually a vector
#'   of normalized pool intensities or binary positive-pool indicators.
#' @param tau Numeric non-negative regularization scalar or \code{NULL}.
#'   When \code{NULL} (default), \eqn{\tau} is computed from the data as
#'   \deqn{\tau = 0.005 \cdot \max|\mathbf{A}^T \mathbf{y}|}
#'   using the package's default GPSR scaling. Supply a positive numeric value
#'   to override the data-adaptive default. Larger values promote sparser
#'   solutions.
#' @param max_iter Integer scalar giving the maximum number of projected-gradient
#'   iterations. The default is \code{500L}.
#' @param tol Numeric scalar giving the convergence tolerance for the relative
#'   change in objective value. The default is \code{1e-5}. Convergence is
#'   declared when
#'   \eqn{|F(\mathbf{z}_{\text{new}}) - F(\mathbf{z})| / \max(1, |F(\mathbf{z})|) < \texttt{tol}}.
#' @param continuation_steps Integer scalar giving the number of warm-start
#'   continuation steps before the main loop. The default is \code{5L}.
#' @param first_tau_factor Numeric scalar giving the scale factor that determines
#'   the starting regularization in the continuation schedule. The default is
#'   \code{0.8}. It sets the continuation starting value through
#'   \eqn{\tau_{\text{start}} = \texttt{first\_tau\_factor} \cdot \max|\mathbf{A}^T\mathbf{y}|}.
#'
#' @return A list with components:
#' \describe{
#'   \item{x}{Numeric vector; the estimated sparse solution}
#'   \item{iterations}{Integer; number of iterations performed}
#'   \item{converged}{Logical; TRUE if converged before max_iter}
#' }
#'
#' @examples
#' # Simple example
#' A <- matrix(c(1, 1, 0,
#'               0, 1, 1,
#'               1, 0, 1), nrow = 3, byrow = TRUE)
#' A
#' y <- c(1, 0.8, 0.5)
#' y
#'
#' result <- pp_gpsr(A, y)          # tau computed from data
#' result
#' result$x
#' result$converged
#'
#' # Typical use: score candidates from COMP
#' M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#' M
#' x_true <- rep(0L, 384)
#' x_true[c(72, 142, 200)] <- 1L  # 3 positives
#' x_true[c(72, 142, 200)]
#' y <- as.integer(as.vector(M %*% x_true) > 0)
#' y
#'
#' comp_result <- pp_comp(M, y)
#' comp_result
#' candidates <- comp_result$candidates
#' candidates
#'
#' # Score candidates using GPSR
#' pos_pools <- which(y == 1)
#' A_sub <- M[pos_pools, candidates]
#' A_sub
#' y_sub <- rep(1, length(pos_pools))
#' y_sub
#'
#' gpsr_result <- pp_gpsr(A_sub, y_sub)
#' gpsr_result
#' scores <- gpsr_result$x
#' scores
#' top_candidates <- candidates[order(scores, decreasing = TRUE)]
#' top_candidates
#'
#' @references
#' Figueiredo, M. A. T., Nowak, R. D., & Wright, S. J. (2007). Gradient
#' Projection for Sparse Reconstruction. \emph{IEEE Journal of Selected Topics
#' in Signal Processing}, 1(4), 586-597.
#'
#' Kim, H. Y., Hudgens, M. G., Dreyfuss, J. M., Westreich, D. J., & Pilcher,
#' C. D. (2007). Comparison of group testing algorithms for case identification
#' in the presence of test error. \emph{Biometrics}, 63(4), 1152-1163.
#'
#' @family decoding
#' @seealso
#' \code{\link{pp_comp}}, \code{\link{pp_decode}}, \code{\link{pp_matrix}}
#'
#' @export
pp_gpsr <- function(A, y, tau = NULL, max_iter = 500L, tol = 1e-5,
                    continuation_steps = 5L, first_tau_factor = 0.8) {

  y <- as.numeric(y)

  m <- nrow(A)
  n <- ncol(A)

  if (length(y) != m) {
    stop("Length of y (", length(y), ") must equal nrow(A) (", m, ")")
  }

  # Data-adaptive tau: 0.005 * max|A'y| - applyGPSR.m:7 / example_PBEST.m line 20
  Aty     <- as.vector(Matrix::crossprod(A, y))
  aty_max <- max(abs(Aty))
  if (is.null(tau)) {
    tau <- 0.005 * aty_max
    if (tau == 0) tau <- 1e-4
  }

  # Precompute AtA (used by both continuation and main loop)
  AtA <- crossprod(A)

  # Gradient and objective parameterised by tau_k (supports continuation)
  grad_k <- function(z, tau_k) as.vector(AtA %*% z) - Aty + tau_k
  obj_k  <- function(z, tau_k) {
    r <- as.vector(A %*% z) - y
    0.5 * sum(r * r) + tau_k * sum(z)
  }

  # Backtracking parameters (Armijo / step 4 of Algorithm F.2)
  beta  <- 0.5
  sigma <- 0.01

  # Initialization
  z     <- rep(0, n)
  alpha <- 1

  # Continuation warm-start: applyGPSR.m:8-9 (Continuation=1, ContinuationSteps=5,
  # FirstTauFactor=0.8 -> tau_start = first_tau_factor * max|A'y|, geometric to tau)
  if (continuation_steps > 0L && aty_max > 0 && first_tau_factor * aty_max > tau) {
    tau_start <- first_tau_factor * aty_max
    tau_seq   <- exp(seq(log(tau_start), log(tau), length.out = continuation_steps + 1L))
    tau_seq   <- tau_seq[seq_len(continuation_steps)]   # exclude endpoint (= target tau)
    for (tau_k in tau_seq) {
      g      <- grad_k(z, tau_k)
      Fz_k   <- obj_k(z, tau_k)
      delta  <- pmax(z - alpha * g, 0) - z
      gTd    <- sum(g * delta)
      lambda <- 1
      for (bt in seq_len(50L)) {
        if (obj_k(z + lambda * delta, tau_k) <= Fz_k + sigma * lambda * gTd) break
        lambda <- lambda * beta
      }
      z_new <- pmax(z + lambda * delta, 0)
      g_new <- grad_k(z_new, tau_k)
      dTd   <- sum(delta * delta)
      if (dTd > 1e-14) {
        alpha <- max(min(abs(sum(delta * g_new) / dTd), 1e10), 1e-10)
      }
      z <- z_new
    }
  }

  # Main GPSR-BB loop at target tau
  for (iter in seq_len(max_iter)) {

    g   <- grad_k(z, tau)
    Fz  <- obj_k(z, tau)

    # Projected direction: delta = P_+(z - alpha * g) - z
    delta <- pmax(z - alpha * g, 0) - z
    gTd   <- sum(g * delta)

    # Backtracking line search (step 4 of Algorithm F.2)
    lambda <- 1
    for (bt in seq_len(50L)) {
      if (obj_k(z + lambda * delta, tau) <= Fz + sigma * lambda * gTd) break
      lambda <- lambda * beta
    }

    z_new  <- pmax(z + lambda * delta, 0)
    g_new  <- grad_k(z_new, tau)

    # BB step update: alpha = (delta^T g_new) / ||delta||^2
    dTd <- sum(delta * delta)
    if (dTd > 1e-14) {
      alpha <- max(min(abs(sum(delta * g_new) / dTd), 1e10), 1e-10)
    }

    Fz_new <- obj_k(z_new, tau)
    z      <- z_new

    # Convergence: relative change in objective - applyGPSR.m:11 (stopCri=3, tolA=1e-5)
    if (abs(Fz_new - Fz) / max(1, abs(Fz)) < tol) {
      return(list(x = z, iterations = iter, converged = TRUE))
    }
  }

  list(x = z, iterations = max_iter, converged = FALSE)
}


#' @title Complete PP Decoding Pipeline
#'
#' @description
#' Full detection pipeline combining COMP filtering, GPSR scoring, exhaustive
#' verification, and an optional three-way probabilistic classification
#' (positive / suspected / negative) following Zismanov et al. (2024).
#'
#' @details
#' The current package workflow applies COMP filtering, optional GPSR scoring,
#' exact subset verification, and an optional three-way classification layer.
#'
#' @section Algorithm Steps:
#' \enumerate{
#'   \item \strong{COMP filtering}: identify sure negatives and candidates.
#'   \item \strong{Early exit}: if the candidate count is
#'   \eqn{\leq k_{max}}, accept all candidates as decoded positives.
#'   \item \strong{GPSR scoring}: if there are too many candidates, rank them by
#'   GPSR scores.
#'   \item \strong{Exhaustive verification}: search the top candidate sets for
#'   the best subset under
#'   \eqn{\|\psi(\mathbf{M}\hat{\mathbf{x}}) - \mathbf{y}\|_1}.
#'   \item \strong{Three-way classification}: optionally split decoded
#'   candidates into confirmed positives and suspected samples that require
#'   retesting.
#' }
#'
#' @section Handling Ct Values:
#' Pool results in \code{y} can be supplied either as binary outcomes or as
#' Ct-like values interpreted with \code{ct_threshold}.
#' \itemize{
#'   \item \strong{Binary}: \eqn{Z_j \in \{0,1\}}.
#'   \item \strong{Ct values}: 0 indicates an undetected pool, and
#'   \eqn{0 < \mathrm{Ct} < 40} indicates a detected pool.
#' }
#' When quantitative Ct-like values are available, the package GPSR scoring
#' layer uses them as a proxy for signal strength.
#'
#' @section Noise Model:
#' When per-sample Ct values are supplied through \code{ct}, the current
#' package implementation models the false-negative probability as
#' \deqn{P(\text{FN} \mid \text{Ct}) =
#'   \frac{1}{1 + \exp\bigl(-(\text{Ct} - 29) \times 0.5\bigr)}.}
#' This places \eqn{P(\text{FN}) = 0.5} at \eqn{\text{Ct} = 29} and
#' \eqn{P(\text{FN}) \approx 0.95} near \eqn{\text{Ct} \approx 35}, which is
#' consistent with the paper's strong-positive threshold
#' \eqn{\text{Ct} \leq 36}.
#'
#' @section Three-Way Classification:
#' Each decoded candidate is assigned to one of three categories:
#' \describe{
#'   \item{\strong{Positive}}{All identifying pools are positive; Ct is
#'   \eqn{\leq 36} when available; \eqn{P(\text{FN} \mid \text{Ct}) \leq 0.95};
#'   and the sample has at least one unique identifying pool.}
#'   \item{\strong{Suspected}}{Pool results are compatible with positivity, but
#'   Ct is above 36, no unique identifying pool is present, or
#'   \eqn{P(\text{FN}) > 0.95}.}
#'   \item{\strong{Negative}}{The sample is eliminated by COMP because at least
#'   one identifying pool is negative.}
#' }
#' A unique identifying pool is a positive pool that contains the sample but no
#' other decoded positive, matching criterion (b) in the paper. When
#' \code{ct = NULL}, the Ct and false-negative criteria are skipped and the
#' uniqueness rule still applies.
#'
#' @param M Binary pooling matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}},
#'   typically created by \code{\link{pp_matrix}}. Its rows index pools
#'   \eqn{j \in \{1,\ldots,J\}} and its columns index individuals
#'   \eqn{i \in \{1,\ldots,N\}}.
#' @param y Numeric vector of length \eqn{J} containing observed pool outcomes.
#'   Allowed inputs are either binary observations \eqn{Z_j \in \{0,1\}} or
#'   Ct-like values interpreted using \code{ct_threshold}.
#' @param k_max Integer scalar or \code{NULL} giving the maximum number of
#'   positives to consider during exact verification.
#'   When \code{NULL} (default), the value is read from
#'   \code{attr(M, "pp")$k_max} if present.  A hard numeric value
#'   supplied by the caller always takes precedence.  If neither source
#'   provides a value, a fallback of 5 is used with a warning.
#' @param ct_threshold Numeric scalar giving the positivity threshold used when
#'   \code{y} contains Ct-like values. The default is \code{40}.
#' @param verbose Logical scalar; if \code{TRUE} (default), print progress
#'   messages describing the decoding path taken.
#' @param ct Numeric vector of length \eqn{N} containing per-individual Ct
#'   values, or \code{NULL}. Named or unnamed; NAs are treated as unknown (Ct criterion
#'   skipped for that sample).  When non-\code{NULL}, enables the noise model
#'   and three-way classification.  When \code{NULL} (default), the function
#'   behaves identically to previous versions except that the uniqueness check
#'   still populates \code{$suspected}.
#'
#' @return A list with components:
#' \describe{
#'   \item{positives}{Integer vector of confirmed positive sample indices
#'     (1-indexed).  These meet all three-way classification criteria.}
#'   \item{n_positives}{Integer; number of confirmed positives.}
#'   \item{suspected}{Integer vector of suspected sample indices (require
#'     retesting).  Empty when all decoded positives pass classification.}
#'   \item{n_suspected}{Integer; number of suspected samples.}
#'   \item{candidates}{Integer vector of all COMP candidates.}
#'   \item{n_candidates}{Integer; number of COMP candidates.}
#'   \item{error}{Integer; number of pool mismatches (0 = perfect match).}
#'   \item{method}{Character; description of method used.}
#' }
#'
#' @examples
#' # Generate P-BEST matrix
#' M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#' M
#'
#' # Simulate 2 positives (within k_max)
#' true_pos <- c(72, 142)
#' x <- rep(0L, 384)
#' x[true_pos] <- 1L
#' y <- as.integer(as.vector(M %*% x) > 0)
#' y
#'
#' # Decode (binary, no Ct) - backward-compatible call
#' result <- pp_decode(M, y)
#' result
#' result$positives   # 72, 142
#' result$suspected   # integer(0) when both have unique pools
#' result$error       # 0
#'
#' # Validate result
#' pp_validate(result$positives, true_pos)
#'
#' # With per-sample Ct values: enables noise model + three-way classification
#' ct_vals <- rep(NA_real_, 384)
#' ct_vals[true_pos] <- c(24, 31)   # sample 72 strong, sample 142 moderate
#' ct_vals[true_pos]
#' result_ct <- pp_decode(M, y, ct = ct_vals, verbose = FALSE)
#' result_ct
#' result_ct$positives   # 72 (Ct=24, strong positive)
#' result_ct$suspected   # 142 (Ct=31 > 29: P(FN) > 0.5, borderline)
#'
#' # With simulated pool Ct values (passed as y)
#' ct_pool <- ifelse(y == 1, 30 + rnorm(48, 0, 2), 0)
#' result_pool_ct <- pp_decode(M, ct_pool, verbose = FALSE)
#' result_pool_ct
#'
#' # Decode with more positives (may need GPSR)
#' true_pos_3 <- c(72, 142, 200)
#' x3 <- rep(0L, 384)
#' x3[true_pos_3] <- 1L
#' y3 <- as.integer(as.vector(M %*% x3) > 0)
#' result3 <- pp_decode(M, y3, k_max = 5)
#' result3
#'
#' # Julia example
#' M_julia <- pp_matrix(q = 4, d = 3, nl = 5)
#' M_julia
#' x_julia <- rep(0L, 64)
#' x_julia[c(2, 13)] <- 1L
#' y_julia <- as.integer(as.vector(M_julia %*% x_julia) > 0)
#' pp_decode(M_julia, y_julia, verbose = FALSE)$positives  # 2, 13
#'
#' @references
#' Zismanov, S., Yelin, I., Klochendler, A., et al. (2024). High capacity
#' clinical SARS-CoV-2 molecular testing using combinatorial pooling.
#' \emph{Communications Medicine}, 4, Article 121.
#' \doi{10.1038/s43856-024-00531-w}. Decoding algorithm and sample
#' classification sections.
#'
#' Figueiredo, M. A. T., Nowak, R. D., & Wright, S. J. (2007). Gradient
#' Projection for Sparse Reconstruction. \emph{IEEE Journal of Selected Topics
#' in Signal Processing}, 1(4), 586-597.
#'
#' @family decoding
#' @seealso
#' \code{\link{pp_comp}}, \code{\link{pp_gpsr}}, \code{\link{pp_validate}},
#' \code{\link{pp_matrix}}
#' @export
pp_decode <- function(M, y, k_max = NULL, ct_threshold = 40, verbose = TRUE,
                      ct = NULL) {

  # Resolve k_max: caller > matrix attr > fallback
  if (is.null(k_max)) {
    pp_attr <- attr(M, "pp")
    if (!is.null(pp_attr) && !is.null(pp_attr$k_max)) {
      k_max <- pp_attr$k_max
    } else {
      warning("k_max not supplied and attr(M, \"pp\")$k_max not found; ",
              "using fallback k_max = 5. ",
              "Pass k_max explicitly or use a matrix created by pp_matrix().")
      k_max <- 5L
    }
  }
  k_max <- as.integer(k_max)

  m <- nrow(M)
  n <- ncol(M)

  if (!is.null(ct) && length(ct) != n) {
    stop("Length of ct (", length(ct), ") must equal ncol(M) (", n, ")")
  }

  # Process y
  if (any(y > 1, na.rm = TRUE)) {
    y_bin   <- as.integer(y > 0 & y < ct_threshold)
    y_quant <- y
  } else {
    y_bin   <- as.integer(y > 0)
    y_quant <- y
  }

  if (verbose) {
    message("\n=== PP Decode ===")
    message("Matrix: ", m, " pools x ", n, " samples")
    message("Positive pools: ", sum(y_bin), "/", m)
  }

  # Step 1: COMP
  comp       <- pp_comp(M, y_bin)
  candidates <- comp$candidates

  if (verbose) {
    message("COMP candidates: ", length(candidates))
  }

  # No positive pools at all: nothing to decode
  if (sum(y_bin) == 0L) {
    if (verbose) message("Result: No positives detected (all pools negative)")
    return(list(
      positives = integer(0L), n_positives = 0L,
      suspected = integer(0L), n_suspected = 0L,
      candidates = candidates, n_candidates = 0L,
      error = 0L, method = "COMP (all pools negative)"
    ))
  }

  # COMP exact: few candidates that perfectly explain the pool outcomes
  if (length(candidates) > 0L && length(candidates) <= k_max) {
    err <- verify_pools_internal(M, y_bin, candidates)
    if (err == 0L) {
      cls <- classify_threeway_internal(M, y_bin, candidates, ct)
      if (verbose) {
        message("Result: ", paste(cls$confirmed, collapse = ", "))
        if (length(cls$suspected) > 0L)
          message("Suspected: ", paste(cls$suspected, collapse = ", "))
        message("Method: COMP (exact)")
      }
      return(list(
        positives  = cls$confirmed,  n_positives  = length(cls$confirmed),
        suspected  = cls$suspected,  n_suspected  = length(cls$suspected),
        candidates = candidates,     n_candidates = length(candidates),
        error = 0L, method = "COMP (exact)"
      ))
    }
  }

  # If COMP returns a modest candidate set, try exhaustive verification over
  # the COMP support before invoking GPSR. This preserves exact recovery when a
  # small zero-error subset exists even if COMP itself returns extra candidates.
  if (length(candidates) > k_max && length(candidates) <= 20L) {
    best_comp <- search_best_internal(M, y_bin, candidates, k_max)
    if (best_comp$error == 0L) {
      cls <- classify_threeway_internal(M, y_bin, best_comp$positives, ct)
      if (verbose) {
        message("Result: ", paste(cls$confirmed, collapse = ", "))
        if (length(cls$suspected) > 0L)
          message("Suspected: ", paste(cls$suspected, collapse = ", "))
        message("Method: COMP + Verification")
      }
      return(list(
        positives  = cls$confirmed,         n_positives  = length(cls$confirmed),
        suspected  = cls$suspected,         n_suspected  = length(cls$suspected),
        candidates = candidates,            n_candidates = length(candidates),
        error = 0L, method = "COMP + Verification"
      ))
    }
  }

  # Step 2: GPSR on the full matrix and full measurement vector.
  # Matches example_PBEST.m lines 19-21:
  #   dt  = max(abs(poolingMatrix' * qMeasurement))
  #   tau = 0.005 * dt
  #   u   = opm(qMeasurement, poolingMatrix, tau, maxNum)
  # Running on the full M (not a COMP-filtered submatrix) ensures that positive
  # samples whose pools all overlap with other positives are still recoverable.
  if (verbose) message("Running GPSR refinement...")

  tau_val <- 0.005 * max(abs(as.vector(Matrix::crossprod(M, y_bin))))
  y_gpsr  <- if (max(abs(y_quant)) > 1) y_quant / max(abs(y_quant)) else as.numeric(y_bin)

  gpsr <- tryCatch(
    pp_gpsr(M, y_gpsr, tau = tau_val),
    error = function(e) list(x = rep(0, n))
  )

  # Top non-zero entries by magnitude - mirrors opm() returning up to maxNum = 20
  scores   <- gpsr$x
  nonzero  <- which(scores > 1e-6)
  n_top    <- min(20L, max(1L, length(nonzero)))
  top_cand <- order(scores, decreasing = TRUE)[seq_len(n_top)]

  # Step 3: Exhaustive verification
  best <- search_best_internal(M, y_bin, top_cand, k_max)

  cls <- classify_threeway_internal(M, y_bin, best$positives, ct)

  if (verbose) {
    message("Result: ", paste(cls$confirmed, collapse = ", "))
    if (length(cls$suspected) > 0L)
      message("Suspected: ", paste(cls$suspected, collapse = ", "))
    message("Method: COMP + GPSR + Verification")
    message("Pool mismatches: ", best$error)
  }

  list(
    positives  = cls$confirmed,       n_positives  = length(cls$confirmed),
    suspected  = cls$suspected,       n_suspected  = length(cls$suspected),
    candidates = candidates,          n_candidates = length(candidates),
    error = best$error, method = "COMP + GPSR + Verification"
  )
}


#' @title Validate Decoding Results
#'
#' @description
#' Compares decoded results against known true positives and computes
#' performance metrics.
#'
#' @details
#' # Metrics Computed
#'
#' - **True Positives (TP)**: Correctly identified positives
#' - **False Positives (FP)**: Incorrectly identified as positive
#' - **False Negatives (FN)**: Missed positives
#' - **Sensitivity**: TP / (TP + FN) = fraction of positives found
#' - **Precision**: TP / (TP + FP) = fraction of identified that are truly positive
#' - **Exact Match**: TRUE if perfect identification (TP = all true, FP = 0)
#'
#' # COMP Properties
#'
#' For COMP decoding with PP matrices:
#' - **FN is always 0** (no false negatives guaranteed)
#' - **Sensitivity is always 1** (100%)
#' - **Exact match** when \eqn{k \leq k_{max}}
#'
#' @param decoded Integer vector containing the decoded positive individual
#'   indices, that is, a subset of \eqn{\{1,\ldots,N\}} returned by a decoder.
#' @param true_pos Integer vector containing the true positive individual
#'   indices, that is, the ground-truth active set \eqn{\mathcal{S}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{TP}{Integer; true positives}
#'   \item{FP}{Integer; false positives}
#'   \item{FN}{Integer; false negatives}
#'   \item{sensitivity}{Numeric; TP / (TP + FN), in \eqn{[0,1]}}
#'   \item{precision}{Numeric; TP / (TP + FP), in \eqn{[0,1]}}
#'   \item{exact_match}{Logical; TRUE if perfect identification}
#' }
#'
#' @examples
#' # Perfect match
#' pp_validate(c(72, 142), c(72, 142))
#'
#' # One false positive
#' pp_validate(c(72, 142, 50), c(72, 142))
#'
#' # One false negative (unusual for COMP)
#' pp_validate(c(72), c(72, 142))
#'
#' # Complete example
#' M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#' true_pos <- c(72, 142)
#' x <- rep(0L, 384)
#' x[true_pos] <- 1L
#' y <- as.integer(as.vector(M %*% x) > 0)
#'
#' result <- pp_decode(M, y, verbose = FALSE)
#' metrics <- pp_validate(result$positives, true_pos)
#' metrics$exact_match  # TRUE
#' metrics$sensitivity  # 1.0
#'
#' @references
#' Altman, D. G., & Bland, J. M. (1994). Diagnostic tests 1: Sensitivity and
#' specificity. \emph{BMJ}, 308(6943), 1552.
#'
#' Zhou, X.-H., Obuchowski, N. A., & McClish, D. K. (2011). \emph{Statistical
#' Methods in Diagnostic Medicine} (2nd ed.). Wiley.
#'
#' @seealso
#' \code{\link{pp_decode}}, \code{\link{pp_comp}}
#'
#' @family decoding
#' @export
pp_validate <- function(decoded, true_pos) {
  TP <- length(intersect(decoded, true_pos))
  FP <- length(setdiff(decoded, true_pos))
  FN <- length(setdiff(true_pos, decoded))

  sens <- if (length(true_pos) > 0) TP / length(true_pos) else 1
  prec <- if (length(decoded) > 0) TP / length(decoded) else 1

  list(
    TP = TP, FP = FP, FN = FN,
    sensitivity = sens,
    precision = prec,
    exact_match = (TP == length(true_pos) && FP == 0)
  )
}


# =============================================================================
# Internal helper functions
# =============================================================================

# Internal: three-way classification - Zismanov et al. 2024, Sample Classification.
# Splits decoded positives into confirmed and suspected.
# Criteria for confirmed positive (all must hold):
#   (a) ct <= ct_strong (if ct supplied and not NA)
#   (b) P(FN|ct) <= fn_threshold, where P(FN|Ct) = 1/(1+exp(-(Ct-29)*0.5))
#   (c) sample has >=1 unique identifying pool: a positive pool shared with no
#       other decoded positive (paper criterion b, overcrowding check)
classify_threeway_internal <- function(M, y_bin, positives,
                                       ct          = NULL,
                                       ct_strong   = 36,
                                       fn_threshold = 0.95) {
  if (length(positives) == 0L) {
    return(list(confirmed = integer(0L), suspected = integer(0L)))
  }

  pos_pools <- which(y_bin == 1L)

  # Per-pool count of how many decoded positives each positive pool contains
  # (used for the uniqueness check)
  if (length(pos_pools) > 0L) {
    pool_pos_count <- as.integer(
      Matrix::rowSums(M[pos_pools, positives, drop = FALSE])
    )
  } else {
    pool_pos_count <- integer(0L)
  }

  has_unique <- logical(length(positives))
  for (jj in seq_along(positives)) {
    j  <- positives[jj]
    # Indices into pos_pools where sample j appears
    pj <- which(as.logical(M[pos_pools, j]))
    has_unique[jj] <- length(pj) > 0L && any(pool_pos_count[pj] == 1L)
  }

  # Ct and FN-probability checks (per sample)
  ct_ok <- rep(TRUE, length(positives))
  fn_ok <- rep(TRUE, length(positives))

  if (!is.null(ct)) {
    ct_j  <- ct[positives]
    valid <- !is.na(ct_j)
    if (any(valid)) {
      ct_ok[valid] <- ct_j[valid] <= ct_strong
      p_fn         <- 1 / (1 + exp(-(ct_j[valid] - 29) * 0.5))
      fn_ok[valid] <- p_fn <= fn_threshold
    }
  }

  is_confirmed <- ct_ok & fn_ok & has_unique

  list(
    confirmed = positives[is_confirmed],
    suspected = positives[!is_confirmed]
  )
}


# Internal: verify pool mismatches
verify_pools_internal <- function(M, y_bin, positives) {
  x <- rep(0L, ncol(M))
  x[positives] <- 1L
  pred <- as.integer(as.vector(M %*% x) > 0)
  sum(abs(pred - y_bin))
}

# Internal: exhaustive search for best solution
# Evaluates all subsets of size 1..k_max over the top 20 candidates.
# Mirrors selectByError.m:8-18 which tries all 2^maxNum - 1 non-empty subsets
# (maxNum = 20).  No subset size is skipped.
search_best_internal <- function(M, y_bin, candidates, k_max) {
  n_cand <- min(20L, length(candidates))
  search <- candidates[seq_len(n_cand)]

  best_err <- Inf
  best_pos <- integer(0)

  for (k in seq_len(min(k_max, n_cand))) {
    combs <- combn(n_cand, k)
    for (j in seq_len(ncol(combs))) {
      test <- search[combs[, j]]
      err  <- verify_pools_internal(M, y_bin, test)

      if (err < best_err) {
        best_err <- err
        best_pos <- test
      }
      if (err == 0L) break
    }
    if (best_err == 0L) break
  }

  list(positives = best_pos, error = best_err)
}
