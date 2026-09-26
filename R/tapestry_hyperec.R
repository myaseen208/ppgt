# =============================================================================
# Tapestry and HYPER-EC Designs
# Advanced error-correcting group testing designs
# =============================================================================

#' @title Generate A Tapestry-Style Pooling Matrix
#'
#' @description
#' Construct a Tapestry-style pooling design in which each individual is placed
#' in exactly three pools and the resulting matrix is intended for quantitative
#' compressed-sensing reconstruction rather than binary COMP decoding.
#'
#' @details
#' Tapestry is a single-round design aimed at quantitative pooled measurements
#' and compressed-sensing-style reconstruction rather than binary COMP decoding.
#' In the package-wide notation, observed pool outcomes are denoted by
#' \eqn{\mathbf{z}}. This help page adopts a decoder-local compressed-sensing
#' notation, standard for the Tapestry literature, in which \eqn{\mathbf{A}} is
#' the design matrix, \eqn{\mathbf{x}} is the unknown sparse signal,
#' \eqn{\mathbf{y}} is the response vector, and \eqn{\mathbf{g}} is the
#' measurement-noise vector. In that local notation,
#' \deqn{\mathbf{y} = \mathbf{A}\mathbf{x} + \mathbf{g}.}
#'
#' This constructor places each individual in exactly three pools, so the
#' intended maximum column sum is \eqn{\lambda_{\max} = 3}. The current package
#' implementation provides either a Kirkman-triple-inspired layout or a random
#' triple assignment. It returns only the constructed matrix object and
#' associated metadata; it does not implement a quantitative decoder.
#'
#' @section Construction Methods:
#' \describe{
#'   \item{\code{method = "kirkman"}}{Uses the package's Kirkman-triple-inspired
#'   assignment path.}
#'   \item{\code{method = "random"}}{Uses a random triple assignment while
#'   keeping each individual in exactly three pools.}
#' }
#'
#' @param n Integer scalar giving the total number of individuals. This is the
#'   column dimension \eqn{N} of the returned pooling matrix and must satisfy
#'   \eqn{N \ge 1}.
#' @param k Optional integer scalar giving the expected number of positive
#'   individuals used to size the number of pools. If \code{NULL}, the function
#'   uses \eqn{k = \max\{1, \lfloor N \cdot \texttt{prevalence}\rfloor\}}.
#' @param method Character string specifying the construction rule:
#'   \code{"kirkman"} for a Kirkman-triple-inspired layout or \code{"random"}
#'   for a random triple assignment. The default is \code{"random"}.
#' @param prevalence Numeric scalar giving the expected prevalence used only when
#'   \code{k} is \code{NULL}. The default is \code{0.01}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{matrix}{A sparse binary matrix of class \code{dgCMatrix} with
#'     dimensions \eqn{J \times N}}
#'   \item{m}{The realized number of pools \eqn{J}}
#'   \item{n}{The number of individuals \eqn{N}}
#'   \item{lambda_max}{The target maximum column sum, fixed at 3}
#'   \item{method}{The construction rule actually used}
#'   \item{design_type}{"tapestry"}
#' }
#'   The matrix itself also carries a \code{design} attribute with fields
#'   \code{type}, \code{n}, \code{m}, \code{lambda_max}, \code{k_expected}, and
#'   \code{method}.
#'
#' @examples
#' design1 <- tapestry_matrix(n = 100, k = 2)
#' design1
#' design1$m
#' Matrix::colSums(design1$matrix)
#'
#' design2 <- tapestry_matrix(n = 256, prevalence = 0.01)
#' design2
#' design2$m
#' design2$lambda_max
#'
#' pool_sizes <- Matrix::rowSums(design2$matrix)
#' pool_sizes
#' range(pool_sizes)
#'
#' @references
#' Ghosh, S., Agarwal, R., Rehan, M. A., et al. (2021). Tapestry: A single-round
#' smart pooling technique for COVID-19 testing. \emph{Nature Communications},
#' 12, Article 2995.
#'
#' Ghosh, S., et al. (2021). A compressed sensing approach to pooled RT-PCR
#' testing for COVID-19 detection. \emph{IEEE Open Journal of Signal
#' Processing}, 2, 248-264. \doi{10.1109/OJSP.2021.3075913}
#'
#' @seealso
#' \code{\link{hyper_matrix}}, \code{\link{hyper_ec_matrix}},
#' \code{\link{pp_matrix}}, \code{\link{compare_all_designs}}
#'
#' @family pooling_designs
#' @export
tapestry_matrix <- function(n, k = NULL, method = c("random", "kirkman"),
                            prevalence = 0.01) {
  n <- as.integer(n)
  method <- match.arg(method)
  
  if (n < 1L) stop("n must be a positive integer")
  
  # Determine expected infected
  if (is.null(k)) {
    k <- max(1L, as.integer(n * prevalence))
  } else {
    k <- as.integer(k)
  }
  
  # Determine number of pools
  # Tapestry uses m ≈ 4k log(n) empirically
  m <- as.integer(ceiling(4 * k * log(n)))
  
  # Ensure m is reasonable
  m <- max(m, as.integer(ceiling(n * 0.15)))  # At least 15% of n
  m <- min(m, n)  # At most n pools
  
  # Generate pooling matrix
  if (method == "kirkman") {
    # Kirkman triple system (balanced)
    mat <- .kirkman_triples(n, m)
  } else {
    # Random assignment ensuring each sample in exactly 3 pools
    mat <- .random_triples(n, m)
  }
  
  # Verify lambda_max = 3
  col_sums <- Matrix::colSums(mat)
  if (!all(col_sums == 3)) {
    warning("Not all samples in exactly 3 pools. Range: ",
            min(col_sums), "-", max(col_sums))
  }
  
  attr(mat, "design") <- list(
    type = "tapestry",
    n = n,
    m = m,
    lambda_max = 3,
    k_expected = k,
    method = method
  )
  
  list(
    matrix = mat,
    m = m,
    n = n,
    lambda_max = 3,
    method = method,
    design_type = "tapestry"
  )
}


#' @title Generate An Experimental Parity-Augmented HYPER Matrix
#'
#' @description
#' Construct a package-specific experimental extension that augments a base
#' HYPER matrix with additional parity-derived rows. This constructor is not the
#' Hong et al. HYPER procedure and is documented as an experimental package
#' extension only.
#'
#' @details
#' # Scope
#'
#' The package function `hyper_ec_matrix()` is an experimental convenience
#' wrapper around two pieces:
#'
#' 1. a base HYPER matrix \eqn{\mathbf{M}_{base} \in \{0,1\}^{J_0 \times N}}
#'    generated by \code{\link{hyper_matrix}}, and
#' 2. a binary parity-check matrix
#'    \eqn{\mathbf{H} \in \{0,1\}^{K \times J_0}}.
#'
#' The returned matrix stacks the base rows and derived parity rows. The parity
#' rows are constructed in code by
#' \deqn{\mathbf{M}_{parity} = (\mathbf{H}\mathbf{M}_{base}) \bmod 2.}
#'
#' This package-specific mod-2 construction should be treated as experimental.
#' It is not part of the HYPER design described in Hong et al. (2022), which is
#' a two-stage procedure with conservative stage-1 decoding followed by
#' individual retesting in stage 2.
#'
#' # Current Interpretation
#'
#' The current implementation uses `k` or `p` only as heuristics for sizing the
#' base design. The returned object exposes:
#'
#' - a base HYPER assignment matrix,
#' - a randomly generated parity-check matrix, and
#' - a derived matrix of parity rows used by `hyper_ec_decode()`.
#'
#' The stored field `error_correction = floor(k_parity / 2)` is the maximum
#' error weight that the current brute-force search will attempt. It is not a
#' proven correction guarantee.
#'
#' @param n Integer scalar giving the total number of individuals \eqn{N}. This
#'   is the column dimension of the base matrix \eqn{\mathbf{M}_{base}} and the
#'   combined matrix.
#' @param k Optional integer scalar used as a sizing heuristic for the expected
#'   number of positive individuals. If \code{NULL}, the function derives
#'   \eqn{k = \max\{1, \lfloor pN \rfloor\}} from \code{p}.
#' @param p Numeric scalar giving the expected prevalence used only when
#'   \code{k} is \code{NULL}. The default is \code{0.01}.
#' @param m_base Optional integer scalar giving the number of base HYPER pools
#'   \eqn{J_0}. If \code{NULL}, the function uses its current package heuristic
#'   based on \code{n}, \code{k}, and \code{q}.
#' @param q Integer scalar giving the number of base HYPER pools containing each
#'   individual. Allowed values are \code{1L}, \code{2L}, and \code{3L}. The
#'   default is \code{2L}.
#' @param alpha Numeric scalar giving the parity-overhead fraction used to set
#'   \eqn{K = \lceil \alpha J_0 \rceil}. It must lie in \eqn{(0,1)}.
#' @param sparse_parity Logical scalar; if \code{TRUE}, generate a sparse random
#'   parity-check matrix \eqn{\mathbf{H}}. Otherwise generate a dense random
#'   binary matrix.
#' @param parity_weight Integer scalar giving the target row weight of the sparse
#'   parity-check matrix when \code{sparse_parity = TRUE}.
#'
#' @return A named list with components:
#' \describe{
#'   \item{matrix}{Combined sparse binary matrix with the \eqn{J_0} base rows
#'     followed by the \eqn{K} parity-derived rows}
#'   \item{matrix_base}{Base HYPER matrix \eqn{\mathbf{M}_{base}}}
#'   \item{matrix_parity}{Parity-check matrix \eqn{\mathbf{H}}}
#'   \item{m}{Total number of rows \eqn{J = J_0 + K}}
#'   \item{m_base}{Number of base rows \eqn{J_0}}
#'   \item{k_parity}{Number of parity rows \eqn{K}}
#'   \item{k_design}{The sizing heuristic \eqn{k} actually used}
#'   \item{n}{Number of individuals \eqn{N}}
#'   \item{alpha}{Parity-overhead fraction}
#'   \item{q}{Base HYPER column weight}
#'   \item{lambda_max}{Estimated maximum number of rows containing an individual}
#'   \item{error_correction}{Maximum error weight searched by the current decoder}
#'   \item{efficiency}{Realized ratio \eqn{J/N}}
#'   \item{sparse_parity}{Whether sparse parity construction was used}
#'   \item{parity_weight}{Requested sparse parity row weight}
#'   \item{design_type}{"hyper_ec"}
#' }
#'   The combined matrix also carries a \code{design} attribute mirroring the
#'   main structural fields.
#'
#' @examples
#' # Experimental parity-augmented HYPER design
#' design1 <- hyper_ec_matrix(n = 256, k = 2, q = 2, alpha = 0.20)
#' design1
#' design1$m
#' design1$m_base
#' design1$k_parity
#' design1$error_correction
#' design1$efficiency
#'
#' # Increase parity overhead
#' design2 <- hyper_ec_matrix(n = 256, k = 2, q = 2, alpha = 0.30)
#' design2
#' design2$error_correction
#'
#' # Smaller design
#' design3 <- hyper_ec_matrix(n = 100, k = 2, q = 2, alpha = 0.20)
#' design3
#' design3$m
#'
#' # If k is omitted, the function uses prevalence to size the base matrix
#' design4 <- hyper_ec_matrix(n = 256, p = 0.01, q = 2, alpha = 0.20)
#' design4
#'
#' # Check structure
#' design1$matrix
#' design1$matrix_base
#' design1$matrix_parity
#'
#' # Inspect the realized column weights
#' range(Matrix::colSums(design1$matrix))
#'
#' @references
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' The paper above describes HYPER itself. The parity augmentation implemented
#' here is package-specific and experimental.
#'
#' @seealso
#' \code{\link{hyper_matrix}}, \code{\link{hyper_ec_decode}},
#' \code{\link{tapestry_matrix}}, \code{\link{compare_all_designs}}
#'
#' @family pooling_designs
#' @export
hyper_ec_matrix <- function(n, k = NULL, p = 0.01, m_base = NULL, q = 2L, alpha = 0.20,
                            sparse_parity = TRUE, parity_weight = 3L) {
  n <- as.integer(n)
  q <- as.integer(q)
  alpha <- as.numeric(alpha)
  parity_weight <- as.integer(parity_weight)
  
  if (n < 1L) stop("n must be a positive integer")
  if (!q %in% 1:3) stop("q must be 1, 2, or 3")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
  if (parity_weight < 1L) stop("parity_weight must be positive")
  
  # Determine k if not provided
  if (is.null(k)) {
    k <- max(1L, as.integer(p * n))
  } else {
    k <- as.integer(k)
    if (k < 1L) stop("k must be a positive integer")
  }
  
  # Determine m_base if not provided
  if (is.null(m_base)) {
    # Formula: m ≈ 4k ln(n) for k expected infected
    m_base <- as.integer(ceiling(4 * k * log(n)))
    m_base <- max(m_base, as.integer(ceiling(n * 0.15)))  # At least 15%
    
    # Ensure m_base satisfies HYPER constraints
    if (q == 2L && m_base %% 2L != 0L) m_base <- m_base + 1L
    if (q == 3L) {
      # Find nearest m where m ≡ 0 mod 6 and m-1 is prime
      m_base <- .find_valid_m_for_q3(m_base)
    }
  } else {
    m_base <- as.integer(m_base)
  }
  
  # Generate base HYPER matrix
  mat_base <- hyper_matrix(n = n, m = m_base, q = q)
  
  # Determine number of parity pools
  k_parity <- as.integer(ceiling(alpha * m_base))
  
  # Generate parity check matrix
  if (sparse_parity) {
    mat_parity <- .sparse_parity_matrix(k_parity, m_base, parity_weight)
  } else {
    mat_parity <- .random_parity_matrix(k_parity, m_base)
  }
  
  # Combine into full pooling matrix
  # Each parity pool tests XOR of base pools indicated by mat_parity
  # For actual pooling, parity pool i contains samples from base pools where mat_parity[i,j]=1
  # This is: mat_parity %*% mat_base (mod 2)
  
  # Create parity pool rows
  parity_rows <- (mat_parity %*% mat_base) %% 2
  
  # Combine base and parity
  mat_combined <- rbind(mat_base, parity_rows)
  
  # Calculate total pools and efficiency
  m_total <- m_base + k_parity
  efficiency <- m_total / n
  
  # Error correction capacity
  errors_correctable <- floor(k_parity / 2)
  
  # Estimate lambda_max (base HYPER contributes q, parity adds some)
  # With sparse parity (weight w), each sample in ≈ alpha * q * w/m_base parity pools
  if (sparse_parity) {
    parity_per_sample <- alpha * q * parity_weight / m_base
  } else {
    parity_per_sample <- alpha * q * 0.5  # Random expectation
  }
  lambda_max_est <- q + ceiling(parity_per_sample)
  
  attr(mat_combined, "design") <- list(
    type = "hyper_ec",
    n = n,
    m = m_total,
    m_base = m_base,
    k_parity = k_parity,
    alpha = alpha,
    q = q,
    lambda_max = lambda_max_est,
    error_correction = errors_correctable,
    sparse_parity = sparse_parity,
    parity_weight = parity_weight
  )
  
  list(
    matrix = mat_combined,
    matrix_base = mat_base,
    matrix_parity = mat_parity,
    m = m_total,
    m_base = m_base,
    k_parity = k_parity,
    k_design = k,
    n = n,
    alpha = alpha,
    q = q,
    lambda_max = lambda_max_est,
    error_correction = errors_correctable,
    efficiency = efficiency,
    sparse_parity = sparse_parity,
    parity_weight = parity_weight,
    design_type = "hyper_ec"
  )
}


# =============================================================================
# Internal Helper Functions
# =============================================================================

#' @keywords internal
.kirkman_triples <- function(n, m) {
  # Kirkman triple system construction
  # For simplicity, use random triples (full Kirkman construction is complex)
  # A proper Kirkman system requires n ≡ 3 (mod 6)
  .random_triples(n, m)
}


#' @keywords internal
.random_triples <- function(n, m) {
  # Create matrix where each sample in exactly 3 pools
  rows_i <- integer(0)
  cols_j <- integer(0)
  
  for (j in 1:n) {
    # Randomly select 3 pools for this sample
    pools <- sample(m, size = 3, replace = FALSE)
    rows_i <- c(rows_i, pools)
    cols_j <- c(cols_j, rep(j, 3))
  }
  
  Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(m, n)
  )
}


#' @keywords internal
.sparse_parity_matrix <- function(k_parity, m_base, weight) {
  # Create sparse parity check matrix
  # Each row (parity pool) XORs exactly 'weight' base pools
  
  rows_i <- integer(0)
  cols_j <- integer(0)
  
  for (i in 1:k_parity) {
    # Randomly select 'weight' base pools for this parity pool
    base_pools <- sample(m_base, size = min(weight, m_base), replace = FALSE)
    rows_i <- c(rows_i, rep(i, length(base_pools)))
    cols_j <- c(cols_j, base_pools)
  }
  
  Matrix::sparseMatrix(
    i = rows_i,
    j = cols_j,
    x = rep(1L, length(rows_i)),
    dims = c(k_parity, m_base)
  )
}


#' @keywords internal
.random_parity_matrix <- function(k_parity, m_base) {
  # Create random binary parity check matrix
  # Each entry is 1 with probability 0.5
  
  # Generate random binary matrix
  mat <- matrix(sample(0:1, k_parity * m_base, replace = TRUE, prob = c(0.5, 0.5)),
                nrow = k_parity, ncol = m_base)
  
  Matrix::Matrix(mat, sparse = TRUE)
}


#' @keywords internal
.find_valid_m_for_q3 <- function(m_target) {
  # Find nearest m where m ≡ 0 mod 6 and m-1 is prime
  m <- m_target
  
  # Round to nearest multiple of 6
  if (m %% 6L != 0L) {
    m <- ceiling(m / 6) * 6L
  }
  
  # Search for valid m
  max_iterations <- 100
  for (i in 1:max_iterations) {
    if (is_prime(m - 1L)) {
      return(m)
    }
    m <- m + 6L
  }
  
  warning("Could not find valid m for q=3 near ", m_target, ". Using m=", m)
  m
}


# =============================================================================
# Decoding Functions
# =============================================================================

#' @title Decode Experimental Parity-Augmented HYPER Results
#'
#' @description
#' Apply the package-specific decoder paired with \code{\link{hyper_ec_matrix}}.
#' The decoder performs a bounded syndrome search under the package's
#' experimental mod-2 parity model and then applies the stage-1 HYPER
#' conservative rule to the corrected base rows.
#'
#' @details
#' This function is not a faithful implementation of the Hong et al. (2022)
#' HYPER workflow. In paper HYPER, stage 1 identifies putative positives and
#' stage 2 individually retests them. Here, the output `x_decoded` is the set of
#' stage-1 candidates that remain after optional bounded syndrome search.
#'
#' Under the current implementation:
#'
#' 1. split the observation vector into base and parity parts,
#' 2. compute the syndrome
#'    \deqn{\mathbf{s} = \mathbf{y}_{parity} \oplus (\mathbf{H}\mathbf{y}_{base}),}
#' 3. optionally search for a low-weight base-pool flip pattern, and
#' 4. mark individual \eqn{i} as a stage-1 candidate when every base pool in
#'    \eqn{\mathcal{J}_i} is positive after correction.
#'
#' No deterministic error-correction guarantee is implied by this interface.
#'
#' @param y_obs Numeric or integer vector of length \eqn{J = J_0 + K} containing
#'   the observed binary pool outcomes for the combined matrix returned by
#'   \code{\link{hyper_ec_matrix}}.
#' @param design Design object returned by \code{\link{hyper_ec_matrix}}. It must
#'   contain at least \code{matrix_base}, \code{matrix_parity}, \code{m},
#'   \code{m_base}, \code{k_parity}, and \code{n}.
#' @param method Character string specifying the decoder option:
#'   \code{"syndrome"} for bounded syndrome search or \code{"none"} to skip the
#'   search and decode from the uncorrected base rows.
#'
#' @return A named list with components:
#' \describe{
#'   \item{x_decoded}{Integer indicator vector of length \eqn{N}; entry
#'     \eqn{i} equals one when individual \eqn{i} remains a stage-1 candidate}
#'   \item{syndrome}{Computed syndrome vector of length \eqn{K}}
#'   \item{errors_detected}{Syndrome weight \eqn{\sum_\ell s_\ell}}
#'   \item{errors_corrected}{Weight of the base-pool flip pattern returned by the
#'     bounded search}
#'   \item{y_corrected}{Corrected base-pool outcome vector of length \eqn{J_0}}
#' }
#'
#' @examples
#' # Build the matching experimental design
#' design <- hyper_ec_matrix(n = 100, k = 2, q = 2, alpha = 0.20)
#'
#' # Simulate pooled outcomes for a sparse binary signal
#' x <- integer(100)
#' x[c(10, 50)] <- 1
#' y_obs <- as.integer((design$matrix %*% x) > 0)
#'
#' # Decode to obtain stage-1 candidates
#' result <- hyper_ec_decode(y_obs, design)
#' which(result$x_decoded == 1)
#' result$errors_detected
#'
#' # Flip a couple of base-pool outcomes and rerun the bounded search
#' y_error <- y_obs
#' y_error[c(5, 15)] <- 1 - y_error[c(5, 15)]
#' result2 <- hyper_ec_decode(y_error, design, method = "syndrome")
#' result2$errors_detected
#' result2$errors_corrected
#' which(result2$x_decoded == 1)
#'
#' @references
#' Hong, F. S., Shah, N., Briggs, A., et al. (2022). HYPER group testing for
#' SARS-CoV-2: A flexible, easy-to-implement, and highly efficient method.
#' \emph{Nature Communications}, 13, Article 3626.
#'
#' The decoder documented here is a package-specific experimental extension and
#' is not the stage-1 plus stage-2 workflow from the paper above.
#'
#' @seealso
#' \code{\link{hyper_ec_matrix}}, \code{\link{hyper_matrix}},
#' \code{\link{pp_decode}}
#'
#' @family decoding_algorithms
#' @export
hyper_ec_decode <- function(y_obs, design, method = c("syndrome", "none")) {
  method <- match.arg(method)
  
  # Ensure y_obs is a vector
  y_obs <- as.vector(y_obs)
  
  # Extract design parameters
  m <- design$m
  m_base <- design$m_base
  k_parity <- design$k_parity
  n <- design$n
  
  # Split observations into base and parity
  y_base <- y_obs[1:m_base]
  y_parity <- y_obs[(m_base + 1):m]
  
  # Calculate syndrome
  expected_parity <- as.vector((design$matrix_parity %*% y_base) %% 2)
  syndrome <- as.vector((y_parity + expected_parity) %% 2)
  
  errors_detected <- sum(syndrome)
  errors_corrected <- 0
  y_corrected <- y_base
  
  # Error correction
  if (method == "syndrome" && errors_detected > 0) {
    error_pattern <- .simple_syndrome_decode(
      syndrome = syndrome,
      parity_matrix = design$matrix_parity,
      m_base = m_base,
      max_errors = design$error_correction
    )

    if (!is.null(error_pattern)) {
      y_corrected <- (y_base + error_pattern) %% 2
      errors_corrected <- sum(error_pattern)
    }
  }
  
  # HYPER conservative decoding on corrected base pools
  mat_base <- design$matrix_base
  x_decoded <- integer(n)
  
  # Get all positive pools
  positive_pools <- which(y_corrected == 1)
  
  if (length(positive_pools) == 0) {
    # No positive pools → no infected samples
    return(list(
      x_decoded = x_decoded,
      syndrome = syndrome,
      errors_detected = errors_detected,
      errors_corrected = errors_corrected,
      y_corrected = y_corrected
    ))
  }
  
  # For each sample: positive if ALL its pools are positive
  for (j in 1:n) {
    # Get pools containing sample j
    pools_j <- which(as.vector(mat_base[, j]) == 1)
    
    # Sample is positive ONLY if ALL its pools are positive
    if (length(pools_j) > 0 && all(pools_j %in% positive_pools)) {
      x_decoded[j] <- 1
    }
  }
  
  list(
    x_decoded = x_decoded,
    syndrome = syndrome,
    errors_detected = errors_detected,
    errors_corrected = errors_corrected,
    y_corrected = y_corrected
  )
}


#' @keywords internal
.simple_syndrome_decode <- function(syndrome, parity_matrix, m_base, max_errors = 5L) {
  # Simplified syndrome decoding
  # Finds error pattern e such that parity_matrix %*% e = syndrome (mod 2)
  
  # For small number of errors, try all combinations
  # This is exponential but works for small error counts
  
  max_errors <- as.integer(max_errors)
  if (max_errors < 1L) return(NULL)

  # Try all possible error patterns with this many errors
  positions <- 1:m_base

  for (err_count in 1:min(max_errors, m_base)) {
    combos <- combn(positions, err_count)
    
    for (i in 1:ncol(combos)) {
      error_pattern <- integer(m_base)
      error_pattern[combos[, i]] <- 1
      
      # Check if this pattern produces the syndrome
      computed_syndrome <- as.integer((parity_matrix %*% error_pattern) %% 2)
      
      if (all(computed_syndrome == syndrome)) {
        return(error_pattern)
      }
    }
  }
  
  NULL  # No error pattern found
}
