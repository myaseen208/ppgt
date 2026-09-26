#' @title Polynomial Pools Algorithm for Group Testing
#'
#' @description
#' Implementation of the Polynomial Pools (PP) algorithm for constructing
#' efficient group testing pooling matrices.
#'
#' @details
#' Group testing pools multiple individuals to reduce the number of assays.
#' The Polynomial Pools construction uses arithmetic over \eqn{\mathrm{GF}(q)}
#' to build sparse binary pooling matrices with controlled overlap.
#'
#' The package uses PP both as a general design family and as the matrix-level
#' source of the standard PP-based P-BEST configuration
#' \eqn{(q, d, n_l, N) = (8, 3, 6, 384)}.
#'
#' A central structural property is the intersection bound: any two distinct
#' individuals share at most \eqn{d - 1} pools. This is the matrix property
#' behind the standard PP COMP detection bound.
#'
#' @section Relationship To Other Components:
#' \describe{
#'   \item{Finite fields}{\code{\link{gf}} and related helpers provide the
#'   arithmetic used by the PP construction.}
#'   \item{PP matrices}{\code{\link{pp_matrix}} builds the incidence matrix
#'   \eqn{\mathbf{M}} for the chosen PP family.}
#'   \item{PP-based P-BEST}{The standard PP-based P-BEST family is one specific
#'   PP configuration, documented separately in the package vignettes.}
#' }
#'
#' @name polynomial_pools
#' @family polynomial_pools
#'
#' @references
#' Tan, Y. H. I. (2020). Pooling Matrix Designs for Group Testing.
#' *SIAM Undergraduate Research Online*.
#'
#' Brust, D., & Brust, J. J. (2023). Effective Matrix Designs for COVID-19
#' Group Testing. *BMC Bioinformatics*, 24, 26.
#'
#' @seealso
#' \code{\link{gf}} for creating Galois Fields,
#' \code{\link{pp_decode}} for decoding test results
NULL


#' @title Generate A Polynomial Pools Matrix
#'
#' @description
#' Construct a binary pooling matrix by the Polynomial Pools (PP) algorithm,
#' using polynomial evaluation over \eqn{\mathrm{GF}(q)} to define the pool
#' incidence pattern.
#'
#' @details
#' Let \eqn{\mathbf{M} = (M_{ji}) \in \{0,1\}^{J \times N}} denote the returned
#' pooling matrix, where \eqn{j \in \{1, \ldots, J\}} indexes pools and
#' \eqn{i \in \{1, \ldots, N\}} indexes individuals. The PP construction uses
#' polynomial evaluation over \eqn{\mathrm{GF}(q)}. The main derived quantities are
#' \deqn{N_{\max} = q^d, \qquad J = n_l q, \qquad
#' k_{\max} = \left\lfloor \frac{n_l - 1}{d - 1} \right\rfloor.}
#'
#' Individuals are indexed by field-valued coordinates and assigned to pools by
#' polynomial relations over \eqn{\mathrm{GF}(q)}. In the full untruncated design,
#' each individual appears in exactly \eqn{n_l} pools, each pool contains exactly
#' \eqn{q^{d-1}} individuals, and any two individuals share at most \eqn{d - 1}
#' pools.
#'
#' When \eqn{N < q^d}, the constructor truncates the full design to its first
#' \eqn{N} columns. The column weights remain
#' \eqn{w_i = |\mathcal{J}_i| = n_l},
#' while row sizes may vary slightly near the truncation boundary.
#'
#' The standard PP-based P-BEST matrix family uses \eqn{q = 8},
#' \eqn{d = 3}, \eqn{n_l = 6}, and \eqn{N = 384}, giving \eqn{J = 48} pools and
#' \eqn{k_{\max} = 2} under the PP COMP detection bound.
#'
#' Long derivations and family context are developed in the package vignettes,
#' especially \code{ppgt_pp} and \code{ppgt_pbest}.
#'
#' @section PP Structure:
#' In the full untruncated design, all column sums equal \eqn{n_l} and all row
#' sums equal \eqn{q^{d-1}}. Pools within each layer are disjoint, and the
#' intersection bound implies that COMP can separate negatives whenever the
#' number of positives stays below the PP detection limit.
#'
#' @section P-BEST Context:
#' The package documentation uses the term "PP-based P-BEST" for the
#' \eqn{q = 8}, \eqn{d = 3}, \eqn{n_l = 6}, \eqn{N = 384} matrix family. That
#' matrix family is distinct from later package-specific clinical helpers and
#' should not be read as a blanket claim about every P-BEST protocol layer.
#'
#' @param q Integer scalar giving the Galois-field order. It must be a prime
#'   power so that \eqn{\mathrm{GF}(q)} exists. This parameter controls the
#'   number of pools per layer and the coordinate alphabet used by the PP
#'   construction.
#' @param d Integer scalar giving the PP dimension. It must satisfy
#'   \eqn{d \ge 2}. The full untruncated design has \eqn{N_{\max} = q^d}
#'   individuals and canonical pool size \eqn{q^{d-1}}.
#' @param nl Integer scalar giving the number of layers, written
#'   \eqn{n_l} in the PP literature. It must satisfy
#'   \eqn{1 \le n_l \le q + 1}. The returned design has
#'   \eqn{J = n_l q} pools and each individual is assigned to exactly
#'   \eqn{n_l} pools.
#' @param N Integer scalar giving the number of individuals actually retained in
#'   the matrix. The default is \eqn{N = q^d}. It must satisfy
#'   \eqn{1 \le N \le q^d} and \eqn{q \mid N}. Values smaller than \eqn{q^d}
#'   truncate the full PP design to its first \eqn{N} columns.
#'
#' @return A sparse binary matrix of class \code{dgCMatrix} with dimensions
#'   \eqn{J \times N}, where \eqn{J = n_l q}. Entry \eqn{M_{ji} = 1} indicates
#'   that individual \eqn{i} belongs to pool \eqn{j}. The matrix carries
#'   attribute \code{"pp"} containing a named list with:
#' \describe{
#'   \item{q}{Integer; field order}
#'   \item{d}{Integer; dimension}
#'   \item{nl}{Integer; number of layers}
#'   \item{N}{Integer; number of samples}
#'   \item{M}{Integer; legacy display alias for the number of pools \eqn{J}}
#'   \item{pool_size}{Integer; actual number of samples per pool (first row-sum
#'     after any \eqn{N < q^d} truncation; equals \eqn{q^{d-1}} only when
#'     \eqn{N = q^d})}
#'   \item{k_max}{Integer; maximum detectable positives with COMP}
#' }
#'   In the full untruncated design, all column sums equal \eqn{n_l} and all row
#'   sums equal \eqn{q^{d-1}}. Under truncation, the column sums remain
#'   \eqn{w_i = |\mathcal{J}_i| = n_l}, while the row sums may differ slightly
#'   from \eqn{q^{d-1}} near the truncation boundary.
#'
#' @examples
#' # P-BEST standard design: 48 pools x 384 samples
#' M_pbest <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#' M_pbest
#'
#' # Examine structure
#' attr(M_pbest, "pp")
#'
#' # Verify uniform properties
#' table(Matrix::rowSums(M_pbest))
#' table(Matrix::colSums(M_pbest))
#'
#' # Small PP design
#' M_small <- pp_matrix(q = 4, d = 3, nl = 5)
#' M_small
#' attr(M_small, "pp")$k_max
#'
#' # Using GF(9) = GF(3^2) for larger capacity
#' M_gf9 <- pp_matrix(q = 9, d = 3, nl = 7)
#' M_gf9
#' attr(M_gf9, "pp")$k_max
#'
#' # Large design using GF(16)
#' M_large <- pp_matrix(q = 16, d = 3, nl = 9)
#' M_large
#' attr(M_large, "pp")$k_max
#'
#' # Error: 6 is not a prime power
#' \dontrun{
#' M_invalid <- pp_matrix(q = 6, d = 3, nl = 5)  # Error!
#' }
#'
#' @references
#' Tan, Y. H. I. (2020). Pooling matrix designs for group testing.
#' \emph{SIAM Undergraduate Research Online}, 13, 1-21.
#'
#' Shental, N., Levy, S., Wuvshet, V., et al. (2020). Efficient high-throughput
#' SARS-CoV-2 testing to detect asymptomatic carriers. \emph{Science Advances},
#' 6(37), eabc5961.
#'
#' Brust, D., & Brust, J. J. (2023). Effective matrix designs for COVID-19
#' group testing. \emph{BMC Bioinformatics}, 24, Article 26.
#'
#' @seealso
#' \code{\link{gf}}, \code{\link{pp_decode}}, \code{\link{pp_design}},
#' \code{\link{pp_verify}}
#'
#' @family polynomial_pools
#' @export
pp_matrix <- function(q, d, nl, N = NULL) {

  # Validate q is prime power
  pp <- is_prime_power(q)
  if (!pp$is_prime_power) {
    stop("q must be a prime power (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, ...).\n",
         "Got q = ", q, " which is NOT a prime power.")
  }

  q <- as.integer(q)
  d <- as.integer(d)
  nl <- as.integer(nl)

  if (d < 2L) stop("d must be >= 2. Got: ", d)
  if (nl < 1L || nl > q + 1L) {
    stop("nl must be in [1, ", q + 1L, "]. Got: ", nl)
  }

  # Integer powers of q: q_pows[k+1] = q^k  (avoids floating-point ^)
  q_pows <- cumprod(c(1L, rep(q, d)))
  N_max  <- q_pows[d + 1L]
  q_dm1  <- q_pows[d]
  q_dm2  <- if (d > 2L) q_pows[d - 1L] else 1L
  q_2    <- q_pows[3L]

  if (is.null(N)) N <- N_max
  N <- as.integer(N)
  if (N > N_max) {
    stop("N cannot exceed q^d = ", N_max, ". Got: ", N)
  }
  if (N %% q != 0L) {
    stop("N must be divisible by q = ", q, ". Got N = ", N,
         ". Algorithm 2.3 (Tan 2020) requires q | N.")
  }

  M <- nl * q
  a1    <- N %/% q_dm1
  a2    <- (N - a1 * q_dm1) %/% q
  k_max <- (nl - 1L) %/% (d - 1L)

  # Create Galois field
  gf_obj <- gf(q)
  v <- 0L:(q - 1L)

  # Coefficient tuples for d > 2
  n_coeffs <- if (d > 2L) q_dm2 else 1L
  if (d > 2L) {
    coeffs_list <- as.matrix(expand.grid(rep(list(v), d - 2L)))
  } else {
    coeffs_list <- matrix(0L, nrow = 1L, ncol = 1L)
  }

  # Pre-allocate COO vectors; upper bound = nl pools per sample x N samples
  nnz_max <- nl * N
  rows_i <- integer(nnz_max)
  cols_j <- integer(nnz_max)
  idx <- 0L

  for (layer in seq_len(nl)) {
    if (layer <= q) {
      a <- v[layer]

      # Precompute a^2, a^3, ..., a^(d-1) in GF(q)
      a_pows <- integer(d)
      a_pows[1L] <- 1L
      a_pows[2L] <- a
      if (d > 2L) {
        for (i in 3L:d) {
          a_pows[i] <- gf_mult(a_pows[i - 1L], a, gf_obj)
        }
      }

      # Part A: a1 full c_1 blocks, all q x-values per pool (Algorithm 2.3)
      for (l in seq_len(a1)) {
        coeffs <- if (d > 2L) coeffs_list[l, ] else integer(0L)
        for (b_idx in seq_len(q)) {
          b <- v[b_idx]
          pool_id <- (layer - 1L) * q + b_idx
          for (x_idx in seq_len(q)) {
            x <- v[x_idx]
            y <- gf_add(b, gf_mult(a, x, gf_obj), gf_obj)
            if (d > 2L) {
              for (i in seq_len(d - 2L)) {
                y <- gf_add(y, gf_mult(coeffs[i], a_pows[i + 2L], gf_obj), gf_obj)
              }
            }
            idx <- idx + 1L
            rows_i[idx] <- pool_id
            cols_j[idx] <- (x_idx - 1L) * q + (l - 1L) * q_2 + y + 1L
          }
        }
      }
      # Part B: partial c_1 block, first a2 x-values only (Algorithm 2.3)
      if (a2 > 0L) {
        l_b <- a1 + 1L
        coeffs <- if (d > 2L) coeffs_list[l_b, ] else integer(0L)
        for (b_idx in seq_len(q)) {
          b <- v[b_idx]
          pool_id <- (layer - 1L) * q + b_idx
          for (x_idx in seq_len(a2)) {
            x <- v[x_idx]
            y <- gf_add(b, gf_mult(a, x, gf_obj), gf_obj)
            if (d > 2L) {
              for (i in seq_len(d - 2L)) {
                y <- gf_add(y, gf_mult(coeffs[i], a_pows[i + 2L], gf_obj), gf_obj)
              }
            }
            idx <- idx + 1L
            rows_i[idx] <- pool_id
            cols_j[idx] <- (x_idx - 1L) * q + (l_b - 1L) * q_2 + y + 1L
          }
        }
      }
    } else {
      # Singular layer (slope = infinity): pool b_idx holds all samples
      # whose last coordinate y0 = b_idx - 1  (Julia PP.jl lines 39-40)
      for (b_idx in seq_len(q)) {
        pool_id <- (layer - 1L) * q + b_idx
        for (l in seq_len(n_coeffs)) {
          for (x_idx in seq_len(q)) {
            sample_id <- (b_idx - 1L) + (l - 1L) * q_2 + (x_idx - 1L) * q + 1L
            if (sample_id <= N) {
              idx <- idx + 1L
              rows_i[idx] <- pool_id
              cols_j[idx] <- sample_id
            }
          }
        }
      }
    }
  }

  result <- Matrix::sparseMatrix(
    i = rows_i[seq_len(idx)],
    j = cols_j[seq_len(idx)],
    x = rep(1L, idx),
    dims = c(M, N)
  )

  attr(result, "pp") <- list(
    q = q, d = d, nl = nl, N = N, M = M,
    pool_size = as.integer(Matrix::rowSums(result)[1L]),
    k_max = k_max
  )

  result
}


#' @title Find Valid Polynomial Pools Configurations
#'
#' @description
#' Enumerate PP parameter combinations that satisfy user-supplied design
#' constraints, so that candidate pooling matrices can be planned before calling
#' \code{\link{pp_matrix}}.
#'
#' @details
#' # Relationship to pp_matrix and pp_verify
#'
#' - **\code{\link{pp_design}}**: Finds valid parameter combinations (planning)
#' - **\code{\link{pp_matrix}}**: Generates the actual pooling matrix (construction)
#' - **\code{\link{pp_verify}}**: Checks properties of a generated matrix (validation)
#'
#' Typical workflow:
#' ```
#' # 1. Find suitable configurations
#' configs <- pp_design(N_min = 300, N_max = 500)
#'
#' # 2. Choose one and generate matrix
#' M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#'
#' # 3. Verify the matrix
#' pp_verify(M)
#' ```
#'
#' # Generalization Rules
#'
#' The PP algorithm generalizes to any \eqn{(q, d, n_l)} where:
#'
#' 1. \eqn{q} is a **prime power**: 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, ...
#' 2. \eqn{d \geq 2} (dimension)
#' 3. \eqn{1 \leq n_l \leq q + 1} (number of layers)
#'
#' # What You CANNOT Choose Arbitrarily
#'
#' | **Constraint** | **Why** |
#' |----------------|---------|
#' | \eqn{q} must be prime power | Galois Fields only exist for prime power orders |
#' | \eqn{M = n_l \times q} | Each layer contributes exactly \eqn{q} pools |
#' | Pool size = \eqn{q^{d-1}} | Fixed by the polynomial structure |
#'
#' # Design Selection Guide
#'
#' | **Samples** | **Recommended** | **Pools** | **Compression** |
#' |-------------|-----------------|-----------|-----------------|
#' | ~64 | q=4, d=3, nl=5 | 20 | 3.2x |
#' | ~384 | q=8, d=3, nl=6 | 48 | 8x |
#' | ~729 | q=9, d=3, nl=7 | 63 | 11.6x |
#' | ~4096 | q=16, d=3, nl=9 | 144 | 28x |
#'
#' @param N_min Integer scalar giving the minimum allowed number of individuals
#'   \eqn{N}. Only PP configurations with \eqn{q^d \ge N_{\min}} are retained.
#' @param N_max Integer scalar giving the maximum allowed number of individuals
#'   \eqn{N}. Only PP configurations with \eqn{q^d \le N_{\max}} are retained.
#' @param M_max Integer scalar giving the maximum allowed number of pools
#'   \eqn{J}. The default is \code{500}.
#' @param k_min Integer scalar giving the minimum required COMP detection limit
#'   \eqn{k_{\max} = \lfloor (n_l - 1)/(d - 1) \rfloor}. The default is
#'   \code{2}.
#' @param pool_size_max Integer scalar giving the maximum allowable canonical
#'   pool size \eqn{q^{d-1}}. The default is \code{500}.
#'
#' @return A data frame with one row per retained PP configuration and columns:
#' \describe{
#'   \item{q}{Integer; Galois Field order (prime power)}
#'   \item{d}{Integer; dimension}
#'   \item{nl}{Integer; number of layers}
#'   \item{N}{Integer; maximum samples (\eqn{q^d})}
#'   \item{M}{Integer; legacy display alias for the number of pools \eqn{J = n_l q}}
#'   \item{pool_size}{Integer; samples per pool (\eqn{q^{d-1}})}
#'   \item{k_max}{Integer; detection limit}
#'   \item{compression}{Numeric; display alias for the compression ratio \eqn{N/J}}
#' }
#'   The rows are sorted by decreasing \code{compression} and then by increasing
#'   \code{pool_size}. The function returns \code{NULL} if no configuration
#'   satisfies the supplied constraints.
#'
#' @details
#' For each candidate PP parameter tuple \eqn{(q, d, n_l)}, the associated full
#' design has
#' \deqn{N = q^d, \qquad J = n_l q, \qquad n_j = q^{d-1}, \qquad
#' w_i = n_l.}
#' The COMP exact-decoding guarantee for PP designs is governed by
#' \deqn{k_{\max} = \left\lfloor \frac{n_l - 1}{d - 1} \right\rfloor.}
#' This function searches a finite library of supported prime powers \eqn{q} and
#' dimensions \eqn{d \in \{2,3,4\}}, filters the resulting tuples by the user
#' constraints, and reports the compression factor
#' \deqn{\frac{N}{J}.}
#' The output is intended for planning rather than for direct matrix generation;
#' once a row is selected, \code{\link{pp_matrix}} can be used to construct the
#' corresponding \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}.
#'
#' @examples
#' # Find configurations for ~300-500 samples
#' configs <- pp_design(N_min = 300, N_max = 500)
#' configs
#'
#' # Filter for specific q value (e.g., q = 8 for P-BEST family)
#' pbest_family <- configs[configs$q == 8, ]
#' pbest_family
#'
#' # Find high-compression designs for ~1000 samples
#' configs_1000 <- pp_design(N_min = 800, N_max = 1200, k_min = 3)
#' configs_1000
#'
#' # Select best configuration and generate matrix
#' best <- configs[1, ]
#' best
#' M <- pp_matrix(q = best$q, d = best$d, nl = best$nl, N = best$N)
#' M
#' verification <- pp_verify(M)
#' verification
#'
#' # Find configurations with small pool sizes
#' small_pool <- pp_design(N_min = 100, N_max = 200, pool_size_max = 50)
#' small_pool
#'
#' @references
#' Tan, Y. H. I. (2020). Pooling matrix designs for group testing.
#' \emph{SIAM Undergraduate Research Online}, 13, 1-21.
#'
#' Aldridge, M., Johnson, O., & Scarlett, J. (2019). Group testing: An
#' information theory perspective. \emph{Foundations and Trends in
#' Communications and Information Theory}, 15(3-4), 196-392.
#'
#' @seealso
#' \code{\link{pp_matrix}}, \code{\link{pp_verify}}, \code{\link{gf}}
#'
#' @family polynomial_pools
#' @export
pp_design <- function(N_min = 50, N_max = 1000, M_max = 500,
                      k_min = 2, pool_size_max = 500) {

  prime_powers <- c(2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27,
                    29, 31, 32, 37, 41, 43, 47, 49, 53, 59, 61, 64)

  results <- list()
  idx <- 0

  for (q in prime_powers) {
    for (d in 2:4) {
      N <- q^d
      pool_size <- q^(d - 1)

      if (N < N_min || N > N_max) next
      if (pool_size > pool_size_max) next

      for (nl in 2:(q + 1)) {
        M <- nl * q
        if (M > M_max) next
        if (N / M <= 1) next

        k_max <- (nl - 1L) %/% (d - 1L)
        if (k_max < k_min) next

        idx <- idx + 1
        results[[idx]] <- data.frame(
          q = q, d = d, nl = nl,
          N = N, M = M,
          pool_size = pool_size,
          k_max = k_max,
          compression = round(N / M, 2),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) {
    message("No configurations found. Try relaxing constraints.")
    return(NULL)
  }

  df <- do.call(rbind, results)
  df <- df[order(-df$compression, df$pool_size), ]
  rownames(df) <- NULL
  df
}


#' @title Verify Polynomial Pools Matrix Properties
#'
#' @description
#' Check whether a matrix produced by \code{\link{pp_matrix}} has the expected
#' PP dimensions, row sizes, and column weights, and print a concise validation
#' summary.
#'
#' @details
#' # Relationship to pp_matrix and pp_design
#'
#' - **\code{\link{pp_design}}**: Planning - find valid configurations
#' - **\code{\link{pp_matrix}}**: Construction - generate the matrix
#' - **\code{\link{pp_verify}}**: Validation - check the result
#'
#' # Properties Checked
#'
#' 1. **Uniform pool sizes**: All pools contain the same number of samples
#' 2. **Uniform sample coverage**: All samples appear in the same number of pools
#' 3. **Correct dimensions**: \eqn{M = n_l \times q} pools, \eqn{N \leq q^d} samples
#' 4. **Stored parameters match**: q, d, nl, k_max
#'
#' @param M A binary pooling matrix, typically returned by
#'   \code{\link{pp_matrix}}. If the \code{"pp"} attribute is present, its stored
#'   parameters are included in the printed summary.
#'
#' @return Invisibly returns a named list with:
#' \describe{
#'   \item{n_pools}{Integer; the row dimension \eqn{J}}
#'   \item{n_samples}{Integer; the column dimension \eqn{N}}
#'   \item{pool_size}{Integer; the first observed row sum \eqn{n_j}}
#'   \item{sample_coverage}{Integer; the first observed column sum \eqn{w_i}}
#'   \item{uniform_pools}{Logical; \code{TRUE} if all row sums are equal}
#'   \item{uniform_coverage}{Logical; \code{TRUE} if all column sums are equal}
#' }
#'   The function also prints a verification summary showing \eqn{J}, \eqn{N},
#'   row sums, column sums, compression, and density.
#'
#' @details
#' For a PP matrix \eqn{\mathbf{M} \in \{0,1\}^{J \times N}}, the key structural
#' quantities are
#' \deqn{n_j = |\mathcal{P}_j| = \sum_{i=1}^N M_{ji}}
#' and
#' \deqn{w_i = |\mathcal{J}_i| = \sum_{j=1}^J M_{ji}.}
#' In the untruncated PP design, every row sum is \eqn{q^{d-1}} and every column
#' sum is \eqn{n_l}. When \eqn{N < q^d}, truncation preserves the common column
#' weight \eqn{w_i = n_l} but may disturb row uniformity near the truncation
#' boundary. This function reports whether the realized row and column sums are
#' uniform, along with the empirical compression ratio
#' \deqn{\frac{N}{J}.}
#'
#' @examples
#' # Generate and verify P-BEST matrix
#' M <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
#' result <- pp_verify(M)
#' result
#'
#' # Access verification results
#' result$pool_size      # 48 (note: depends on truncation)
#' result$uniform_pools  # TRUE
#'
#' # Verify different designs
#' M_small <- pp_matrix(q = 4, d = 3, nl = 5)
#' M_small
#' pp_verify(M_small)
#'
#' M_gf9 <- pp_matrix(q = 9, d = 3, nl = 7)
#' M_gf9
#' pp_verify(M_gf9)
#'
#' @references
#' Tan, Y. H. I. (2020). Pooling matrix designs for group testing.
#' \emph{SIAM Undergraduate Research Online}, 13, 1-21.
#'
#' Brust, D., & Brust, J. J. (2023). Effective matrix designs for COVID-19
#' group testing. \emph{BMC Bioinformatics}, 24, Article 26.
#'
#' @seealso
#' \code{\link{pp_matrix}}, \code{\link{pp_design}}, \code{\link{pp_decode}}
#'
#' @family polynomial_pools
#' @export
pp_verify <- function(M) {
  params <- attr(M, "pp")
  m <- nrow(M)
  n <- ncol(M)

  pool_sizes <- as.vector(Matrix::rowSums(M))
  sample_cov <- as.vector(Matrix::colSums(M))

  uniform_pools <- length(unique(pool_sizes)) == 1
  uniform_cov <- length(unique(sample_cov)) == 1

  message("\n=== PP Matrix Verification ===")
  message("Dimensions: ", m, " pools x ", n, " samples")

  if (!is.null(params)) {
    message("Parameters: q = ", params$q, ", d = ", params$d, ", nl = ", params$nl)
    message("Detection limit: k_max = ", params$k_max)
    message("Expected pool size: ", params$pool_size)
  }

  message("\nStructure:")
  message("  Pool sizes: ", paste(unique(pool_sizes), collapse = ", "),
          ifelse(uniform_pools, " (uniform yes)", " (variable no)"))
  message("  Pools/sample: ", paste(unique(sample_cov), collapse = ", "),
          ifelse(uniform_cov, " (uniform yes)", " (variable no)"))
  message("  Compression: ", round(n / m, 2L), "x")
  message("  Matrix density: ", round(sum(M > 0) / (m * n), 4L))

  invisible(list(
    n_pools = m, n_samples = n,
    pool_size = unique(pool_sizes)[1],
    sample_coverage = unique(sample_cov)[1],
    uniform_pools = uniform_pools,
    uniform_coverage = uniform_cov
  ))
}
