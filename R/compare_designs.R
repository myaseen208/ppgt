#' @title Compare Built-In Group Testing Designs
#'
#' @description
#' Generate several built-in pooling designs and summarize their structural
#' properties in a common comparison object. The current implementation compares
#' the package constructors for P-BEST, Tapestry, HYPER, and the package-specific
#' experimental parity-augmented \code{hyper_ec_matrix()} extension.
#'
#' @param n Integer scalar giving the target number of individuals. This is the
#'   intended column dimension \eqn{N} for the Tapestry, HYPER, and experimental
#'   HYPER-EC comparison branches. The default is \code{256}.
#' @param prevalence Numeric scalar giving the assumed common prevalence
#'   \eqn{p_i = p} used to define the approximate expected number of positives
#'   \eqn{k = \max\{1, \lfloor Np \rfloor\}} for design constructors that require
#'   an input sparsity level. It must lie in \eqn{(0,1)}; the default is
#'   \code{0.01}.
#' @param alpha Numeric scalar giving the parity-overhead factor for the
#'   package-specific experimental \code{hyper_ec_matrix()} branch. The default
#'   is \code{0.20}. In that branch the number of added parity pools is driven by
#'   \eqn{\alpha}.
#'
#' @return A named list with five components:
#'   \describe{
#'   \item{\code{pbest}}{A list with components \code{matrix} and \code{metrics}
#'   for the P-BEST branch.}
#'   \item{\code{tapestry}}{A list with components \code{matrix} and
#'   \code{metrics} for the Tapestry branch.}
#'   \item{\code{hyper}}{A list with components \code{matrix} and \code{metrics}
#'   for the HYPER branch.}
#'   \item{\code{hyper_ec}}{A list with components \code{matrix} and
#'   \code{metrics} for the experimental parity-augmented HYPER-EC branch.}
#'   \item{\code{comparison_table}}{A data frame with one row per design and
#'   columns \code{Design}, \code{Samples_n}, \code{Pools_m},
#'   \code{Efficiency_m_n}, \code{Lambda_max}, and \code{Error_Correction}.}
#'   }
#'   Each \code{metrics} component is itself a named list. For successful design
#'   generation it contains summary statistics such as the achieved number of
#'   individuals \eqn{N}, the number of pools \eqn{J}, the ratio \eqn{J/N}, the
#'   maximum column sum \eqn{\max_i |\mathcal{J}_i|}, and row-size summaries
#'   based on \eqn{n_j = |\mathcal{P}_j|}. If a branch fails, its metrics entry
#'   contains a design label and a failure note. The
#'   \code{Error_Correction} field is a package-specific comparison label for
#'   this helper; it is not a formal paper-backed guarantee.
#'
#' @details
#' Let \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denote a pooling matrix, where
#' \eqn{N} is the number of individuals and \eqn{J} is the number of pools. This
#' function builds one design object for each supported branch and reports a
#' common set of structural summaries. For every successfully generated design,
#' the comparison metrics are derived from
#' \deqn{n_j = |\mathcal{P}_j| = \sum_{i=1}^N M_{ji}}
#' and
#' \deqn{w_i = |\mathcal{J}_i| = \sum_{j=1}^J M_{ji}.}
#' The reported efficiency ratio is
#' \deqn{\frac{J}{N},}
#' which is stored in the output as \code{Efficiency_m_n}. The reported
#' \code{Lambda_max} is a display label for the maximum number of pools
#' containing any single individual,
#' \deqn{\lambda_{\max} = \max_{1 \le i \le N} |\mathcal{J}_i| = \max_i w_i.}
#' The row-size summaries are the empirical mean and standard deviation of
#' \eqn{n_j}.
#'
#' The comparison procedure uses the following package-specific
#' constructor choices.
#' For Tapestry and the experimental HYPER-EC branch, the expected number of
#' positives is approximated by
#' \deqn{k = \max\left\{1, \left\lfloor N p \right\rfloor\right\}.}
#' For HYPER, the baseline number of pools is set heuristically to
#' \deqn{J_{\mathrm{base}} = \max\left\{\left\lceil 4 k \log(N) \right\rceil,
#' \left\lceil 0.15 N \right\rceil\right\},}
#' followed by an adjustment to the nearest even integer so that the
#' \code{q = 2} constructor call is admissible. The P-BEST branch always uses the
#' package's fixed standard constructor call
#' \code{pp_matrix(q = 8, d = 3, nl = 6, N = 384)}, so its achieved column
#' dimension is not forced to equal the input \code{n}. The experimental
#' HYPER-EC branch uses \code{hyper_ec_matrix(n = n, k = k, m_base = m_hyper,
#' q = 2, alpha = alpha)} and should be interpreted as a package-specific
#' extension rather than the Hong et al. HYPER paper design.
#'
#' @examples
#' comparison <- compare_all_designs(n = 128, prevalence = 0.01, alpha = 0.20)
#' comparison
#' names(comparison)
#' comparison$comparison_table
#'
#' comparison$hyper$metrics$lambda_max
#' comparison$tapestry$metrics$pool_size_mean
#'
#' @references
#' Shental, N., Levy, S., Wuvshet, V., et al. (2020). Efficient high-throughput
#' SARS-CoV-2 testing to detect asymptomatic carriers. \emph{Science Advances},
#' 6(37), eabc5961.
#'
#' Ghosh, S., Agarwal, R., Rehan, M. A., et al. (2021). Tapestry: A single-round
#' smart pooling technique for COVID-19 testing. \emph{Nature Communications},
#' 12, Article 2995.
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
#' \code{\link{pp_matrix}}, \code{\link{tapestry_matrix}},
#' \code{\link{hyper_matrix}}, \code{\link{hyper_ec_matrix}},
#' \code{\link{compare_designs}}
#'
#' @family comparison_functions
#' @export
compare_all_designs <- function(n = 256, prevalence = 0.01, alpha = 0.20) {
  n <- as.integer(n)
  
  # Expected infected
  k <- max(1L, as.integer(n * prevalence))
  
  # =============================
  # 1. P-BEST Design
  # =============================
  message("Generating P-BEST design...")
  
  # P-BEST package reference configuration: N=384, M=48 (Shental et al., 2020)
  # Note: q^d = 8^3 = 512, but paper uses N=384
  pbest <- tryCatch({
    pp_matrix(q = 8, d = 3, nl = 6, N = 384)  # ALWAYS use standard config
  }, error = function(e) {
    warning("P-BEST generation failed: ", e$message)
    NULL
  })
  
  if (!is.null(pbest)) {
    pbest_metrics <- list(
      n = ncol(pbest),
      m = nrow(pbest),
      efficiency = nrow(pbest) / ncol(pbest),
      lambda_max = max(Matrix::colSums(pbest)),
      avg_lambda = mean(Matrix::colSums(pbest)),
      pool_size_mean = mean(Matrix::rowSums(pbest)),
      pool_size_sd = sd(Matrix::rowSums(pbest)),
      error_correction = "Strong (Reed-Solomon)",
      design_type = "P-BEST"
    )
  } else {
    pbest_metrics <- list(design_type = "P-BEST", note = "Generation failed")
  }
  
  # =============================
  # 2. Tapestry Design
  # =============================
  message("Generating Tapestry design...")
  
  tapestry <- tryCatch({
    tapestry_matrix(n = n, k = k, method = "random")
  }, error = function(e) {
    warning("Tapestry generation failed: ", e$message)
    NULL
  })
  
  if (!is.null(tapestry)) {
    tapestry_metrics <- list(
      n = tapestry$n,
      m = tapestry$m,
      efficiency = tapestry$m / tapestry$n,
      lambda_max = tapestry$lambda_max,
      avg_lambda = mean(Matrix::colSums(tapestry$matrix)),
      pool_size_mean = mean(Matrix::rowSums(tapestry$matrix)),
      pool_size_sd = sd(Matrix::rowSums(tapestry$matrix)),
      error_correction = "Implicit (Compressed Sensing)",
      design_type = "Tapestry"
    )
  } else {
    tapestry_metrics <- list(design_type = "Tapestry", note = "Generation failed")
  }
  
  # =============================
  # 3. HYPER Design
  # =============================
  message("Generating HYPER design...")
  
  # Determine m for HYPER
  m_hyper <- as.integer(ceiling(4 * k * log(n)))
  m_hyper <- max(m_hyper, as.integer(ceiling(n * 0.15)))
  
  # Adjust for q=2 constraint (m must be even)
  if (m_hyper %% 2L != 0L) m_hyper <- m_hyper + 1L
  
  hyper <- tryCatch({
    hyper_matrix(n = n, m = m_hyper, q = 2)
  }, error = function(e) {
    warning("HYPER generation failed: ", e$message)
    NULL
  })
  
  if (!is.null(hyper)) {
    hyper_metrics <- list(
      n = ncol(hyper),
      m = nrow(hyper),
      efficiency = nrow(hyper) / ncol(hyper),
      lambda_max = max(Matrix::colSums(hyper)),
      avg_lambda = mean(Matrix::colSums(hyper)),
      pool_size_mean = mean(Matrix::rowSums(hyper)),
      pool_size_sd = sd(Matrix::rowSums(hyper)),
      error_correction = "None",
      design_type = "HYPER"
    )
  } else {
    hyper_metrics <- list(design_type = "HYPER", note = "Generation failed")
  }
  
  # =============================
  # 4. HYPER-EC Design
  # =============================
  message("Generating HYPER-EC design...")
  
  hyper_ec <- tryCatch({
    hyper_ec_matrix(n = n, k = k, m_base = m_hyper, q = 2, alpha = alpha)
  }, error = function(e) {
    warning("HYPER-EC generation failed: ", e$message)
    NULL
  })
  
  if (!is.null(hyper_ec)) {
    hyper_ec_metrics <- list(
      n = hyper_ec$n,
      m = hyper_ec$m,
      m_base = hyper_ec$m_base,
      k_parity = hyper_ec$k_parity,
      efficiency = hyper_ec$efficiency,
      lambda_max = hyper_ec$lambda_max,
      avg_lambda = mean(Matrix::colSums(hyper_ec$matrix)),
      pool_size_mean = mean(Matrix::rowSums(hyper_ec$matrix)),
      pool_size_sd = sd(Matrix::rowSums(hyper_ec$matrix)),
      error_correction = paste0("Medium (XOR, corrects ", hyper_ec$error_correction, " errors)"),
      design_type = "HYPER-EC"
    )
  } else {
    hyper_ec_metrics <- list(design_type = "HYPER-EC", note = "Generation failed")
  }
  
  # =============================
  # Create Comparison Table
  # =============================
  comparison_df <- data.frame(
    Design = c("P-BEST", "Tapestry", "HYPER", "HYPER-EC"),
    Samples_n = c(
      ifelse(!is.null(pbest), pbest_metrics$n, NA),
      ifelse(!is.null(tapestry), tapestry_metrics$n, NA),
      ifelse(!is.null(hyper), hyper_metrics$n, NA),
      ifelse(!is.null(hyper_ec), hyper_ec_metrics$n, NA)
    ),
    Pools_m = c(
      ifelse(!is.null(pbest), pbest_metrics$m, NA),
      ifelse(!is.null(tapestry), tapestry_metrics$m, NA),
      ifelse(!is.null(hyper), hyper_metrics$m, NA),
      ifelse(!is.null(hyper_ec), hyper_ec_metrics$m, NA)
    ),
    Efficiency_m_n = c(
      ifelse(!is.null(pbest), round(pbest_metrics$efficiency, 3), NA),
      ifelse(!is.null(tapestry), round(tapestry_metrics$efficiency, 3), NA),
      ifelse(!is.null(hyper), round(hyper_metrics$efficiency, 3), NA),
      ifelse(!is.null(hyper_ec), round(hyper_ec_metrics$efficiency, 3), NA)
    ),
    Lambda_max = c(
      ifelse(!is.null(pbest), pbest_metrics$lambda_max, NA),
      ifelse(!is.null(tapestry), tapestry_metrics$lambda_max, NA),
      ifelse(!is.null(hyper), hyper_metrics$lambda_max, NA),
      ifelse(!is.null(hyper_ec), hyper_ec_metrics$lambda_max, NA)
    ),
    Error_Correction = c(
      ifelse(!is.null(pbest), pbest_metrics$error_correction, "N/A"),
      ifelse(!is.null(tapestry), tapestry_metrics$error_correction, "N/A"),
      ifelse(!is.null(hyper), hyper_metrics$error_correction, "N/A"),
      ifelse(!is.null(hyper_ec), hyper_ec_metrics$error_correction, "N/A")
    ),
    stringsAsFactors = FALSE
  )
  
  message("\nComparison complete!")
  
  list(
    pbest = list(matrix = pbest, metrics = pbest_metrics),
    tapestry = list(matrix = tapestry, metrics = tapestry_metrics),
    hyper = list(matrix = hyper, metrics = hyper_metrics),
    hyper_ec = list(matrix = hyper_ec, metrics = hyper_ec_metrics),
    comparison_table = comparison_df
  )
}


#' @keywords internal
.find_nearest_prime_power <- function(target) {
  # Find nearest prime power to target
  # Check small primes first
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
  
  best_q <- target
  best_dist <- Inf
  
  for (p in primes) {
    # Try powers of this prime
    q <- p
    while (q <= target * 2) {
      dist <- abs(q - target)
      if (dist < best_dist) {
        best_dist <- dist
        best_q <- q
      }
      q <- q * p
    }
  }
  
  best_q
}
