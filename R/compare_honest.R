#' @title Compare Designs At Multiple Reference Sample Sizes
#'
#' @description
#' Generate three side-by-side comparison tables for the package constructors of
#' P-BEST, Tapestry, HYPER, and the package-specific experimental
#' \code{hyper_ec_matrix()} branch. The function reports each design at its
#' chosen reference configuration and at two common sample sizes.
#'
#' @return A named list with four components:
#'   \describe{
#'   \item{\code{standard}}{A data frame comparing the package's chosen standard
#'   constructor settings for each design.}
#'   \item{\code{n256}}{A data frame comparing all four designs at
#'   \eqn{N = 256}.}
#'   \item{\code{n384}}{A data frame comparing all four designs at
#'   \eqn{N = 384}.}
#'   \item{\code{summary}}{A named list of short textual summaries describing the
#'   branch identified by the current implementation as most efficient, lowest
#'   column weight, or most flexible.}
#'   }
#'   Each comparison data frame has four rows, one per design, and columns
#'   describing the achieved number of individuals \eqn{N}, the number of pools
#'   \eqn{J}, the ratio \eqn{J/N}, the maximum column weight
#'   \eqn{\lambda_{\max} = \max_i |\mathcal{J}_i|}, and configuration labels.
#'   The \code{Configuration} column and the entries of \code{summary} are
#'   package-specific labels for this comparison helper; they are not paper
#'   claims or formal optimality guarantees.
#'
#' @details
#' Let \eqn{\mathbf{M} \in \{0,1\}^{J \times N}} denote a pooling matrix, where
#' \eqn{N} is the number of individuals and \eqn{J} is the number of pools. This
#' function reports the comparison metric
#' \deqn{\frac{J}{N}}
#' together with the maximum column weight
#' \deqn{\lambda_{\max} = \max_{1 \le i \le N} |\mathcal{J}_i| = \max_{1 \le i \le
#' N} \sum_{j=1}^J M_{ji}.}
#' These quantities are extracted directly from the matrices returned by the
#' package constructors, or from the corresponding constructor objects when a
#' design returns additional metadata.
#'
#' The function produces three tables. The \code{standard} table compares the
#' specific constructor settings hard-coded as package reference points:
#' P-BEST at \eqn{N = 384}, Tapestry at \eqn{N = 256}, HYPER at \eqn{N = 256},
#' and experimental HYPER-EC at \eqn{N = 256}. The \code{n256} table recomputes
#' all branches at the common target size \eqn{N = 256}, and the \code{n384}
#' table recomputes all branches at the common target size \eqn{N = 384}. The
#' purpose is not to prove a universally optimal design, but to make the
#' dependence of \eqn{J/N} and \eqn{\lambda_{\max}} on the chosen operating size
#' explicit.
#'
#' For a successful branch, the table entry \code{m_over_N} is the realized ratio
#' \eqn{J/N}. When a branch fails to generate, the corresponding table entries are
#' \code{NA}. The experimental HYPER-EC branch should be interpreted as a
#' package-specific extension rather than a faithful implementation of the Hong
#' et al. paper design. Likewise, the printed "summary" text is package-specific
#' reporting for the current helper and should not be interpreted as a literature
#' result.
#'
#' @examples
#' comp <- compare_all_designs_honest()
#' comp
#' names(comp)
#' comp$standard
#' comp$n256
#' comp$n384
#' comp$summary
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
#' \code{\link{compare_all_designs}}, \code{\link{pp_matrix}},
#' \code{\link{tapestry_matrix}}, \code{\link{hyper_matrix}},
#' \code{\link{hyper_ec_matrix}}
#'
#' @export
compare_all_designs_honest <- function() {
  
  message("========================================")
  message("COMPREHENSIVE DESIGN COMPARISON")
  message("========================================\n")
  
  # =============================
  # Comparison 1: Standard Configurations
  # =============================
  message("1. Generating designs at their STANDARD configurations...")
  
  pbest_std <- tryCatch(pp_matrix(q = 8, d = 3, nl = 6, N = 384), error = function(e) NULL)
  tapestry_std <- tryCatch(tapestry_matrix(n = 256, k = 2), error = function(e) NULL)
  hyper_std <- tryCatch(hyper_matrix(n = 256, m = 44, q = 2), error = function(e) NULL)
  hyperec_std <- tryCatch(hyper_ec_matrix(n = 256, q = 2, alpha = 0.20), error = function(e) NULL)
  
  standard_df <- data.frame(
    Design = c("P-BEST", "Tapestry", "HYPER", "HYPER-EC"),
    N_Standard = c(384, 256, 256, 256),
    M = c(
      ifelse(!is.null(pbest_std), nrow(pbest_std), NA),
      ifelse(!is.null(tapestry_std), tapestry_std$m, NA),
      ifelse(!is.null(hyper_std), nrow(hyper_std), NA),
      ifelse(!is.null(hyperec_std), hyperec_std$m, NA)
    ),
    m_over_N = c(
      ifelse(!is.null(pbest_std), round(nrow(pbest_std)/384, 3), NA),
      ifelse(!is.null(tapestry_std), round(tapestry_std$m/256, 3), NA),
      ifelse(!is.null(hyper_std), round(nrow(hyper_std)/256, 3), NA),
      ifelse(!is.null(hyperec_std), round(hyperec_std$m/256, 3), NA)
    ),
    Lambda_max = c(
      ifelse(!is.null(pbest_std), max(Matrix::colSums(pbest_std)), NA),
      ifelse(!is.null(tapestry_std), max(Matrix::colSums(tapestry_std$matrix)), NA),
      ifelse(!is.null(hyper_std), max(Matrix::colSums(hyper_std)), NA),
      ifelse(!is.null(hyperec_std), hyperec_std$lambda_max, NA)
    ),
    Configuration = c("Standard (paper)", "Standard", "Standard", "Standard"),
    stringsAsFactors = FALSE
  )
  
  # =============================
  # Comparison 2: All at N=256
  # =============================
  message("2. Generating all designs at N=256 (common comparison)...")
  
  pbest_256 <- tryCatch(pp_matrix(q = 8, d = 3, nl = 6, N = 256), error = function(e) NULL)
  tapestry_256 <- tryCatch(tapestry_matrix(n = 256, k = 2), error = function(e) NULL)
  hyper_256 <- tryCatch(hyper_matrix(n = 256, m = 44, q = 2), error = function(e) NULL)
  hyperec_256 <- tryCatch(hyper_ec_matrix(n = 256, q = 2, alpha = 0.20), error = function(e) NULL)
  
  n256_df <- data.frame(
    Design = c("P-BEST", "Tapestry", "HYPER", "HYPER-EC"),
    N = rep(256, 4),
    M = c(
      ifelse(!is.null(pbest_256), nrow(pbest_256), NA),
      ifelse(!is.null(tapestry_256), tapestry_256$m, NA),
      ifelse(!is.null(hyper_256), nrow(hyper_256), NA),
      ifelse(!is.null(hyperec_256), hyperec_256$m, NA)
    ),
    m_over_N = c(
      ifelse(!is.null(pbest_256), round(nrow(pbest_256)/256, 3), NA),
      ifelse(!is.null(tapestry_256), round(tapestry_256$m/256, 3), NA),
      ifelse(!is.null(hyper_256), round(nrow(hyper_256)/256, 3), NA),
      ifelse(!is.null(hyperec_256), round(hyperec_256$m/256, 3), NA)
    ),
    Lambda_max = c(
      ifelse(!is.null(pbest_256), max(Matrix::colSums(pbest_256)), NA),
      ifelse(!is.null(tapestry_256), max(Matrix::colSums(tapestry_256$matrix)), NA),
      ifelse(!is.null(hyper_256), max(Matrix::colSums(hyper_256)), NA),
      ifelse(!is.null(hyperec_256), hyperec_256$lambda_max, NA)
    ),
    Configuration = c("Sub-optimal", "Optimal", "Optimal", "Optimal"),
    stringsAsFactors = FALSE
  )
  n256_df$Rank_Efficiency <- rank(
    ifelse(is.na(n256_df$m_over_N), Inf, n256_df$m_over_N),
    ties.method = "min"
  )
  
  # =============================
  # Comparison 3: All at N=384
  # =============================
  message("3. Generating all designs at N=384 (P-BEST optimal)...")
  
  pbest_384 <- tryCatch(pp_matrix(q = 8, d = 3, nl = 6, N = 384), error = function(e) NULL)
  tapestry_384 <- tryCatch(tapestry_matrix(n = 384, k = 4), error = function(e) NULL)
  hyper_384 <- tryCatch(hyper_matrix(n = 384, m = 72, q = 2), error = function(e) NULL)
  hyperec_384 <- tryCatch(hyper_ec_matrix(n = 384, m_base = 72, q = 2, alpha = 0.20), error = function(e) NULL)
  
  n384_df <- data.frame(
    Design = c("P-BEST", "Tapestry", "HYPER", "HYPER-EC"),
    N = rep(384, 4),
    M = c(
      ifelse(!is.null(pbest_384), nrow(pbest_384), NA),
      ifelse(!is.null(tapestry_384), tapestry_384$m, NA),
      ifelse(!is.null(hyper_384), nrow(hyper_384), NA),
      ifelse(!is.null(hyperec_384), hyperec_384$m, NA)
    ),
    m_over_N = c(
      ifelse(!is.null(pbest_384), round(nrow(pbest_384)/384, 3), NA),
      ifelse(!is.null(tapestry_384), round(tapestry_384$m/384, 3), NA),
      ifelse(!is.null(hyper_384), round(nrow(hyper_384)/384, 3), NA),
      ifelse(!is.null(hyperec_384), round(hyperec_384$m/384, 3), NA)
    ),
    Lambda_max = c(
      ifelse(!is.null(pbest_384), max(Matrix::colSums(pbest_384)), NA),
      ifelse(!is.null(tapestry_384), max(Matrix::colSums(tapestry_384$matrix)), NA),
      ifelse(!is.null(hyper_384), max(Matrix::colSums(hyper_384)), NA),
      ifelse(!is.null(hyperec_384), hyperec_384$lambda_max, NA)
    ),
    Configuration = c("Optimal", "Scaled", "Scaled", "Scaled"),
    stringsAsFactors = FALSE
  )
  n384_df$Rank_Efficiency <- rank(
    ifelse(is.na(n384_df$m_over_N), Inf, n384_df$m_over_N),
    ties.method = "min"
  )
  
  message("\n========================================")
  message("COMPARISON COMPLETE")
  message("========================================\n")

  .print_df_as_message <- function(df) {
    message(paste(utils::capture.output(print(df)), collapse = "\n"))
  }

  # Best (lowest m_over_N) design(s) at each common N, computed from the
  # actual tables rather than hard-coded, so the printed narrative and the
  # returned $summary always agree with $n256 / $n384.
  .best_at <- function(df) {
    best_ratio <- min(df$m_over_N, na.rm = TRUE)
    winners <- df$Design[!is.na(df$m_over_N) & df$m_over_N == best_ratio]
    list(
      label = paste0(paste(winners, collapse = " & "),
                      ifelse(length(winners) > 1L, " (tied)", "")),
      ratio = best_ratio
    )
  }
  best_256 <- .best_at(n256_df)
  best_384 <- .best_at(n384_df)

  lowest_lambda_idx <- which.min(n256_df$Lambda_max)
  lowest_lambda_label <- sprintf(
    "%s (lambda=%d)",
    n256_df$Design[lowest_lambda_idx],
    n256_df$Lambda_max[lowest_lambda_idx]
  )

  message("\nTABLE 1: STANDARD CONFIGURATIONS")
  message("(Each design at its package reference N)")
  message("----------------------------------------")
  .print_df_as_message(standard_df)

  message("\n\nTABLE 2: ALL DESIGNS AT N=256")
  message("(Fair to HYPER, Tapestry, HYPER-EC)")
  message("----------------------------------------")
  .print_df_as_message(n256_df)
  message(sprintf(
    "\nPACKAGE-SPECIFIC SUMMARY (efficiency): %s (m/N = %.3f)",
    best_256$label, best_256$ratio
  ))

  message("\n\nTABLE 3: ALL DESIGNS AT N=384")
  message("(Fair to P-BEST)")
  message("----------------------------------------")
  .print_df_as_message(n384_df)
  message(sprintf(
    "\nPACKAGE-SPECIFIC SUMMARY (efficiency): %s (m/N = %.3f)",
    best_384$label, best_384$ratio
  ))

  message("\n\nPACKAGE-SPECIFIC FINDINGS:")
  message("============================================")
  message(sprintf(
    "1. Package-specific summary: %s most efficient at N=384 (%.3f)",
    best_384$label, best_384$ratio
  ))
  message(sprintf(
    "2. Package-specific summary: %s at N=256 (%.3f)",
    best_256$label, best_256$ratio
  ))
  message("3. Package-specific summary: HYPER-EC trades efficiency for error correction")
  message("4. Package-specific summary: different designs excel at different N\n")

  list(
    standard = standard_df,
    n256 = n256_df,
    n384 = n384_df,
    summary = list(
      best_at_n256 = best_256$label,
      best_at_n384 = best_384$label,
      best_error_correction = "P-BEST (Strong RS)",
      lowest_lambda = lowest_lambda_label,
      most_flexible = "HYPER, Tapestry, HYPER-EC"
    )
  )
}
