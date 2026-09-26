test_that("compare_all_designs returns exact metrics for successful designs", {
  n <- 128L
  prevalence <- 0.01
  alpha <- 0.20
  k <- max(1L, as.integer(n * prevalence))
  m_hyper <- max(
    as.integer(ceiling(4 * k * log(n))),
    as.integer(ceiling(n * 0.15))
  )
  if (m_hyper %% 2L != 0L) {
    m_hyper <- m_hyper + 1L
  }

  set.seed(42)
  expect_message(
    comp <- compare_all_designs(n = n, prevalence = prevalence, alpha = alpha),
    "Comparison complete!"
  )

  expect_named(comp, c("pbest", "tapestry", "hyper", "hyper_ec", "comparison_table"))
  expect_equal(
    comp$comparison_table$Design,
    c("P-BEST", "Tapestry", "HYPER", "HYPER-EC")
  )
  expect_equal(
    names(comp$comparison_table),
    c("Design", "Samples_n", "Pools_m", "Efficiency_m_n", "Lambda_max", "Error_Correction")
  )
  expect_equal(dim(comp$comparison_table), c(4, 6))

  set.seed(42)
  pbest <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  tapestry <- tapestry_matrix(n = n, k = k, method = "random")
  hyper <- hyper_matrix(n = n, m = m_hyper, q = 2)
  hyper_ec <- hyper_ec_matrix(n = n, k = k, m_base = m_hyper, q = 2, alpha = alpha)

  expect_equal(comp$pbest$metrics$design_type, "P-BEST")
  expect_equal(comp$pbest$metrics$n, ncol(pbest))
  expect_equal(comp$pbest$metrics$m, nrow(pbest))
  expect_equal(comp$pbest$metrics$efficiency, nrow(pbest) / ncol(pbest))
  expect_equal(comp$pbest$metrics$lambda_max, max(Matrix::colSums(pbest)))
  expect_equal(comp$pbest$metrics$avg_lambda, mean(Matrix::colSums(pbest)))
  expect_equal(comp$pbest$metrics$pool_size_mean, mean(Matrix::rowSums(pbest)))
  expect_equal(comp$pbest$metrics$pool_size_sd, sd(Matrix::rowSums(pbest)))
  expect_equal(comp$pbest$metrics$error_correction, "Strong (Reed-Solomon)")

  expect_equal(comp$tapestry$metrics$design_type, "Tapestry")
  expect_equal(comp$tapestry$metrics$n, tapestry$n)
  expect_equal(comp$tapestry$metrics$m, tapestry$m)
  expect_equal(comp$tapestry$metrics$efficiency, tapestry$m / tapestry$n)
  expect_equal(comp$tapestry$metrics$lambda_max, tapestry$lambda_max)
  expect_equal(comp$tapestry$metrics$avg_lambda, mean(Matrix::colSums(tapestry$matrix)))
  expect_equal(comp$tapestry$metrics$pool_size_mean, mean(Matrix::rowSums(tapestry$matrix)))
  expect_equal(comp$tapestry$metrics$pool_size_sd, sd(Matrix::rowSums(tapestry$matrix)))
  expect_equal(comp$tapestry$metrics$error_correction, "Implicit (Compressed Sensing)")

  expect_equal(comp$hyper$metrics$design_type, "HYPER")
  expect_equal(comp$hyper$metrics$n, ncol(hyper))
  expect_equal(comp$hyper$metrics$m, nrow(hyper))
  expect_equal(comp$hyper$metrics$efficiency, nrow(hyper) / ncol(hyper))
  expect_equal(comp$hyper$metrics$lambda_max, max(Matrix::colSums(hyper)))
  expect_equal(comp$hyper$metrics$avg_lambda, mean(Matrix::colSums(hyper)))
  expect_equal(comp$hyper$metrics$pool_size_mean, mean(Matrix::rowSums(hyper)))
  expect_equal(comp$hyper$metrics$pool_size_sd, sd(Matrix::rowSums(hyper)))
  expect_equal(comp$hyper$metrics$error_correction, "None")

  expect_equal(comp$hyper_ec$metrics$design_type, "HYPER-EC")
  expect_equal(comp$hyper_ec$metrics$n, hyper_ec$n)
  expect_equal(comp$hyper_ec$metrics$m, hyper_ec$m)
  expect_equal(comp$hyper_ec$metrics$m_base, hyper_ec$m_base)
  expect_equal(comp$hyper_ec$metrics$k_parity, hyper_ec$k_parity)
  expect_equal(comp$hyper_ec$metrics$efficiency, hyper_ec$efficiency)
  expect_equal(comp$hyper_ec$metrics$lambda_max, hyper_ec$lambda_max)
  expect_equal(comp$hyper_ec$metrics$avg_lambda, mean(Matrix::colSums(hyper_ec$matrix)))
  expect_equal(comp$hyper_ec$metrics$pool_size_mean, mean(Matrix::rowSums(hyper_ec$matrix)))
  expect_equal(comp$hyper_ec$metrics$pool_size_sd, sd(Matrix::rowSums(hyper_ec$matrix)))
  expect_equal(
    comp$hyper_ec$metrics$error_correction,
    paste0("Medium (XOR, corrects ", hyper_ec$error_correction, " errors)")
  )

  expect_equal(comp$comparison_table$Samples_n, c(
    ncol(pbest), tapestry$n, ncol(hyper), hyper_ec$n
  ))
  expect_equal(comp$comparison_table$Pools_m, c(
    nrow(pbest), tapestry$m, nrow(hyper), hyper_ec$m
  ))
  expect_equal(comp$comparison_table$Efficiency_m_n, c(
    round(nrow(pbest) / ncol(pbest), 3),
    round(tapestry$m / tapestry$n, 3),
    round(nrow(hyper) / ncol(hyper), 3),
    round(hyper_ec$efficiency, 3)
  ))
  expect_equal(comp$comparison_table$Lambda_max, c(
    max(Matrix::colSums(pbest)),
    tapestry$lambda_max,
    max(Matrix::colSums(hyper)),
    hyper_ec$lambda_max
  ))
  expect_equal(comp$comparison_table$Error_Correction, c(
    "Strong (Reed-Solomon)",
    "Implicit (Compressed Sensing)",
    "None",
    paste0("Medium (XOR, corrects ", hyper_ec$error_correction, " errors)")
  ))
})

test_that("compare_all_designs records failed design generation exactly", {
  local_mocked_bindings(
    pp_matrix = function(...) stop("pp fail"),
    tapestry_matrix = function(...) stop("tap fail"),
    hyper_matrix = function(...) stop("hyper fail"),
    hyper_ec_matrix = function(...) stop("hyper-ec fail")
  )

  expect_warning(
    expect_warning(
      expect_warning(
        expect_warning(
          comp <- compare_all_designs(n = 128, prevalence = 0.01, alpha = 0.2),
          "P-BEST generation failed: pp fail"
        ),
        "Tapestry generation failed: tap fail"
      ),
      "HYPER generation failed: hyper fail"
    ),
    "HYPER-EC generation failed: hyper-ec fail"
  )

  expect_null(comp$pbest$matrix)
  expect_null(comp$tapestry$matrix)
  expect_null(comp$hyper$matrix)
  expect_null(comp$hyper_ec$matrix)

  expect_equal(comp$pbest$metrics, list(design_type = "P-BEST", note = "Generation failed"))
  expect_equal(comp$tapestry$metrics, list(design_type = "Tapestry", note = "Generation failed"))
  expect_equal(comp$hyper$metrics, list(design_type = "HYPER", note = "Generation failed"))
  expect_equal(comp$hyper_ec$metrics, list(design_type = "HYPER-EC", note = "Generation failed"))

  expect_true(all(is.na(comp$comparison_table$Samples_n)))
  expect_true(all(is.na(comp$comparison_table$Pools_m)))
  expect_true(all(is.na(comp$comparison_table$Efficiency_m_n)))
  expect_true(all(is.na(comp$comparison_table$Lambda_max)))
  expect_equal(comp$comparison_table$Error_Correction, rep("N/A", 4))
})

test_that(".find_nearest_prime_power finds the closest supported prime power", {
  expect_equal(.find_nearest_prime_power(2L), 2L)
  expect_equal(.find_nearest_prime_power(6L), 5L)
  expect_equal(.find_nearest_prime_power(15L), 16L)
  expect_equal(.find_nearest_prime_power(28L), 27L)
  expect_equal(.find_nearest_prime_power(50L), 49L)
})
