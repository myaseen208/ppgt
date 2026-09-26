test_that("compare_all_designs_honest returns exact comparison tables", {
  # All of compare_all_designs_honest()'s output goes through message(), not
  # cat()/print() (project convention), so capture_messages() both captures
  # the narrative and lets `comp <-` propagate normally.
  msgs <- testthat::capture_messages(comp <- compare_all_designs_honest())
  expect_true(any(grepl("COMPREHENSIVE DESIGN COMPARISON", msgs, fixed = TRUE)))

  expect_named(comp, c("standard", "n256", "n384", "summary"))
  expect_named(
    comp$summary,
    c(
      "best_at_n256",
      "best_at_n384",
      "best_error_correction",
      "lowest_lambda",
      "most_flexible"
    )
  )

  # comp$summary and the printed narrative are derived from the actual
  # computed m_over_N / Lambda_max columns (not hard-coded), so these
  # expectations are cross-checked against comp$n256 / comp$n384 below
  # rather than pinned to literal strings that could silently drift out of
  # sync with the tables.
  expect_equal(comp$summary$best_at_n256, "HYPER")
  expect_equal(comp$summary$best_at_n384, "P-BEST")
  expect_equal(comp$summary$best_error_correction, "P-BEST (Strong RS)")
  expect_equal(comp$summary$lowest_lambda, "HYPER (lambda=2)")
  expect_equal(comp$summary$most_flexible, "HYPER, Tapestry, HYPER-EC")

  expect_true(any(grepl("TABLE 1: STANDARD CONFIGURATIONS", msgs, fixed = TRUE)))
  expect_true(any(grepl("TABLE 2: ALL DESIGNS AT N=256", msgs, fixed = TRUE)))
  expect_true(any(grepl("TABLE 3: ALL DESIGNS AT N=384", msgs, fixed = TRUE)))

  expected_designs <- c("P-BEST", "Tapestry", "HYPER", "HYPER-EC")
  expect_equal(comp$standard$Design, expected_designs)
  expect_equal(comp$n256$Design, expected_designs)
  expect_equal(comp$n384$Design, expected_designs)

  expect_equal(dim(comp$standard), c(4, 6))
  expect_equal(dim(comp$n256), c(4, 7))
  expect_equal(dim(comp$n384), c(4, 7))

  expect_equal(comp$standard$N_Standard, c(384, 256, 256, 256))
  expect_equal(comp$n256$N, rep(256, 4))
  expect_equal(comp$n384$N, rep(384, 4))
  # Ranks are computed from the actual m_over_N column (1 = most efficient,
  # i.e. lowest J/N), not hard-coded.
  expect_equal(
    comp$n256$Rank_Efficiency,
    rank(comp$n256$m_over_N, ties.method = "min")
  )
  expect_equal(
    comp$n384$Rank_Efficiency,
    rank(comp$n384$m_over_N, ties.method = "min")
  )

  pbest_std <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  tapestry_std <- tapestry_matrix(n = 256, k = 2)
  hyper_std <- hyper_matrix(n = 256, m = 44, q = 2)
  hyperec_std <- hyper_ec_matrix(n = 256, q = 2, alpha = 0.20)

  expect_equal(comp$standard$M, c(
    nrow(pbest_std),
    tapestry_std$m,
    nrow(hyper_std),
    hyperec_std$m
  ))
  expect_equal(comp$standard$m_over_N, c(
    round(nrow(pbest_std) / 384, 3),
    round(tapestry_std$m / 256, 3),
    round(nrow(hyper_std) / 256, 3),
    round(hyperec_std$m / 256, 3)
  ))
  expect_equal(comp$standard$Lambda_max, c(
    max(Matrix::colSums(pbest_std)),
    max(Matrix::colSums(tapestry_std$matrix)),
    max(Matrix::colSums(hyper_std)),
    hyperec_std$lambda_max
  ))
  expect_equal(
    comp$standard$Configuration,
    c("Standard (paper)", "Standard", "Standard", "Standard")
  )

  pbest_256 <- pp_matrix(q = 8, d = 3, nl = 6, N = 256)
  tapestry_256 <- tapestry_matrix(n = 256, k = 2)
  hyper_256 <- hyper_matrix(n = 256, m = 44, q = 2)
  hyperec_256 <- hyper_ec_matrix(n = 256, q = 2, alpha = 0.20)

  expect_equal(comp$n256$M, c(
    nrow(pbest_256),
    tapestry_256$m,
    nrow(hyper_256),
    hyperec_256$m
  ))
  expect_equal(comp$n256$m_over_N, c(
    round(nrow(pbest_256) / 256, 3),
    round(tapestry_256$m / 256, 3),
    round(nrow(hyper_256) / 256, 3),
    round(hyperec_256$m / 256, 3)
  ))
  expect_equal(comp$n256$Lambda_max, c(
    max(Matrix::colSums(pbest_256)),
    max(Matrix::colSums(tapestry_256$matrix)),
    max(Matrix::colSums(hyper_256)),
    hyperec_256$lambda_max
  ))
  expect_equal(
    comp$n256$Configuration,
    c("Sub-optimal", "Optimal", "Optimal", "Optimal")
  )

  pbest_384 <- pp_matrix(q = 8, d = 3, nl = 6, N = 384)
  tapestry_384 <- tapestry_matrix(n = 384, k = 4)
  hyper_384 <- hyper_matrix(n = 384, m = 72, q = 2)
  hyperec_384 <- hyper_ec_matrix(n = 384, m_base = 72, q = 2, alpha = 0.20)

  expect_equal(comp$n384$M, c(
    nrow(pbest_384),
    tapestry_384$m,
    nrow(hyper_384),
    hyperec_384$m
  ))
  expect_equal(comp$n384$m_over_N, c(
    round(nrow(pbest_384) / 384, 3),
    round(tapestry_384$m / 384, 3),
    round(nrow(hyper_384) / 384, 3),
    round(hyperec_384$m / 384, 3)
  ))
  expect_equal(comp$n384$Lambda_max, c(
    max(Matrix::colSums(pbest_384)),
    max(Matrix::colSums(tapestry_384$matrix)),
    max(Matrix::colSums(hyper_384)),
    hyperec_384$lambda_max
  ))
  expect_equal(
    comp$n384$Configuration,
    c("Optimal", "Scaled", "Scaled", "Scaled")
  )
})
