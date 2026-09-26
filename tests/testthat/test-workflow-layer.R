testthat::test_that("simulate_group_testing returns latent and observed outcomes", {
  M <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L)
  Y_tilde <- integer(64L)
  Y_tilde[c(2L, 13L)] <- 1L

  sim <- ppgt::simulate_group_testing(M, Y_tilde)

  testthat::expect_true(is.list(sim))
  testthat::expect_equal(sim$matrix, M)
  testthat::expect_equal(sim$Y_tilde, Y_tilde)
  testthat::expect_equal(sim$z, sim$z_tilde)
  testthat::expect_equal(length(sim$z), nrow(M))
  testthat::expect_identical(sim$mode, "noiseless")
})

testthat::test_that("simulate_group_testing validates dimensions and probabilities", {
  M <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L)

  testthat::expect_error(
    ppgt::simulate_group_testing(M, integer(10L)),
    "Length of `Y_tilde`"
  )
  testthat::expect_error(
    ppgt::simulate_group_testing(M, integer(64L), s_e = 1.2),
    "`s_e` must lie in \\[0, 1\\]"
  )
})

testthat::test_that("simulate_group_testing supports noisy pool outcomes with a reproducible seed", {
  M <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L)
  Y_tilde <- integer(64L)
  Y_tilde[c(2L, 13L)] <- 1L

  sim1 <- ppgt::simulate_group_testing(M, Y_tilde, s_e = 0.95, s_p = 0.99, seed = 42L)
  sim2 <- ppgt::simulate_group_testing(M, Y_tilde, s_e = 0.95, s_p = 0.99, seed = 42L)

  testthat::expect_identical(sim1$mode, "noisy")
  testthat::expect_equal(sim1$z, sim2$z)  # same seed -> reproducible noisy outcomes
  testthat::expect_length(sim1$s_e, nrow(M))
  testthat::expect_length(sim1$s_p, nrow(M))
  testthat::expect_true(all(sim1$s_e == 0.95))
  testthat::expect_true(all(sim1$s_p == 0.99))

  # A .Random.seed already present in .GlobalEnv is saved and restored, not
  # merely removed.
  set.seed(1L)
  before <- get(".Random.seed", envir = .GlobalEnv)
  invisible(ppgt::simulate_group_testing(M, Y_tilde, s_e = 0.95, seed = 7L))
  after <- get(".Random.seed", envir = .GlobalEnv)
  testthat::expect_identical(before, after)
})

testthat::test_that("run_testing_workflow supports non-adaptive pp_decode stage 1", {
  M <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L)
  Y_tilde <- integer(64L)
  Y_tilde[c(2L, 13L)] <- 1L
  sim <- ppgt::simulate_group_testing(M, Y_tilde)

  wf <- ppgt::run_testing_workflow(
    design = M,
    z = sim$z,
    workflow = "non_adaptive",
    decoder = "pp_decode",
    decoder_args = list(verbose = FALSE),
    stage2 = "none"
  )

  testthat::expect_s3_class(wf, "ppgt_workflow")
  testthat::expect_equal(wf$workflow, "non_adaptive")
  testthat::expect_true(is.list(wf$stage1))
  testthat::expect_true(is.list(wf$final_calls))
  testthat::expect_true(all(c(2L, 13L) %in% wf$final_calls$positive))
})

testthat::test_that("run_testing_workflow supports adaptive individual retest consolidation", {
  M <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L)
  Y_tilde <- integer(64L)
  Y_tilde[c(2L, 13L)] <- 1L
  sim <- ppgt::simulate_group_testing(M, Y_tilde)

  individual_results <- integer(64L)
  individual_results[c(2L, 13L)] <- 1L

  wf <- ppgt::run_testing_workflow(
    design = M,
    z = sim$z,
    workflow = "adaptive_retest",
    decoder = "pp_decode",
    decoder_args = list(verbose = FALSE),
    stage2 = "individual_retest",
    individual_results = individual_results
  )

  testthat::expect_equal(wf$workflow, "adaptive_retest")
  testthat::expect_true(length(wf$stage2$planned) >= 2L)
  testthat::expect_true(all(c(2L, 13L) %in% wf$final_calls$positive))
  testthat::expect_true(length(wf$final_calls$negative) >= 0L)
  testthat::expect_equal(
    sort(unique(c(
      wf$final_calls$positive,
      wf$final_calls$negative,
      wf$final_calls$unresolved
    ))),
    seq_len(ncol(M))
  )
})

testthat::test_that("run_testing_workflow validates workflow-stage2 combinations", {
  M <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L)
  z <- integer(nrow(M))

  testthat::expect_error(
    ppgt::run_testing_workflow(
      design = M,
      z = z,
      workflow = "non_adaptive",
      decoder = "pp_decode",
      stage2 = "individual_retest"
    ),
    "requires `stage2 = \"none\"`"
  )
  testthat::expect_error(
    ppgt::run_testing_workflow(
      design = M,
      z = z,
      workflow = "adaptive_retest",
      decoder = "pp_decode",
      stage2 = "none"
    ),
    "requires `stage2 = \"individual_retest\"`"
  )
  testthat::expect_error(
    ppgt::run_testing_workflow(
      design = M,
      z = z,
      workflow = "non_adaptive",
      decoder = "pp_decode",
      stage2 = "none",
      individual_results = integer(ncol(M))
    ),
    "`individual_results` may be supplied only"
  )
})

testthat::test_that("protocol_summary returns compact workflow counts", {
  M <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L)
  Y_tilde <- integer(64L)
  Y_tilde[2L] <- 1L
  sim <- ppgt::simulate_group_testing(M, Y_tilde)

  wf <- ppgt::run_testing_workflow(
    design = M,
    z = sim$z,
    workflow = "non_adaptive",
    decoder = "pp_decode",
    decoder_args = list(verbose = FALSE),
    stage2 = "none"
  )

  summary_obj <- ppgt::protocol_summary(wf)

  testthat::expect_true(is.list(summary_obj))
  testthat::expect_equal(summary_obj$n_pools, nrow(M))
  testthat::expect_equal(summary_obj$n_individuals, ncol(M))
  testthat::expect_named(
    summary_obj,
    c(
      "workflow",
      "n_pools",
      "n_individuals",
      "n_stage1_positive",
      "n_stage1_suspected",
      "n_stage1_negative",
      "n_stage2_planned",
      "n_final_positive",
      "n_final_negative",
      "n_final_unresolved"
    )
  )
})
