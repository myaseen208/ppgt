testthat::test_that("PoolMatrix dispatches to pp_matrix", {
  M_expected <- ppgt::pp_matrix(q = 4L, d = 3L, nl = 5L, N = 64L)
  M_actual <- ppgt::PoolMatrix(family = "pp", q = 4L, d = 3L, nl = 5L, N = 64L)

  testthat::expect_s4_class(M_actual, "dgCMatrix")
  testthat::expect_equal(as.matrix(M_actual), as.matrix(M_expected))
  testthat::expect_equal(attr(M_actual, "poolmatrix")$family, "pp")
  testthat::expect_false("stage1" %in% names(attributes(M_actual)))
  testthat::expect_false("stage2" %in% names(attributes(M_actual)))
  testthat::expect_false("final_calls" %in% names(attributes(M_actual)))
})

testthat::test_that("PoolMatrix can infer PP nl from J and q", {
  M_actual <- ppgt::PoolMatrix(family = "pp", q = 4L, d = 3L, J = 20L, N = 64L)

  testthat::expect_equal(dim(M_actual), c(20L, 64L))
  testthat::expect_equal(attr(M_actual, "pp")$nl, 5L)
})

testthat::test_that("PoolMatrix dispatches to PP-based and clinical PBEST paths", {
  M_pp <- ppgt::PoolMatrix(family = "pbest")
  M_expected <- ppgt::pp_matrix(q = 8L, d = 3L, nl = 6L, N = 384L)

  testthat::expect_equal(as.matrix(M_pp), as.matrix(M_expected))
  testthat::expect_equal(attr(M_pp, "poolmatrix")$method, "pp")

  M_clinical <- ppgt::PoolMatrix(
    family = "pbest",
    method = "clinical",
    N = 384L,
    J = 94L,
    seed = 42L
  )

  testthat::expect_equal(dim(M_clinical), c(94L, 384L))
  testthat::expect_equal(attr(M_clinical, "poolmatrix")$method, "clinical")
})

testthat::test_that("PoolMatrix dispatches to hyper_matrix", {
  M_expected <- ppgt::hyper_matrix(n = 12L, m = 6L, q = 2L, reorder = FALSE)
  M_actual <- ppgt::PoolMatrix(
    family = "hyper",
    N = 12L,
    J = 6L,
    q = 2L,
    reorder = FALSE
  )

  testthat::expect_equal(as.matrix(M_actual), as.matrix(M_expected))
  testthat::expect_equal(attr(M_actual, "poolmatrix")$family, "hyper")
})

testthat::test_that("PoolMatrix dispatches to separable and disjunct constructors", {
  M_sep <- ppgt::PoolMatrix(
    family = "separable",
    N = 16L,
    d = 2L,
    method = "reed-solomon",
    verify = TRUE
  )
  M_disj <- ppgt::PoolMatrix(
    family = "disjunct",
    N = 16L,
    d = 2L,
    method = "reed-solomon",
    verify = TRUE
  )

  testthat::expect_s4_class(M_sep, "dgCMatrix")
  testthat::expect_s4_class(M_disj, "dgCMatrix")
  testthat::expect_true(isTRUE(as.logical(attr(M_sep, "separable_verified"))))
  testthat::expect_true(isTRUE(as.logical(attr(M_disj, "disjunct_verified"))))
})

testthat::test_that("PoolMatrix dispatches to additional deterministic families", {
  M_dorfman <- ppgt::PoolMatrix(family = "dorfman", N = 12L, g = 3L)
  M_array <- ppgt::PoolMatrix(family = "array", method = "2d", r = 4L, c = 5L)
  M_array3 <- ppgt::PoolMatrix(
    family = "array",
    method = "3d",
    d1 = 2L,
    d2 = 2L,
    d3 = 2L
  )
  M_bibd <- ppgt::PoolMatrix(family = "bibd", N = 7L, k = 3L)
  M_cube <- ppgt::PoolMatrix(family = "hypercube", q = 2L, d = 3L)
  M_kirkman <- ppgt::PoolMatrix(family = "kirkman", N = 9L)
  M_pg <- ppgt::PoolMatrix(family = "pg", q = 2L)

  testthat::expect_equal(dim(M_dorfman), c(4L, 12L))
  testthat::expect_equal(dim(M_array), c(9L, 20L))
  testthat::expect_equal(dim(M_array3), c(6L, 8L))
  testthat::expect_equal(dim(M_bibd), c(7L, 7L))
  testthat::expect_equal(dim(M_cube), c(6L, 8L))
  testthat::expect_equal(dim(M_kirkman), c(12L, 9L))
  testthat::expect_equal(dim(M_pg), c(7L, 7L))
})

testthat::test_that("PoolMatrix errors on missing family-specific arguments", {
  testthat::expect_error(
    ppgt::PoolMatrix(family = "pp", q = 4L, d = 3L),
    "requires `nl`, or `J`"
  )
  testthat::expect_error(
    ppgt::PoolMatrix(family = "hyper", N = 12L, q = 2L),
    "requires `J`"
  )
  testthat::expect_error(
    ppgt::PoolMatrix(family = "dorfman", N = 12L),
    "requires `g`"
  )
})

testthat::test_that("PoolMatrix validates implied dimensions for fixed-size families", {
  testthat::expect_error(
    ppgt::PoolMatrix(family = "hypercube", q = 2L, d = 3L, N = 9L),
    "N must equal"
  )
  testthat::expect_error(
    ppgt::PoolMatrix(family = "pg", q = 2L, J = 8L),
    "J must equal"
  )
})
