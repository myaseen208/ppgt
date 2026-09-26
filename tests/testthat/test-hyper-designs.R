test_that("hyper_matrix uses upstream reordered layout by default", {
  mat_default <- hyper_matrix(n = 20L, m = 6L, q = 2L)
  mat_explicit <- hyper_matrix(n = 20L, m = 6L, q = 2L, reorder = TRUE)
  mat_unordered <- hyper_matrix(n = 20L, m = 6L, q = 2L, reorder = FALSE)

  expect_identical(as.matrix(mat_default), as.matrix(mat_explicit))
  expect_false(identical(as.matrix(mat_default), as.matrix(mat_unordered)))
  expect_true(all(Matrix::colSums(mat_default) == 2))
})

test_that("hyper_matrix matches the exact small q = 2 canonical layout", {
  mat <- hyper_matrix(n = 12L, m = 6L, q = 2L, reorder = FALSE)

  expected <- matrix(
    c(
      1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L,
      0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L,
      0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 0L,
      0L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L,
      0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 0L,
      1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L
    ),
    nrow = 6L,
    byrow = TRUE
  )

  expect_equal(as.matrix(mat), expected)
  expect_true(all(Matrix::rowSums(mat) == 4L))
  expect_true(all(Matrix::colSums(mat) == 2L))
})

test_that("hyper_matrix cycles canonical factors when n exceeds factor count", {
  mat <- hyper_matrix(n = 8L, m = 4L, q = 1L, reorder = FALSE)

  expect_equal(dim(mat), c(4L, 8L))
  expect_true(all(Matrix::colSums(mat) == 1))
  expect_equal(as.integer(Matrix::rowSums(mat)), c(2L, 2L, 2L, 2L))
})

test_that("hyper_matrix exact q = 1 cycling layout is stable", {
  mat <- hyper_matrix(n = 5L, m = 4L, q = 1L, reorder = FALSE)

  expected <- matrix(
    c(
      1L, 0L, 0L, 0L, 1L,
      0L, 1L, 0L, 0L, 0L,
      0L, 0L, 1L, 0L, 0L,
      0L, 0L, 0L, 1L, 0L
    ),
    nrow = 4L,
    byrow = TRUE
  )

  expect_equal(as.matrix(mat), expected)
})

test_that("hyper_matrix supports q = 3 reference constraints", {
  mat <- hyper_matrix(n = 12L, m = 12L, q = 3L)

  expect_equal(dim(mat), c(12L, 12L))
  expect_true(all(Matrix::colSums(mat) == 3))
  expect_true(all(Matrix::rowSums(mat) >= 1))
})

test_that("hyper_matrix exact q = 3 small layout has balanced structure", {
  mat <- hyper_matrix(n = 12L, m = 12L, q = 3L, reorder = FALSE)

  expect_equal(dim(mat), c(12L, 12L))
  expect_true(all(Matrix::colSums(mat) == 3L))
  expect_true(all(Matrix::rowSums(mat) == 3L))
  expect_equal(
    as.matrix(mat[, 1:6]),
    matrix(
      c(
        1L, 0L, 0L, 0L, 0L, 1L,
        0L, 0L, 0L, 1L, 1L, 0L,
        0L, 1L, 0L, 0L, 0L, 0L,
        0L, 0L, 1L, 0L, 0L, 1L,
        0L, 0L, 1L, 0L, 0L, 0L,
        0L, 0L, 0L, 1L, 0L, 0L,
        0L, 0L, 0L, 1L, 0L, 0L,
        0L, 1L, 0L, 0L, 0L, 0L,
        0L, 0L, 1L, 0L, 0L, 1L,
        1L, 0L, 0L, 0L, 0L, 0L,
        0L, 1L, 0L, 0L, 1L, 0L,
        1L, 0L, 0L, 0L, 1L, 0L
      ),
      nrow = 12L,
      byrow = TRUE
    )
  )
})

test_that("hyper_matrix validates q-specific arithmetic constraints", {
  expect_error(hyper_matrix(n = 10L, m = 5L, q = 2L), "m must be even")
  expect_error(hyper_matrix(n = 12L, m = 10L, q = 3L), "divisible by 6")
  expect_error(hyper_matrix(n = 12L, m = 36L, q = 3L), "m-1 must be prime")
  expect_error(hyper_matrix(n = 12L, m = 6L, q = 4L), "q must be 1, 2, or 3")
})

test_that("HyperDesign is an exact compatibility wrapper around hyper_matrix", {
  wrapped <- HyperDesign(n = 12L, m = 6L, q = 2L, reorder = FALSE)
  direct <- hyper_matrix(n = 12L, m = 6L, q = 2L, reorder = FALSE)

  expect_s4_class(wrapped, "dgCMatrix")
  expect_identical(as.matrix(wrapped), as.matrix(direct))
  expect_identical(attr(wrapped, "design"), attr(direct, "design"))
})

test_that("simple syndrome decoder recovers multi-error patterns", {
  parity_matrix <- Matrix::Matrix(
    matrix(
      c(
        1L, 0L, 0L, 0L,
        0L, 1L, 0L, 0L,
        0L, 0L, 1L, 0L,
        0L, 0L, 0L, 1L
      ),
      nrow = 4L,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  syndrome <- c(1L, 0L, 1L, 0L)

  error_pattern <- .simple_syndrome_decode(
    syndrome = syndrome,
    parity_matrix = parity_matrix,
    m_base = 4L,
    max_errors = 2L
  )

  expect_equal(error_pattern, c(1L, 0L, 1L, 0L))
})

test_that("hyper_ec_decode applies syndrome correction before decoding", {
  matrix_base <- Matrix::Matrix(
    matrix(
      c(
        1L, 1L, 0L,
        1L, 0L, 1L,
        0L, 1L, 1L
      ),
      nrow = 3L,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  matrix_parity <- Matrix::Diagonal(3L)
  matrix_full <- rbind(matrix_base, matrix_base)
  design <- list(
    matrix = matrix_full,
    matrix_base = matrix_base,
    matrix_parity = matrix_parity,
    m = 6L,
    m_base = 3L,
    k_parity = 3L,
    n = 3L,
    error_correction = 1L
  )
  x_true <- c(1L, 0L, 0L)
  y_true <- as.integer((matrix_full %*% x_true) > 0)
  y_obs <- y_true
  y_obs[2] <- 1L - y_obs[2]

  result <- hyper_ec_decode(y_obs, design, method = "syndrome")

  expect_equal(result$errors_corrected, 1)
  expect_equal(result$y_corrected, y_true[1:3])
  expect_equal(result$x_decoded, x_true)
})
