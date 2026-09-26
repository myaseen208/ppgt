test_that("ppgt_notations returns the expected registry structure", {
  notations <- ppgt_notations()

  expect_s3_class(notations, "data.frame")
  expect_named(
    notations,
    c(
      "symbol",
      "render",
      "meaning",
      "category",
      "status",
      "scope",
      "family",
      "alias_of",
      "source",
      "notes"
    )
  )
  expect_true(nrow(notations) >= 20L)
})


test_that("ppgt_notations contains canonical package-wide symbols", {
  notations <- ppgt_notations()

  expect_true(any(
    notations$symbol == "N" &
      notations$category == "global_canonical" &
      notations$status == "canonical"
  ))
  expect_true(any(
    notations$symbol == "J" &
      notations$category == "global_canonical" &
      notations$status == "canonical"
  ))
  expect_true(any(
    notations$symbol == "\\mathbf{M}" &
      notations$category == "global_canonical"
  ))
  expect_true(any(
    notations$symbol == "\\widetilde{\\mathbf{y}}" &
      notations$category == "global_canonical"
  ))
  expect_true(any(
    notations$symbol == "w_i" &
      notations$category == "global_canonical"
  ))
  expect_true(any(
    notations$symbol == "\\widetilde{\\mathbf{z}}" &
      notations$category == "global_canonical"
  ))
  expect_true(any(
    notations$symbol == "\\mathbf{z}" &
      notations$category == "global_canonical"
  ))
})


test_that("ppgt_notations records aliases and deprecated forms", {
  notations <- ppgt_notations()

  expect_true(any(
    notations$symbol == "pool_size" &
      notations$category == "display_alias" &
      notations$alias_of == "n_j"
  ))
  expect_true(any(
    notations$symbol == "pools_per_sample" &
      notations$category == "display_alias" &
      notations$alias_of == "w_i"
  ))
  expect_true(any(
    notations$symbol == "M" &
      notations$category == "display_alias" &
      notations$alias_of == "J"
  ))
  expect_true(any(
    notations$symbol == "n_i" &
      notations$category == "forbidden_alias" &
      notations$alias_of == "w_i"
  ))
  expect_true(any(
    notations$symbol == "\\widetilde{\\mathbf{Y}}" &
      notations$category == "forbidden_alias" &
      notations$alias_of == "\\widetilde{\\mathbf{y}}"
  ))
})
