read_utf8_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}


escape_regex <- function(x) {
  gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", x)
}


collect_matches <- function(paths, pattern) {
  hits <- list()

  for (path in paths) {
    if (!file.exists(path)) {
      next
    }

    txt <- read_utf8_text(path)
    if (grepl(pattern, txt, perl = TRUE)) {
      hits[[path]] <- regmatches(
        txt,
        gregexpr(pattern, txt, perl = TRUE)
      )[[1L]]
    }
  }

  hits
}


pkg_root <- normalizePath(testthat::test_path("..", ".."))

source_doc_files <- c(
  list.files(file.path(pkg_root, "R"), pattern = "[.]R$", full.names = TRUE),
  list.files(
    file.path(pkg_root, "vignettes"),
    pattern = "[.]Rmd$",
    full.names = TRUE
  ),
  list.files(
    file.path(pkg_root, "vignettes-src"),
    pattern = "[.]qmd$",
    full.names = TRUE
  ),
  file.path(pkg_root, "README.md")
)

man_files <- list.files(
  file.path(pkg_root, "man"),
  pattern = "[.]Rd$",
  full.names = TRUE
)


test_that("forbidden notation aliases do not appear outside the registry", {
  source_exclusions <- normalizePath(
    c(
      file.path(pkg_root, "R", "notations.R")
    ),
    mustWork = FALSE
  )

  source_targets <- setdiff(normalizePath(source_doc_files, mustWork = FALSE), source_exclusions)

  forbidden_patterns <- c(
    "\\\\widetilde\\{\\\\mathbf\\{Y\\}\\}",
    "\\bn_i\\b"
  )

  for (pattern in forbidden_patterns) {
    hits <- collect_matches(source_targets, pattern)
    if (length(hits) != 0L) {
      fail(paste(
        "Unexpected forbidden notation pattern in source docs:",
        pattern,
        paste(names(hits), collapse = ", ")
      ))
    }
  }

  succeed()
})


test_that("scalar M is not used as an unlabeled canonical pool-count symbol", {
  source_exclusions <- normalizePath(
    c(
      file.path(pkg_root, "R", "notations.R")
    ),
    mustWork = FALSE
  )

  source_targets <- setdiff(normalizePath(source_doc_files, mustWork = FALSE), source_exclusions)

  bad_hits <- character(0L)

  for (path in source_targets) {
    if (!file.exists(path)) {
      next
    }

    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    flagged <- grepl(
      "\\\\eqn\\{M\\}.*number of pools|number of pools.*\\\\eqn\\{M\\}|\\bM\\s*=\\s*number of pools\\b|\\bnumber of pools\\s*=\\s*M\\b",
      lines,
      perl = TRUE
    )
    flagged <- flagged & !grepl("legacy display alias|legacy display label", lines, perl = TRUE)

    if (any(flagged)) {
      bad_hits <- c(
        bad_hits,
        paste0(path, ":", which(flagged))
      )
    }
  }

  if (length(bad_hits) != 0L) {
    fail(paste("Unlabeled scalar M usage:", paste(bad_hits, collapse = ", ")))
  }

  succeed()
})


test_that("decoder-local notation stays in decoder-local files", {
  decoder_patterns <- c(
    "\\\\mathbf\\{A\\}",
    "\\\\mathbf\\{x\\}",
    "(?<!widetilde\\{)\\\\mathbf\\{y\\}",
    "\\\\hat\\{\\\\mathbf\\{x\\}\\}",
    "\\\\mathbf\\{g\\}_k",
    "\\\\tau\\b"
  )

  allowed_source_files <- normalizePath(
    c(
      file.path(pkg_root, "R", "notations.R"),
      file.path(pkg_root, "R", "pp_decode.R"),
      file.path(pkg_root, "R", "tapestry_hyperec.R"),
      file.path(pkg_root, "vignettes", "ppgt_pp.Rmd"),
      file.path(pkg_root, "vignettes-src", "ppgt_pp.qmd")
    ),
    mustWork = FALSE
  )

  source_targets <- setdiff(normalizePath(source_doc_files, mustWork = FALSE), allowed_source_files)

  for (pattern in decoder_patterns) {
    hits <- collect_matches(source_targets, pattern)
    if (length(hits) != 0L) {
      fail(paste(
        "Unexpected decoder-local notation outside allowed source files:",
        pattern,
        paste(names(hits), collapse = ", ")
      ))
    }
  }

  succeed()
})


test_that("generated man pages do not reintroduce forbidden aliases outside registry topics", {
  allowed_man_files <- normalizePath(
    c(
      file.path(pkg_root, "man", "ppgt_notations.Rd")
    ),
    mustWork = FALSE
  )

  man_targets <- setdiff(normalizePath(man_files, mustWork = FALSE), allowed_man_files)

  forbidden_patterns <- c(
    "\\\\widetilde\\{\\\\mathbf\\{Y\\}\\}",
    "\\bn_i\\b"
  )

  for (pattern in forbidden_patterns) {
    hits <- collect_matches(man_targets, pattern)
    if (length(hits) != 0L) {
      fail(paste(
        "Unexpected forbidden notation pattern in generated man pages:",
        pattern,
        paste(names(hits), collapse = ", ")
      ))
    }
  }

  succeed()
})
