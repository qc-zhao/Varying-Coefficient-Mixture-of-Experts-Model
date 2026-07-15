#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !file.exists(args[[1L]])) {
  stop("Usage: audit_release_tarball.R path/to/VCMoE_0.2.0.tar.gz", call. = FALSE)
}

tarball <- normalizePath(args[[1L]], mustWork = TRUE)
members <- utils::untar(tarball, list = TRUE)
result_dir <- file.path(getwd(), "release-validation", "0.2.0", "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

checks <- list()
record <- function(check, passed, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

forbidden_members <- c(
  "/.git",
  "/dev/", "/cluster/", "/release-validation/", "/simulation_results/",
  "/results/", ".h5ad", ".RData", ".rda", ".tar", ".zip"
)
for (pattern in forbidden_members) {
  hits <- members[grepl(pattern, members, fixed = TRUE)]
  record(
    paste0("member:", pattern),
    !length(hits),
    paste(utils::head(hits, 10L), collapse = ";")
  )
}

rds_members <- members[grepl("\\.rds$", members, ignore.case = TRUE)]
nonstandard_rds <- setdiff(rds_members, "VCMoE/build/vignette.rds")
record(
  "member:nonstandard_rds",
  !length(nonstandard_rds),
  paste(utils::head(nonstandard_rds, 10L), collapse = ";")
)

extract_dir <- file.path(
  tempdir(),
  paste0("VCMoE-0.2.0-tarball-audit-", Sys.getpid())
)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
utils::untar(tarball, exdir = extract_dir)
root_candidates <- list.dirs(extract_dir, recursive = FALSE, full.names = TRUE)
if (length(root_candidates) != 1L) {
  stop("Tarball must contain exactly one package root.", call. = FALSE)
}
package_root <- root_candidates[[1L]]

surface <- c(
  file.path(package_root, "R"),
  file.path(package_root, "man"),
  file.path(package_root, "vignettes"),
  file.path(package_root, "NAMESPACE"),
  file.path(package_root, "README.md"),
  file.path(package_root, "NEWS.md"),
  file.path(package_root, "DESCRIPTION")
)
surface_files <- unlist(lapply(surface, function(path) {
  if (dir.exists(path)) {
    list.files(path, recursive = TRUE, full.names = TRUE)
  } else if (file.exists(path)) {
    path
  } else {
    character(0L)
  }
}), use.names = FALSE)
surface_text <- unlist(lapply(surface_files, readLines, warn = FALSE), use.names = FALSE)

removed_api_tokens <- c(
  "composite_alpha",
  "response_weights",
  "confirm_experimental",
  "experimental_joint_path_em",
  "allow_experimental_joint_path_inference",
  "ordered_u_starts",
  "legacy_gaussian_raw"
)
for (token in removed_api_tokens) {
  hits <- surface_files[vapply(surface_files, function(path) {
    any(grepl(token, readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1L))]
  record(
    paste0("removed_api:", token),
    !length(hits),
    paste(basename(hits), collapse = ";")
  )
}

all_text_files <- list.files(
  package_root,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.(R|Rd|Rmd|md|txt|yml|yaml|dcf|csv)$"
)
private_tokens <- c(
  "/Users/", "/scratch/", "GoogleDrive-", "q83zhao", "VCMoE-work-private"
)
for (token in private_tokens) {
  hits <- all_text_files[vapply(all_text_files, function(path) {
    any(grepl(token, readLines(path, warn = FALSE), fixed = TRUE))
  }, logical(1L))]
  record(
    paste0("private_text:", token),
    !length(hits),
    paste(basename(hits), collapse = ";")
  )
}

record(
  "version_0.2.0",
  identical(unname(read.dcf(file.path(package_root, "DESCRIPTION"))[1L, "Version"]), "0.2.0")
)
record(
  "joint_path_documented",
  any(grepl("joint_path_em", surface_text, fixed = TRUE))
)

audit <- do.call(rbind, checks)
audit$tarball <- basename(tarball)
output <- file.path(result_dir, "release_tarball_audit.csv")
utils::write.csv(audit, output, row.names = FALSE)
cat("Tarball audit:", output, "\n")
cat("Passed:", sum(audit$passed), "/", nrow(audit), "\n")
if (!all(audit$passed)) {
  print(audit[!audit$passed, , drop = FALSE])
  quit(status = 1L)
}
