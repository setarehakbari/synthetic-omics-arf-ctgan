############################################################
# lock_dataset_provenance.R
#
# Verifies the exact synthetic datasets used for the revised
# manuscript and writes a provenance manifest.
#
# Run this separately; it does NOT rerun any evaluation.
############################################################

base_dir <- path.expand("~/Desktop/My paper")
out_dir <- file.path(base_dir, "reviewer_revision", "all_generators", "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  ARF   = file.path(base_dir, "synthetic_rna_arf.csv"),
  CTGAN = file.path(base_dir, "synthetic_ctgan.csv"),
  TVAE  = file.path(base_dir, "synthetic_rna_tvae.csv")
)

expected_md5 <- c(
  ARF   = "24ee0d628ff7a1dfd5f6ebf66a2bd358",
  CTGAN = "630050778c8a36eb65f507d42c607c7c",
  TVAE  = "279052ce8e8e9bd154e526c6d7814a94"
)

expected_rows <- c(ARF = 200, CTGAN = 200, TVAE = 200)
expected_cols <- c(ARF = 1784, CTGAN = 1784, TVAE = 1784)

missing <- names(files)[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing canonical file(s): ", paste(missing, collapse = ", "))
}

manifest <- do.call(
  rbind,
  lapply(names(files), function(nm) {
    f <- files[[nm]]
    d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

    if (!"class" %in% names(d)) {
      stop(nm, " file has no 'class' column.")
    }

    cls <- table(as.character(d$class))
    md5_now <- unname(tools::md5sum(f))

    data.frame(
      dataset = nm,
      filename = basename(f),
      rows = nrow(d),
      columns = ncol(d),
      class_0 = if ("0" %in% names(cls)) unname(cls["0"]) else NA_integer_,
      class_1 = if ("1" %in% names(cls)) unname(cls["1"]) else NA_integer_,
      md5 = md5_now,
      expected_md5 = expected_md5[[nm]],
      md5_match = identical(md5_now, expected_md5[[nm]]),
      dimensions_match =
        nrow(d) == expected_rows[[nm]] &&
        ncol(d) == expected_cols[[nm]],
      stringsAsFactors = FALSE
    )
  })
)

print(manifest)

if (!all(manifest$md5_match)) {
  bad <- manifest$dataset[!manifest$md5_match]
  stop(
    "PROVENANCE CHECK FAILED: MD5 mismatch for: ",
    paste(bad, collapse = ", "),
    "\nDo not use the affected dataset for manuscript results until resolved."
  )
}

if (!all(manifest$dimensions_match)) {
  bad <- manifest$dataset[!manifest$dimensions_match]
  stop(
    "PROVENANCE CHECK FAILED: unexpected dimensions for: ",
    paste(bad, collapse = ", ")
  )
}

write.csv(
  manifest,
  file.path(out_dir, "dataset_provenance_manifest.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("PROVENANCE CHECK PASSED.\n")
cat("Canonical manuscript datasets are locked and verified.\n")
cat("Manifest written to:\n")
cat(file.path(out_dir, "dataset_provenance_manifest.csv"), "\n")
cat("============================================================\n")
