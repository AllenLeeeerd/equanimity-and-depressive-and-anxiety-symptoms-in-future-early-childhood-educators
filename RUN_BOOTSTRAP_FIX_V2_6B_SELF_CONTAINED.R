# ============================================================
# BOOTSTRAP FIX v2.6B — SELF-CONTAINED RUNNER
#
# PURPOSE
# -------
# The prior v2.6 repair was not actually loaded in the user's final run:
# updated bootstrap files were produced, but no v2.6 integrity-audit files
# existed, which means 00_config.R reloaded the old R/helpers_clpn.R.
#
# This script is deliberately self-contained:
# 1) it sources 00_config.R once;
# 2) it overrides the RNG-sensitive bootstrap functions IN MEMORY;
# 3) it runs the existing Phase-2 scripts WITHOUT letting them re-source
#    00_config.R and overwrite the patch;
# 4) it archives the known-invalid old bootstrap outputs;
# 5) it performs automatic bootstrap-integrity checks.
#
# Put THIS ONE FILE in the root equanimity_revision_R_v2 folder.
# Then run:
# source("RUN_BOOTSTRAP_FIX_V2_6B_SELF_CONTAINED.R")
# ============================================================

source("00_config.R")

if (isTRUE(TEST_MODE)) {
  stop(
    "Open 00_config.R and set TEST_MODE <- FALSE before running v2.6B."
  )
}
if (N_BOOT < 1000L) {
  stop(
    "N_BOOT is below 1000. Set N_BOOT <- 1000L in 00_config.R before running v2.6B."
  )
}

PATCH_VERSION <- "v2.6B-self-contained"

# ------------------------------------------------------------
# 1. Override make_foldid() without changing the OUTER RNG state
# ------------------------------------------------------------
make_foldid <- function(n, k = N_FOLDS, seed = SEED) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed + n))
  sample(rep(seq_len(k), length.out = n))
}

# Preflight: make_foldid must leave the outer RNG stream untouched.
set.seed(84721)
dummy1 <- runif(1)
dummy_fold <- make_foldid(381, 10, SEED)
after_patch <- runif(1)

set.seed(84721)
dummy2 <- runif(1)
expected_after <- runif(1)

if (!isTRUE(all.equal(dummy1, dummy2, tolerance = 0)) ||
    !isTRUE(all.equal(after_patch, expected_after, tolerance = 0))) {
  stop(
    "v2.6B preflight failed: make_foldid() is still altering the outer RNG state."
  )
}

message("v2.6B RNG preflight: PASS")

# ------------------------------------------------------------
# 2. Bootstrap integrity audit
# ------------------------------------------------------------
export_bootstrap_integrity_audit <- function(bootobj, out_dir, suffix) {
  bt <- bootobj$bootTable %>% dplyr::filter(type == "edge")

  per_edge <- bt %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_unique_values = dplyr::n_distinct(round(value, 12)),
      sd_boot = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      bootobj$sampleTable %>%
        dplyr::filter(type == "edge") %>%
        dplyr::select(id, sample_weight = value),
      by = "id"
    )

  nonzero <- per_edge %>% dplyr::filter(abs(sample_weight) > 1e-12)

  summary_tab <- data.frame(
    patch_version = PATCH_VERSION,
    suffix = suffix,
    n_edges = nrow(per_edge),
    n_nonzero_sample_edges = nrow(nonzero),
    median_unique_boot_values_nonzero =
      if (nrow(nonzero) > 0) median(nonzero$n_unique_values) else NA_real_,
    min_unique_boot_values_nonzero =
      if (nrow(nonzero) > 0) min(nonzero$n_unique_values) else NA_real_,
    max_unique_boot_values_nonzero =
      if (nrow(nonzero) > 0) max(nonzero$n_unique_values) else NA_real_,
    proportion_nonzero_edges_with_gt10_unique =
      if (nrow(nonzero) > 0) mean(nonzero$n_unique_values > 10) else NA_real_,
    median_boot_sd_nonzero =
      if (nrow(nonzero) > 0) median(nonzero$sd_boot, na.rm = TRUE) else NA_real_
  )

  write.csv(
    per_edge,
    file.path(out_dir, paste0("bootstrap_integrity_per_edge_", suffix, ".csv")),
    row.names = FALSE
  )
  write.csv(
    summary_tab,
    file.path(out_dir, paste0("bootstrap_integrity_summary_", suffix, ".csv")),
    row.names = FALSE
  )

  if (nrow(nonzero) > 0 &&
      summary_tab$proportion_nonzero_edges_with_gt10_unique < 0.50) {
    stop(
      "BOOTSTRAP INTEGRITY CHECK FAILED for ", suffix,
      ". Less than 50% of sample-nonzero edges have >10 unique bootstrap values. ",
      "Do not interpret bootstrap outputs. Send the integrity CSV to ChatGPT."
    )
  }

  message(
    "Bootstrap integrity PASS for ", suffix,
    ": median unique values among sample-nonzero edges = ",
    summary_tab$median_unique_boot_values_nonzero
  )

  invisible(summary_tab)
}

# ------------------------------------------------------------
# 3. Override run_network_bootstrap()
# ------------------------------------------------------------
run_network_bootstrap <- function(result, out_dir, suffix = "main",
                                  n_boot = N_BOOT, n_cores = N_CORES,
                                  seed_offset = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stats <- c(
    "edge",
    "outExpectedInfluence", "inExpectedInfluence",
    "outStrength", "inStrength",
    "bridgeExpectedInfluence", "bridgeStrength"
  )

  if (is.null(seed_offset)) {
    seed_offset <- sum(utf8ToInt(as.character(suffix))) * 100L
  }

  set.seed(as.integer(SEED + seed_offset))
  b_nonparam <- bootnet::bootnet(
    result$network,
    nBoots = n_boot,
    nCores = n_cores,
    type = "nonparametric",
    directed = TRUE,
    includeDiagonal = TRUE,
    communities = result$communities,
    statistics = stats,
    memorysaver = FALSE
  )

  # CRITICAL: audit immediately after nonparametric bootstrap, before
  # interpreting/exporting it.
  export_bootstrap_integrity_audit(
    b_nonparam, out_dir, suffix
  )

  set.seed(as.integer(SEED + seed_offset + 1L))
  b_case <- bootnet::bootnet(
    result$network,
    nBoots = n_boot,
    nCores = n_cores,
    type = "case",
    directed = TRUE,
    includeDiagonal = TRUE,
    communities = result$communities,
    statistics = c(
      "outExpectedInfluence", "inExpectedInfluence",
      "outStrength", "inStrength",
      "bridgeExpectedInfluence", "bridgeStrength"
    ),
    memorysaver = FALSE
  )

  saveRDS(
    b_nonparam,
    file.path(out_dir, paste0("bootstrap_nonparametric_", suffix, ".rds"))
  )
  saveRDS(
    b_case,
    file.path(out_dir, paste0("bootstrap_case_", suffix, ".rds"))
  )

  cs <- bootnet::corStability(b_case)
  write.csv(
    data.frame(metric = names(cs), cs_coefficient = as.numeric(cs)),
    file.path(out_dir, paste0("cs_coefficients_", suffix, ".csv")),
    row.names = FALSE
  )

  pdf(
    file.path(out_dir, paste0("edge_accuracy_", suffix, ".pdf")),
    width = 11, height = 9
  )
  print(plot(b_nonparam, "edge", labels = FALSE, order = "sample"))
  dev.off()

  pdf(
    file.path(out_dir, paste0("centrality_case_stability_", suffix, ".pdf")),
    width = 12, height = 9
  )
  print(plot(
    b_case,
    c(
      "outExpectedInfluence", "inExpectedInfluence",
      "outStrength", "inStrength",
      "bridgeExpectedInfluence", "bridgeStrength"
    ),
    facet = TRUE
  ))
  dev.off()

  write.csv(
    b_nonparam$sampleTable,
    file.path(out_dir, paste0("bootstrap_sample_table_", suffix, ".csv")),
    row.names = FALSE
  )
  write.csv(
    b_nonparam$bootTable,
    file.path(out_dir, paste0("bootstrap_full_table_", suffix, ".csv")),
    row.names = FALSE
  )

  export_bootstrap_ci(b_nonparam, out_dir, suffix)
  export_selection_frequency(b_nonparam, out_dir, suffix)
  export_equanimity_node_difference_tests(b_nonparam, out_dir, suffix)

  invisible(list(nonparametric = b_nonparam, case = b_case, cs = cs))
}

# ------------------------------------------------------------
# 4. Archive old invalid bootstrap-derived outputs
# ------------------------------------------------------------
archive_root <- file.path(
  OUTPUT_DIR,
  paste0("_invalid_bootstrap_archive_pre_", PATCH_VERSION)
)
dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)

archive_patterns <- c(
  "^bootstrap_",
  "^cs_coefficients_",
  "^edge_accuracy_",
  "^centrality_case_stability_",
  "^equanimity_node_difference_tests_",
  "^ModelA_ModelB_independent_CS_audit\\.csv$",
  "^ModelA_ModelB_CS_identity_check\\.csv$"
)

archive_folder <- function(subdir) {
  src <- file.path(OUTPUT_DIR, subdir)
  if (!dir.exists(src)) return(invisible(NULL))

  files <- list.files(src, full.names = TRUE)
  keep <- vapply(
    basename(files),
    function(x) any(vapply(archive_patterns, grepl, logical(1), x = x)),
    logical(1)
  )
  files <- files[keep]
  if (length(files) == 0) return(invisible(NULL))

  dst <- file.path(archive_root, subdir)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)

  for (f in files) {
    ok <- file.copy(f, file.path(dst, basename(f)), overwrite = TRUE)
    if (ok) unlink(f)
  }
}

archive_folder("05_combined_clpn_primary")
archive_folder("06_scale_level_clpn")
archive_folder("09_operationalization_specific_sensitivity")

rep_file <- file.path(
  OUTPUT_DIR, "10_replicability_and_reporting_tables",
  "bootstrap_edge_replicability_summary.csv"
)
if (file.exists(rep_file)) {
  dst <- file.path(archive_root, "10_replicability_and_reporting_tables")
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  if (file.copy(rep_file, file.path(dst, basename(rep_file)), overwrite = TRUE)) {
    unlink(rep_file)
  }
}

# ------------------------------------------------------------
# 5. Run existing scripts WITHOUT re-sourcing 00_config.R
#    (which would reload the old helper and undo this patch)
# ------------------------------------------------------------
run_script_without_config <- function(file) {
  if (!file.exists(file)) stop("Required script not found: ", file)

  txt <- readLines(file, warn = FALSE, encoding = "UTF-8")
  txt <- txt[
    !grepl(
      '^\\s*source\\s*\\(\\s*["\\\']00_config\\.R["\\\']\\s*\\)\\s*$',
      txt
    )
  ]

  message("\n--- Running ", file, " under ", PATCH_VERSION, " ---")
  eval(parse(text = txt, keep.source = TRUE), envir = .GlobalEnv)
}

# Write a marker so later audits can verify the correct patch was loaded.
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    paste0("patch_version=", PATCH_VERSION),
    paste0("run_time=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "rng_preflight=PASS",
    "runner=self-contained; component scripts executed without re-sourcing 00_config.R"
  ),
  file.path(OUTPUT_DIR, "BOOTSTRAP_FIX_V2_6B_MARKER.txt")
)

message(
  "\n============================================================\n",
  "BOOTSTRAP FIX v2.6B STARTING\n",
  "Primary / Scale / Model A / Model B bootstrap outputs will be rebuilt.\n",
  "Old invalid bootstrap outputs were archived under:\n",
  archive_root, "\n",
  "============================================================\n"
)

run_script_without_config("05_combined_clpn_primary.R")
run_script_without_config("06_scale_level_clpn.R")
run_script_without_config("09_operationalization_specific_sensitivity.R")
run_script_without_config("10_replicability_and_reporting_tables.R")

message(
  "\n============================================================\n",
  "BOOTSTRAP FIX v2.6B FINISHED SUCCESSFULLY.\n",
  "Confirm that revision_outputs_v2 contains:\n",
  "  BOOTSTRAP_FIX_V2_6B_MARKER.txt\n",
  "and bootstrap_integrity_summary_*.csv in 05, 06, and 09.\n",
  "Then ask Codex to rebuild the master Excel.\n",
  "============================================================\n"
)
