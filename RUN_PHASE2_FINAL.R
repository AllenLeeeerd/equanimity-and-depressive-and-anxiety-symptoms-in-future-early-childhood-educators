# ============================================================
# FINAL Phase 2 runner
# ============================================================

source("00_config.R")

if (isTRUE(TEST_MODE)) {
  stop(
    "Before final Phase 2, open 00_config.R and set TEST_MODE <- FALSE.\n",
    "Then rerun: source('RUN_PHASE2_FINAL.R')"
  )
}

if (N_BOOT < 1000L) {
  warning(
    "N_BOOT is below 1000. For the final manuscript run, 1000 is recommended."
  )
}

message(
  "\nFINAL PHASE 2 STARTING\n",
  "Expected primary matched N = ", EXPECTED_MATCHED_ANALYTIC_N, "\n",
  "Expected exact-match N = ", EXPECTED_EXACT_ANALYTIC_N, "\n",
  "N_BOOT = ", N_BOOT, "\n"
)

source("05_combined_clpn_primary.R")
source("06_scale_level_clpn.R")
source("08_exact_match_sensitivity.R")
source("09_operationalization_specific_sensitivity.R")
source("10_replicability_and_reporting_tables.R")

message(
  "\nFINAL PHASE 2 FINISHED.\n",
  "Zip the entire 'revision_outputs_v2' folder and send it back to ChatGPT."
)
