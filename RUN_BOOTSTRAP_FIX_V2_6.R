# ============================================================
# Bootstrap correction rerun v2.6
#
# Reruns ONLY analyses whose bootnet resampling was affected by
# the RNG-reset bug:
#   05 primary combined CLPN + its bootstraps
#   06 scale-level CLPN + its bootstraps
#   09 Model A / Model B + independent bootstraps
#   10 reporting/replicability tables
#
# Does NOT rerun:
#   data cleaning, attrition/IPCW diagnostics, MI, EGA,
#   latent SEMs, or exact-match sensitivity.
# ============================================================

source("00_config.R")

if (isTRUE(TEST_MODE)) {
  stop(
    "Open 00_config.R and set TEST_MODE <- FALSE before the corrected final bootstrap rerun."
  )
}

if (N_BOOT < 1000L) {
  warning("N_BOOT < 1000. Final manuscript run should use at least 1000.")
}

message(
  "\nBOOTSTRAP CORRECTION v2.6 STARTING\n",
  "This rerun replaces the invalid bootstrap CI/selection/CS outputs produced ",
  "before the RNG fix.\n"
)

source("05_combined_clpn_primary.R")
source("06_scale_level_clpn.R")
source("09_operationalization_specific_sensitivity.R")
source("10_replicability_and_reporting_tables.R")

message(
  "\nBOOTSTRAP CORRECTION v2.6 FINISHED.\n",
  "Do NOT rerun Phase 1, EGA, latent SEMs, IPCW diagnostics, or exact-match.\n",
  "Please send back the updated result folders 05, 06, 09, and 10 ",
  "(or regenerate the Codex master Excel from revision_outputs_v2)."
)
