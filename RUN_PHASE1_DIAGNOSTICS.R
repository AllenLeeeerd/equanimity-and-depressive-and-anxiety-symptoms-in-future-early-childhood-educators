# ============================================================
# v2 Phase 1A: diagnostics / assumptions
# Run this first with TEST_MODE <- TRUE in 00_config.R
# ============================================================

source("01_data_preparation.R")
source("02_attrition_ipcw.R")
source("03_measurement_invariance.R")
source("04_ega_community_check.R")
source("RUN_LATENT_CLPM.R")
source("RUN_SIMPLIFIED_LATENT_SENSITIVITY.R")
source("RUN_PHASE2_FINAL.R")
source("RUN_BOOTSTRAP_FIX_V2_6.R")
source("RUN_BOOTSTRAP_FIX_V2_6B_SELF_CONTAINED.R")
message(
  "/nV2 PHASE 1A FINISHED./n",
  "Do NOT run the latent CLPM or Phase 2 yet./n",
  "Zip 'revision_outputs_v2' and send it back to ChatGPT for review of ",
  "measurement invariance, EGA, straightlining QC, and IPCW."
)
