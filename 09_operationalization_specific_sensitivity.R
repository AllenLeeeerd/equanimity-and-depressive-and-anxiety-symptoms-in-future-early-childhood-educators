source("00_config.R")

if (isTRUE(TEST_MODE)) {
  stop(
    "Final operationalization-specific sensitivity requires TEST_MODE <- FALSE ",
    "so Model A and Model B are each independently bootstrapped with the final N_BOOT."
  )
}

infile <- file.path(
  OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv"
)
if (!file.exists(infile)) stop("Run 01_data_preparation.R first.")

dat <- read.csv(infile, check.names=FALSE)

out <- file.path(
  OUTPUT_DIR, "09_operationalization_specific_sensitivity"
)
dir.create(out, recursive=TRUE, showWarnings=FALSE)

nodes_A <- c("EAC","NRT",SYMPTOM_CODES)
nodes_B <- c("EMS","HIN",SYMPTOM_CODES)

message(
  "\nEstimating operationalization-specific Model A ",
  "(EAC + NRT + symptoms; N=", nrow(dat), ")..."
)

A <- estimate_clpn(
  dat=dat,
  node_codes=nodes_A,
  out_dir=out,
  covariates=PRIMARY_COVARIATES,
  weight_col=NULL,
  suffix="ModelA_EAC_NRT"
)

message(
  "Running an INDEPENDENT final bootstrap for Model A (N_BOOT=",
  N_BOOT, ")..."
)

boot_A <- run_network_bootstrap(
  A,
  out_dir=out,
  suffix="ModelA_EAC_NRT",
  n_boot=N_BOOT,
  n_cores=N_CORES
)

message(
  "\nEstimating operationalization-specific Model B ",
  "(EMS + HIN + symptoms; N=", nrow(dat), ")..."
)

B <- estimate_clpn(
  dat=dat,
  node_codes=nodes_B,
  out_dir=out,
  covariates=PRIMARY_COVARIATES,
  weight_col=NULL,
  suffix="ModelB_EMS_HIN"
)

# Deliberately use a different RNG stream before the second independent bootstrap.
set.seed(SEED + 10000L)

message(
  "Running an INDEPENDENT final bootstrap for Model B (N_BOOT=",
  N_BOOT, ")..."
)

boot_B <- run_network_bootstrap(
  B,
  out_dir=out,
  suffix="ModelB_EMS_HIN",
  n_boot=N_BOOT,
  n_cores=N_CORES
)

# Put the two CS tables side-by-side in one audit file.
cs_A <- data.frame(
  model="Model A (EAC/NRT)",
  metric=names(boot_A$cs),
  cs_coefficient=as.numeric(boot_A$cs)
)

cs_B <- data.frame(
  model="Model B (EMS/HIN)",
  metric=names(boot_B$cs),
  cs_coefficient=as.numeric(boot_B$cs)
)

cs_audit <- dplyr::bind_rows(cs_A, cs_B)

write.csv(
  cs_audit,
  file.path(out, "ModelA_ModelB_independent_CS_audit.csv"),
  row.names=FALSE
)

# Explicit numerical check: identical full CS vectors are flagged.
same_cs <- isTRUE(all.equal(
  as.numeric(boot_A$cs),
  as.numeric(boot_B$cs),
  tolerance=1e-12
))

write.csv(
  data.frame(
    check="Are the independently re-estimated full CS vectors identical?",
    identical=same_cs,
    note=ifelse(
      same_cs,
      paste0(
        "They are numerically identical after independent reruns. ",
        "Do not assume an error solely from equality; inspect the saved ",
        "bootstrap objects and report this transparently."
      ),
      paste0(
        "The independently re-estimated CS coefficients differ across models, ",
        "resolving the duplicated-CS concern."
      )
    )
  ),
  file.path(out, "ModelA_ModelB_CS_identity_check.csv"),
  row.names=FALSE
)

# Save everything
saveRDS(
  list(A=A, B=B, boot_A=boot_A, boot_B=boot_B),
  file.path(out, "operationalization_specific_models_with_bootstraps.rds")
)

message(
  "\n09_operationalization_specific_sensitivity FINAL complete.\n",
  "Model A and Model B have been independently re-estimated and bootstrapped."
)
