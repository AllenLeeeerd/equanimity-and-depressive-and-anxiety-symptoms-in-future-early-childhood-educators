source("00_config.R")

infile <- file.path(OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv")
if (!file.exists(infile)) stop("Run 01_data_preparation.R first.")
dat <- read.csv(infile, check.names=FALSE)

out <- file.path(OUTPUT_DIR, "05_combined_clpn_primary")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

nodes <- c(EQ_CODES, SYMPTOM_CODES)

message("Estimating primary 20-node combined CLPN (N=", nrow(dat), ")...")
primary <- estimate_clpn(
  dat=dat,
  node_codes=nodes,
  out_dir=out,
  covariates=PRIMARY_COVARIATES,
  weight_col=NULL,
  suffix="combined_primary_unweighted"
)

message("Running primary network bootstrap (N_BOOT = ", N_BOOT, ")...")
boot_primary <- run_network_bootstrap(
  primary, out_dir=out,
  suffix="combined_primary_unweighted",
  n_boot=N_BOOT, n_cores=N_CORES
)

ipcw_file <- file.path(OUTPUT_DIR, "02_attrition_ipcw", "matched_analytic_with_ipcw.csv")
if (file.exists(ipcw_file)) {
  datw <- read.csv(ipcw_file, check.names=FALSE)

  message("Estimating IPCW sensitivity CLPN (N=", nrow(datw), ")...")
  weighted <- estimate_clpn(
    dat=datw,
    node_codes=nodes,
    out_dir=out,
    covariates=PRIMARY_COVARIATES,
    weight_col="IPCW_TRIM",
    suffix="combined_ipcw_sensitivity"
  )

  comp <- compare_networks(
    primary$graph, weighted$graph,
    paste0("Primary unweighted N=", nrow(dat)),
    paste0("IPCW sensitivity N=", nrow(datw))
  )
  write.csv(comp, file.path(out, "primary_vs_ipcw_network_comparison.csv"), row.names=FALSE)
  saveRDS(weighted, file.path(out, "combined_ipcw_result.rds"))
} else {
  warning("IPCW file not found. Run 02_attrition_ipcw.R first.")
}

saveRDS(primary, file.path(out, "combined_primary_result.rds"))
message("05_combined_clpn_primary v2 complete: ", normalizePath(out))
