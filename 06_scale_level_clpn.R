source("00_config.R")

infile <- file.path(OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv")
if (!file.exists(infile)) stop("Run 01_data_preparation.R first.")
dat <- read.csv(infile, check.names=FALSE)

out <- file.path(OUTPUT_DIR, "06_scale_level_clpn")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

nodes <- c(EQ_CODES, "PHQ_TOTAL", "GAD_TOTAL")

message("Estimating 6-node scale-level complementary CLPN (N=", nrow(dat), ")...")
res <- estimate_clpn(
  dat=dat,
  node_codes=nodes,
  out_dir=out,
  covariates=PRIMARY_COVARIATES,
  weight_col=NULL,
  suffix="scale_level_6node"
)

message("Running scale-level bootstrap...")
boot <- run_network_bootstrap(
  res, out_dir=out,
  suffix="scale_level_6node",
  n_boot=N_BOOT, n_cores=N_CORES
)

saveRDS(res, file.path(out, "scale_level_result.rds"))
message("06_scale_level_clpn v2 complete: ", normalizePath(out))
