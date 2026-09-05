source("00_config.R")

infile <- file.path(OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv")
if (!file.exists(infile)) stop("Run 01_data_preparation.R first.")
dat <- read.csv(infile, check.names=FALSE)

if (nrow(dat) != EXPECTED_MATCHED_ANALYTIC_N) stop("Analytic N mismatch before invariance.")

out <- file.path(OUTPUT_DIR, "03_measurement_invariance")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

# Preserve the untouched analytic data for all non-invariance analyses.
mi_dat <- dat

# Only for PHQ longitudinal invariance:
# T1 depressed mood has a rarely endorsed category 3 that is absent at T2.
# Harmonize the upper categories to 2+ at BOTH waves. This does NOT alter CLPN data.
if (MI_COLLAPSE_PHQ_DEP_TOP) {
  mi_dat$T1_DEP <- ifelse(mi_dat$T1_DEP >= 2, 2, mi_dat$T1_DEP)
  mi_dat$T2_DEP <- ifelse(mi_dat$T2_DEP >= 2, 2, mi_dat$T2_DEP)
}

# Record the category harmonization explicitly.
dep_counts_before <- bind_rows(
  data.frame(wave="T1", category=sort(unique(dat$T1_DEP)),
             n=sapply(sort(unique(dat$T1_DEP)), function(x) sum(dat$T1_DEP==x))),
  data.frame(wave="T2", category=sort(unique(dat$T2_DEP)),
             n=sapply(sort(unique(dat$T2_DEP)), function(x) sum(dat$T2_DEP==x)))
)
dep_counts_after <- bind_rows(
  data.frame(wave="T1", category=sort(unique(mi_dat$T1_DEP)),
             n=sapply(sort(unique(mi_dat$T1_DEP)), function(x) sum(mi_dat$T1_DEP==x))),
  data.frame(wave="T2", category=sort(unique(mi_dat$T2_DEP)),
             n=sapply(sort(unique(mi_dat$T2_DEP)), function(x) sum(mi_dat$T2_DEP==x)))
)
write.csv(dep_counts_before, file.path(out, "PHQ_DEP_categories_before_MI_harmonization.csv"),
          row.names=FALSE)
write.csv(dep_counts_after, file.path(out, "PHQ_DEP_categories_after_MI_harmonization.csv"),
          row.names=FALSE)

structures <- list(
  ES14 = list(
    EAC=paste0("ES",sprintf("%02d",1:8)),
    NRT=paste0("ES",sprintf("%02d",11:16),"R")
  ),
  EQUAS = list(
    EMS=c("EQ01","EQ02","EQ05","EQ08"),
    HIN=paste0("EQ",sprintf("%02d",9:14),"R")
  ),
  PHQ9 = list(
    PHQF=c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI")
  ),
  GAD7 = list(
    GADF=c("ANX","CTL","WRY","RLX","RST","IRR","AFR")
  )
)

results <- list()
for (nm in names(structures)) {
  message("Running ordered longitudinal invariance: ", nm)
  results[[nm]] <- run_longitudinal_invariance(
    dat=mi_dat,
    factor_items=structures[[nm]],
    scale_name=nm,
    out_dir=out
  )
}

saveRDS(results, file.path(out, "all_measurement_invariance_results.rds"))

message(
  "03_measurement_invariance v2 complete.\n",
  "STOP and review the fit tables before specifying the latent CLPM. ",
  "Do not free parameters or impose partial invariance automatically."
)
