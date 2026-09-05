# ============================================================
# 07_latent_clpm_sensitivity.R
# Full-information latent-variable cross-lagged sensitivity model
# ============================================================

source("00_config.R")

infile <- file.path(
  OUTPUT_DIR, "01_data_preparation",
  "full_T1_eligible_with_T2_missing.csv"
)
if (!file.exists(infile)) {
  stop("Run v2 data preparation first.")
}

dat <- read.csv(infile, check.names=FALSE)
if (nrow(dat) != EXPECTED_BASELINE_ELIGIBLE_N) {
  stop("Expected baseline-eligible N=", EXPECTED_BASELINE_ELIGIBLE_N,
       "; found ", nrow(dat))
}

out <- file.path(OUTPUT_DIR, "07_latent_clpm_sensitivity")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

# ------------------------------------------------------------
# 1. Prepare covariates and invalidate all T2 indicators for
#    participants without a valid T2 wave.
# ------------------------------------------------------------
dat$AGE_Z <- as.numeric(scale(dat$T1_AGE))

g <- as.numeric(dat$T1_GEN)
if (all(sort(unique(g[!is.na(g)])) %in% c(1,2))) {
  dat$GEN01 <- ifelse(g == 2, 1, 0)  # 0 male, 1 female
} else {
  stop("Unexpected gender coding.")
}

# Six latent constructs at both waves
fac_items <- list(
  EAC = c("ES01","ES02","ES03","ES04","ES05","ES06","ES07","ES08"),
  NRT = c("ES11R","ES12R","ES13R","ES14R","ES15R","ES16R"),
  EMS = c("EQ01","EQ02","EQ05","EQ08"),
  HIN = c("EQ09R","EQ10R","EQ11R","EQ12R","EQ13R","EQ14R"),
  DEP = c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI"),
  ANX = c("ANX","CTL","WRY","RLX","RST","IRR","AFR")
)

all_base_items <- unique(unlist(fac_items))
t2_obs <- paste0("T2_", all_base_items)

# Matched participants whose T2 was rejected for strict straightlining,
# plus all other non-retained T1 cases, are treated as having no valid T2.
invalid_t2 <- dat$RETAINED != 1
dat[invalid_t2, t2_obs] <- NA

# ------------------------------------------------------------
# 2. Build a strong-invariance measurement model.
#    The ordered WLSMV invariance analysis showed negligible
#    deterioration when adding longitudinal loading/intercept
#    constraints for ES-14, EQUA-S, PHQ-9 and GAD-7.
#
#    FIML is therefore implemented under MLR with equal loadings
#    and equal item intercepts across waves.
# ------------------------------------------------------------
make_factor_lines <- function(fac, items) {
  first <- items[1]
  rest <- items[-1]

  t1_rhs <- paste(
    c(paste0("1*", "T1_", first),
      paste0("l_", fac, "_", seq_along(rest)+1, "*T1_", rest)),
    collapse=" + "
  )
  t2_rhs <- paste(
    c(paste0("1*", "T2_", first),
      paste0("l_", fac, "_", seq_along(rest)+1, "*T2_", rest)),
    collapse=" + "
  )

  c(
    paste0("T1_", fac, " =~ ", t1_rhs),
    paste0("T2_", fac, " =~ ", t2_rhs)
  )
}

make_intercept_lines <- function(fac, items) {
  unlist(lapply(seq_along(items), function(i) {
    c(
      paste0("T1_", items[i], " ~ int_", fac, "_", i, "*1"),
      paste0("T2_", items[i], " ~ int_", fac, "_", i, "*1")
    )
  }))
}

make_residual_pair_lines <- function(items) {
  paste0("T1_", items, " ~~ T2_", items)
}

measurement_lines <- unlist(lapply(names(fac_items), function(f) {
  make_factor_lines(f, fac_items[[f]])
}))

intercept_lines <- unlist(lapply(names(fac_items), function(f) {
  make_intercept_lines(f, fac_items[[f]])
}))

residual_pair_lines <- unlist(lapply(fac_items, make_residual_pair_lines))

# Baseline factor means fixed to zero; T2 latent intercepts are free.
mean_lines <- c(
  paste0("T1_", names(fac_items), " ~ 0*1"),
  paste0("T2_", names(fac_items), " ~ 1")
)

# ------------------------------------------------------------
# 3. Structural model:
#    every T2 latent construct is regressed on all six T1
#    constructs plus age and gender.
# ------------------------------------------------------------
fac_names <- names(fac_items)

structural_lines <- unlist(lapply(fac_names, function(outcome) {
  rhs <- c(paste0("T1_", fac_names), "AGE_Z", "GEN01")
  paste0("T2_", outcome, " ~ ", paste(rhs, collapse=" + "))
}))

# Allow residual correlations among the six T2 latent factors.
t2_resid_cov <- character()
for (i in 1:(length(fac_names)-1)) {
  for (j in (i+1):length(fac_names)) {
    t2_resid_cov <- c(
      t2_resid_cov,
      paste0("T2_", fac_names[i], " ~~ T2_", fac_names[j])
    )
  }
}

model_syntax <- paste(
  c(
    measurement_lines,
    intercept_lines,
    residual_pair_lines,
    mean_lines,
    structural_lines,
    t2_resid_cov
  ),
  collapse="\n"
)

writeLines(model_syntax, file.path(out, "latent_clpm_model_syntax.txt"))

# ------------------------------------------------------------
# 4. Estimate with robust ML + FIML.
# ------------------------------------------------------------
fit <- lavaan::sem(
  model = model_syntax,
  data = dat,
  estimator = "MLR",
  missing = "fiml",
  meanstructure = TRUE,
  fixed.x = FALSE,
  auto.var = TRUE,
  auto.cov.lv.x = TRUE,
  auto.cov.y = FALSE,
  auto.fix.first = FALSE,
  warn = TRUE
)

saveRDS(fit, file.path(out, "latent_clpm_fiml_mlr.rds"))

conv <- lavaan::lavInspect(fit, "converged")
if (!isTRUE(conv)) {
  writeLines("FALSE", file.path(out, "CONVERGENCE_FAILED.txt"))
  stop(
    "Latent CLPM did not converge. Do not modify the model manually. ",
    "Send the 07_latent_clpm_sensitivity output folder back to ChatGPT."
  )
}

# ------------------------------------------------------------
# 5. Fit indices
# ------------------------------------------------------------
fm <- lavaan::fitMeasures(fit)
getfm <- function(candidates) {
  for (nm in candidates) {
    if (nm %in% names(fm) && is.finite(fm[[nm]])) return(unname(fm[[nm]]))
  }
  NA_real_
}

fit_table <- data.frame(
  n_total = nrow(dat),
  n_with_valid_T2 = sum(dat$RETAINED == 1),
  chisq = getfm(c("chisq.scaled","chisq")),
  df = getfm(c("df.scaled","df")),
  pvalue = getfm(c("pvalue.scaled","pvalue")),
  cfi = getfm(c("cfi.robust","cfi.scaled","cfi")),
  tli = getfm(c("tli.robust","tli.scaled","tli")),
  rmsea = getfm(c("rmsea.robust","rmsea.scaled","rmsea")),
  rmsea_ci_low = getfm(c("rmsea.ci.lower.robust","rmsea.ci.lower.scaled","rmsea.ci.lower")),
  rmsea_ci_high = getfm(c("rmsea.ci.upper.robust","rmsea.ci.upper.scaled","rmsea.ci.upper")),
  srmr = getfm(c("srmr")),
  aic = getfm(c("aic")),
  bic = getfm(c("bic"))
)
write.csv(fit_table, file.path(out, "latent_clpm_fit.csv"), row.names=FALSE)

# ------------------------------------------------------------
# 6. Standardized structural paths
# ------------------------------------------------------------
pe <- lavaan::parameterEstimates(
  fit,
  standardized=TRUE,
  ci=TRUE
)

struct <- pe[
  pe$op == "~" &
    grepl("^T2_(EAC|NRT|EMS|HIN|DEP|ANX)$", pe$lhs),
  c("lhs","op","rhs","est","se","z","pvalue","ci.lower","ci.upper","std.all")
]

struct$path_type <- ifelse(
  grepl("^T1_", struct$rhs),
  ifelse(
    sub("^T2_","",struct$lhs) == sub("^T1_","",struct$rhs),
    "autoregressive",
    "cross_lagged"
  ),
  "covariate"
)

write.csv(
  struct,
  file.path(out, "latent_clpm_structural_paths.csv"),
  row.names=FALSE
)

# Focused table: equanimity <-> distress cross-lagged paths.
eq <- c("EAC","NRT","EMS","HIN")
dist <- c("DEP","ANX")

focus <- struct[
  struct$path_type == "cross_lagged" & (
    (sub("^T2_","",struct$lhs) %in% dist &
       sub("^T1_","",struct$rhs) %in% eq) |
    (sub("^T2_","",struct$lhs) %in% eq &
       sub("^T1_","",struct$rhs) %in% dist)
  ),
]
write.csv(
  focus,
  file.path(out, "latent_clpm_equanimity_distress_crosslagged_paths.csv"),
  row.names=FALSE
)

# R-squared
r2 <- lavaan::inspect(fit, "r2")
r2tab <- data.frame(variable=names(r2), r2=as.numeric(r2))
write.csv(r2tab, file.path(out, "latent_clpm_R2.csv"), row.names=FALSE)

# Full standardized parameter table
write.csv(pe, file.path(out, "latent_clpm_all_parameters.csv"), row.names=FALSE)

# Covariance matrix diagnostics
cov_lv <- tryCatch(
  lavaan::lavInspect(fit, "cor.lv"),
  error=function(e) NULL
)
if (!is.null(cov_lv)) {
  write.csv(cov_lv, file.path(out, "latent_factor_correlations.csv"))
}

# Session info
sink(file.path(out, "sessionInfo.txt"))
print(sessionInfo())
sink()

message(
  "\nLATENT CLPM FINISHED.\n",
  "Please zip only the folder 'revision_outputs_v2/07_latent_clpm_sensitivity' ",
  "and send it back to ChatGPT before running Phase 2."
)
