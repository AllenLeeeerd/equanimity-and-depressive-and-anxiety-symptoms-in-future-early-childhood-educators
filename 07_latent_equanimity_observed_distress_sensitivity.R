# ============================================================
# 07_latent_equanimity_observed_distress_sensitivity.R
# Simplified latent sensitivity model
#
# Four equanimity dimensions are latent factors.
# PHQ-9 and GAD-7 are observed total scores.
# FIML uses the baseline-eligible cohort (expected N = 544).
# ============================================================

source("00_config.R")

infile <- file.path(
  OUTPUT_DIR, "01_data_preparation",
  "full_T1_eligible_with_T2_missing.csv"
)
if (!file.exists(infile)) {
  stop("Run v2 data preparation first.")
}

dat <- read.csv(infile, check.names = FALSE)

if (nrow(dat) != EXPECTED_BASELINE_ELIGIBLE_N) {
  stop(
    "Expected baseline-eligible N=", EXPECTED_BASELINE_ELIGIBLE_N,
    "; found ", nrow(dat)
  )
}

out <- file.path(
  OUTPUT_DIR,
  "07_latent_equanimity_observed_distress_sensitivity"
)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Covariates
# ------------------------------------------------------------
dat$AGE_Z <- as.numeric(scale(dat$T1_AGE))

g <- as.numeric(dat$T1_GEN)
if (all(sort(unique(g[!is.na(g)])) %in% c(1, 2))) {
  dat$GEN01 <- ifelse(g == 2, 1, 0)   # 0 = male, 1 = female
} else {
  stop("Unexpected gender coding.")
}

# ------------------------------------------------------------
# 2. Treat invalid / unavailable T2 data as missing
# ------------------------------------------------------------
fac_items <- list(
  EAC = c("ES01","ES02","ES03","ES04","ES05","ES06","ES07","ES08"),
  NRT = c("ES11R","ES12R","ES13R","ES14R","ES15R","ES16R"),
  EMS = c("EQ01","EQ02","EQ05","EQ08"),
  HIN = c("EQ09R","EQ10R","EQ11R","EQ12R","EQ13R","EQ14R")
)

all_eq_items <- unique(unlist(fac_items))
t2_eq_items <- paste0("T2_", all_eq_items)

invalid_t2 <- dat$RETAINED != 1

dat[invalid_t2, t2_eq_items] <- NA
dat$T2_PHQ_TOTAL[invalid_t2] <- NA
dat$T2_GAD_TOTAL[invalid_t2] <- NA

# ------------------------------------------------------------
# 3. Longitudinal strong-invariance measurement model
#    Equal loadings + equal item intercepts across waves.
#    Residual variances remain free.
# ------------------------------------------------------------
make_factor_lines <- function(fac, items) {
  first <- items[1]
  rest <- items[-1]

  t1_rhs <- paste(
    c(
      paste0("1*T1_", first),
      paste0("l_", fac, "_", seq_along(rest) + 1, "*T1_", rest)
    ),
    collapse = " + "
  )

  t2_rhs <- paste(
    c(
      paste0("1*T2_", first),
      paste0("l_", fac, "_", seq_along(rest) + 1, "*T2_", rest)
    ),
    collapse = " + "
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

# Baseline latent means fixed to zero; T2 latent means free.
mean_lines <- c(
  paste0("T1_", names(fac_items), " ~ 0*1"),
  paste0("T2_", names(fac_items), " ~ 1")
)

# ------------------------------------------------------------
# 4. Structural cross-lagged model
#
# T2 outcomes:
#   EAC, NRT, EMS, HIN, PHQ total, GAD total
#
# Each is predicted by:
#   T1 EAC, NRT, EMS, HIN, PHQ total, GAD total, age, gender
# ------------------------------------------------------------
t1_predictors <- c(
  "T1_EAC","T1_NRT","T1_EMS","T1_HIN",
  "T1_PHQ_TOTAL","T1_GAD_TOTAL",
  "AGE_Z","GEN01"
)

t2_outcomes <- c(
  "T2_EAC","T2_NRT","T2_EMS","T2_HIN",
  "T2_PHQ_TOTAL","T2_GAD_TOTAL"
)

structural_lines <- unlist(lapply(t2_outcomes, function(outcome) {
  paste0(outcome, " ~ ", paste(t1_predictors, collapse = " + "))
}))

# Correlated residuals among the six T2 outcomes.
t2_resid_cov <- character()
for (i in 1:(length(t2_outcomes) - 1)) {
  for (j in (i + 1):length(t2_outcomes)) {
    t2_resid_cov <- c(
      t2_resid_cov,
      paste0(t2_outcomes[i], " ~~ ", t2_outcomes[j])
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
  collapse = "\n"
)

writeLines(
  model_syntax,
  file.path(out, "latent_eq_observed_distress_model_syntax.txt")
)

# ------------------------------------------------------------
# 5. Estimate with robust ML + FIML
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

saveRDS(
  fit,
  file.path(out, "latent_eq_observed_distress_fiml_mlr.rds")
)

conv <- lavaan::lavInspect(fit, "converged")

writeLines(
  ifelse(isTRUE(conv), "TRUE", "FALSE"),
  file.path(out, "convergence.txt")
)

if (!isTRUE(conv)) {
  stop(
    "Simplified latent sensitivity model did not converge. ",
    "Do not modify the model manually. Send the output folder to ChatGPT."
  )
}

# ------------------------------------------------------------
# 6. Fit indices
# ------------------------------------------------------------
fm <- lavaan::fitMeasures(fit)

getfm <- function(candidates) {
  for (nm in candidates) {
    if (nm %in% names(fm) && is.finite(fm[[nm]])) {
      return(unname(fm[[nm]]))
    }
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
  rmsea_ci_low = getfm(c(
    "rmsea.ci.lower.robust",
    "rmsea.ci.lower.scaled",
    "rmsea.ci.lower"
  )),
  rmsea_ci_high = getfm(c(
    "rmsea.ci.upper.robust",
    "rmsea.ci.upper.scaled",
    "rmsea.ci.upper"
  )),
  srmr = getfm(c("srmr")),
  aic = getfm(c("aic")),
  bic = getfm(c("bic"))
)

write.csv(
  fit_table,
  file.path(out, "latent_eq_observed_distress_fit.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Structural paths
# ------------------------------------------------------------
pe <- lavaan::parameterEstimates(
  fit,
  standardized = TRUE,
  ci = TRUE
)

struct <- pe[
  pe$op == "~" &
    pe$lhs %in% t2_outcomes,
  c(
    "lhs","op","rhs","est","se","z","pvalue",
    "ci.lower","ci.upper","std.all"
  )
]

get_construct <- function(x) {
  x <- sub("^T1_", "", x)
  x <- sub("^T2_", "", x)
  x
}

struct$path_type <- "covariate"

is_t1 <- grepl("^T1_", struct$rhs)
same_construct <- is_t1 &
  get_construct(struct$lhs) == get_construct(struct$rhs)

struct$path_type[is_t1] <- "cross_lagged"
struct$path_type[same_construct] <- "autoregressive"

write.csv(
  struct,
  file.path(out, "latent_eq_observed_distress_structural_paths.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Focused equanimity <-> distress cross-lagged paths
# ------------------------------------------------------------
eq_t1 <- c("T1_EAC","T1_NRT","T1_EMS","T1_HIN")
eq_t2 <- c("T2_EAC","T2_NRT","T2_EMS","T2_HIN")

dist_t1 <- c("T1_PHQ_TOTAL","T1_GAD_TOTAL")
dist_t2 <- c("T2_PHQ_TOTAL","T2_GAD_TOTAL")

focus <- struct[
  struct$path_type == "cross_lagged" & (
    (struct$lhs %in% dist_t2 & struct$rhs %in% eq_t1) |
    (struct$lhs %in% eq_t2 & struct$rhs %in% dist_t1)
  ),
]

write.csv(
  focus,
  file.path(out, "latent_eq_observed_distress_focus_paths.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. R-squared and factor correlations
# ------------------------------------------------------------
r2 <- lavaan::inspect(fit, "r2")
r2tab <- data.frame(
  variable = names(r2),
  r2 = as.numeric(r2)
)

write.csv(
  r2tab,
  file.path(out, "latent_eq_observed_distress_R2.csv"),
  row.names = FALSE
)

cor_lv <- tryCatch(
  lavaan::lavInspect(fit, "cor.lv"),
  error = function(e) NULL
)

if (!is.null(cor_lv)) {
  write.csv(
    cor_lv,
    file.path(out, "latent_eq_factor_correlations.csv")
  )
}

write.csv(
  pe,
  file.path(out, "latent_eq_observed_distress_all_parameters.csv"),
  row.names = FALSE
)

sink(file.path(out, "sessionInfo.txt"))
print(sessionInfo())
sink()

message(
  "\nSIMPLIFIED LATENT SENSITIVITY MODEL FINISHED.\n",
  "Please zip only the folder:\n",
  "revision_outputs_v2/07_latent_equanimity_observed_distress_sensitivity\n",
  "and send it back to ChatGPT before running Phase 2."
)
