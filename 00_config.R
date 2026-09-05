# ============================================================
# Equanimity longitudinal network revision
# v2 central configuration
# ============================================================

# ---- 1. Paths ----
DATA_DIR <- "."
RAW_FILE <- file.path(DATA_DIR, "筛选掉无效数据以后的总原始数据.xlsx")
PAIRED_FILE <- file.path(DATA_DIR, "paired_wide_demographics.xlsx")

OUTPUT_DIR <- file.path(DATA_DIR, "revision_outputs_v2")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Confirmed pre-QC sample-flow constants ----
RAW_T1_N <- 611L
INITIAL_VALID_T1_N <- 551L
RAW_T2_N <- 500L
INITIAL_VALID_T2_N <- 448L

MATCHED_PRE_QC_N <- 392L
EXACT_MATCH_PRE_QC_N <- 313L

# Strict cross-instrument straightlining was approved for exclusion.
# Expected counts are asserted after the data are reconstructed.
EXPECTED_T1_STRICT_ALL <- 7L
EXPECTED_T2_STRICT_ALL <- 7L
EXPECTED_T1_STRICT_MATCHED <- 5L
EXPECTED_T2_STRICT_MATCHED <- 6L
EXPECTED_MATCHED_STRICT_UNION <- 11L

EXPECTED_BASELINE_ELIGIBLE_N <- 544L
EXPECTED_MATCHED_ANALYTIC_N <- 381L
EXPECTED_EXACT_ANALYTIC_N <- 303L

# ---- 3. Analysis settings ----
SEED <- 20260831L
N_FOLDS <- 10L
LAMBDA_RULE <- "lambda.min"

# First run diagnostics with TEST_MODE <- TRUE.
TEST_MODE <- FALSE
N_BOOT <- if (TEST_MODE) 50L else 1000L
N_EGA_BOOT <- if (TEST_MODE) 50L else 500L
N_CORES <- 1L
PLOT_MIN_EDGE <- 0.03

# Stabilized IPCW truncation percentiles
IPCW_TRIM <- c(0.01, 0.99)

# ---- 4. Primary analytic decisions ----
# Mindfulness-experience is intentionally excluded from all revised models.
# Primary covariates are age and gender only.
PRIMARY_COVARIATES <- c("AGE_T1", "GEN_T1")

# ---- 5. Ordered longitudinal invariance ----
# The PHQ depressed-mood item has category 3 endorsed at T1 but not T2.
# Only for measurement-invariance analysis, categories 2 and 3 are harmonized
# to a single upper category (0 / 1 / 2+). Network analyses retain 0-3 coding.
MI_COLLAPSE_PHQ_DEP_TOP <- TRUE

# ---- 6. Packages ----
REQUIRED_PACKAGES <- c(
  "readxl", "dplyr", "tidyr", "purrr", "stringr", "tibble",
  "ggplot2", "glmnet", "bootnet", "qgraph", "networktools",
  "psych", "effectsize", "lavaan", "semTools", "EGAnet", "mclust"
)

missing_pkgs <- REQUIRED_PACKAGES[
  !vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))"
  )
}

set.seed(SEED)

source(file.path("R", "helpers_data.R"))
source(file.path("R", "helpers_clpn.R"))
source(file.path("R", "helpers_sem.R"))
