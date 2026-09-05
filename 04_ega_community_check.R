# ============================================================
# EGA community check PATCH v2.2
# Fixes EGAnet seed-type error on Windows / current EGAnet build.
#
# Replace the previous 04_ega_community_check.R with this file.
# Then run ONLY:
# source("04_ega_community_check.R")
# ============================================================

source("00_config.R")

infile <- file.path(OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv")
if (!file.exists(infile)) stop("Run 01_data_preparation.R first.")
dat <- read.csv(infile, check.names = FALSE)

if (nrow(dat) != EXPECTED_MATCHED_ANALYTIC_N) {
  stop("Analytic N mismatch before EGA. Expected ",
       EXPECTED_MATCHED_ANALYTIC_N, ", found ", nrow(dat))
}

out <- file.path(OUTPUT_DIR, "04_ega_community_check")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

nodes <- c(EQ_CODES, SYMPTOM_CODES)

# Theoretical three-community assignment:
# 1 = Equanimity, 2 = Depression, 3 = Anxiety
theory <- c(rep(1, 4), rep(2, 9), rep(3, 7))
names(theory) <- nodes

# EGAnet's reproducible_seeds() can require a base numeric (double),
# not an R integer. Convert both iteration count and seed explicitly.
ega_iter <- as.numeric(N_EGA_BOOT)
seed_T1_ega  <- as.numeric(SEED) + 10
seed_T1_boot <- as.numeric(SEED) + 11
seed_T2_ega  <- as.numeric(SEED) + 20
seed_T2_boot <- as.numeric(SEED) + 21

safe_dimstab_table <- function(ds, label) {
  if (is.null(ds) || is.null(ds$dimension.stability)) {
    return(data.frame(
      structure = label,
      community = NA,
      structural_consistency = NA,
      average_item_stability = NA
    ))
  }

  sc <- ds$dimension.stability$structural.consistency
  av <- ds$dimension.stability$average.item.stability

  data.frame(
    structure = label,
    community = names(sc),
    structural_consistency = as.numeric(sc),
    average_item_stability = as.numeric(av)
  )
}

run_one <- function(wave) {

  message("\n===== Running EGA for ", wave, " =====")

  x <- dat[, paste0(wave, "_", nodes), drop = FALSE]
  names(x) <- nodes

  ega_seed <- if (wave == "T1") seed_T1_ega else seed_T2_ega
  boot_seed <- if (wave == "T1") seed_T1_boot else seed_T2_boot

  # Empirical EGA
  set.seed(as.integer(ega_seed))
  ega <- EGAnet::EGA(
    data = x,
    corr = "auto",
    model = "glasso",
    algorithm = "walktrap",
    plot.EGA = FALSE,
    verbose = FALSE
  )

  empirical_membership <- as.numeric(ega$wc)
  names(empirical_membership) <- names(ega$wc)

  if (length(empirical_membership) != length(nodes)) {
    stop(wave, ": empirical EGA membership length does not match node count.")
  }

  # Record any node not assigned to an empirical dimension (typically coded 0).
  unassigned <- names(empirical_membership)[
    is.na(empirical_membership) | empirical_membership <= 0
  ]
  if (length(unassigned) > 0) {
    writeLines(
      paste(unassigned, collapse = "\n"),
      file.path(out, paste0(wave, "_EGA_unassigned_nodes.txt"))
    )
    warning(
      wave, ": empirical EGA left these node(s) unassigned: ",
      paste(unassigned, collapse = ", "),
      ". They are retained in the theoretical community solution; do not delete them."
    )
  }

  # Bootstrap EGA
  set.seed(as.integer(boot_seed))
  boot <- EGAnet::bootEGA(
    data = x,
    corr = "auto",
    model = "glasso",
    algorithm = "walktrap",
    iter = ega_iter,
    type = "resampling",
    ncores = as.numeric(N_CORES),
    plot.itemStability = FALSE,
    typicalStructure = FALSE,
    plot.typicalStructure = FALSE,
    seed = boot_seed,
    verbose = FALSE
  )

  # Empirical and theory-referenced item stability
  is_emp <- EGAnet::itemStability(
    boot,
    IS.plot = FALSE
  )

  is_theory <- EGAnet::itemStability(
    boot,
    IS.plot = FALSE,
    structure = theory
  )

  # Dimension stability
  ds_emp <- EGAnet::dimensionStability(
    boot,
    IS.plot = FALSE
  )

  ds_theory <- EGAnet::dimensionStability(
    boot,
    IS.plot = FALSE,
    structure = theory
  )

  empirical_stability <- as.numeric(
    is_emp$item.stability$empirical.dimensions
  )
  theoretical_stability <- as.numeric(
    is_theory$item.stability$empirical.dimensions
  )

  membership <- data.frame(
    node = nodes,
    theoretical_community = as.numeric(theory),
    empirical_community = empirical_membership,
    empirical_item_stability = empirical_stability,
    theoretical_item_stability = theoretical_stability
  )

  # ARI remains computable even if EGAnet encodes an unassigned node as 0.
  ari <- mclust::adjustedRandIndex(
    as.numeric(theory),
    empirical_membership
  )

  # Make bootstrap summary extraction robust to minor version differences.
  bt <- boot$summary.table
  get_bt <- function(candidates) {
    for (nm in candidates) {
      if (!is.null(bt[[nm]])) return(as.numeric(bt[[nm]][1]))
    }
    NA_real_
  }

  summary_row <- data.frame(
    wave = wave,
    n = nrow(x),
    empirical_n_communities = ega$n.dim,
    n_unassigned_nodes = length(unassigned),
    unassigned_nodes = ifelse(length(unassigned) == 0, "", paste(unassigned, collapse = ";")),
    adjusted_rand_vs_theory = ari,
    median_bootstrap_dimensions = get_bt(c("median.dim","Median","median")),
    lower_95_bootstrap_dimensions = get_bt(c("Lower.Quantile","lower","lower.quantile")),
    upper_95_bootstrap_dimensions = get_bt(c("Upper.Quantile","upper","upper.quantile"))
  )

  # Bootstrap dimension-frequency table
  freq <- tryCatch({
    fr <- as.data.frame(boot$frequency)
    if (ncol(fr) >= 2) {
      names(fr)[1:2] <- c("n_dimensions", "proportion")
    } else {
      fr <- data.frame(
        n_dimensions = seq_along(as.numeric(boot$frequency)),
        proportion = as.numeric(boot$frequency)
      )
    }
    fr$wave <- wave
    fr[, c("wave", "n_dimensions", "proportion")]
  }, error = function(e) {
    data.frame(
      wave = wave,
      n_dimensions = NA,
      proportion = NA
    )
  })

  dimstab <- dplyr::bind_rows(
    safe_dimstab_table(ds_emp, "empirical"),
    safe_dimstab_table(ds_theory, "theoretical_3_community")
  )
  dimstab$wave <- wave
  dimstab <- dimstab[, c(
    "wave", "structure", "community",
    "structural_consistency", "average_item_stability"
  )]

  write.csv(
    membership,
    file.path(out, paste0(wave, "_EGA_membership.csv")),
    row.names = FALSE
  )

  write.csv(
    summary_row,
    file.path(out, paste0(wave, "_EGA_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    freq,
    file.path(out, paste0(wave, "_bootstrap_dimension_frequency.csv")),
    row.names = FALSE
  )

  write.csv(
    dimstab,
    file.path(out, paste0(wave, "_dimension_stability.csv")),
    row.names = FALSE
  )

  saveRDS(
    list(
      ega = ega,
      boot = boot,
      item_stability_empirical = is_emp,
      item_stability_theoretical = is_theory,
      dimension_stability_empirical = ds_emp,
      dimension_stability_theoretical = ds_theory
    ),
    file.path(out, paste0(wave, "_EGA_objects.rds"))
  )

  message(
    wave, " done: empirical dimensions = ", ega$n.dim,
    "; unassigned = ",
    ifelse(length(unassigned) == 0, "none", paste(unassigned, collapse = ",")),
    "; ARI vs theory = ", round(ari, 3)
  )

  list(
    membership = membership,
    summary = summary_row,
    frequency = freq,
    dimstab = dimstab
  )
}

r1 <- run_one("T1")
r2 <- run_one("T2")

t1t2_ari <- mclust::adjustedRandIndex(
  r1$membership$empirical_community,
  r2$membership$empirical_community
)

write.csv(
  data.frame(
    metric = "Adjusted Rand Index between T1 and T2 empirical EGA memberships",
    value = t1t2_ari
  ),
  file.path(out, "T1_T2_EGA_similarity.csv"),
  row.names = FALSE
)

message(
  "\nEGA PATCH v2.2 FINISHED.\n",
  "Do not run Phase 2 yet. Zip revision_outputs_v2 and send it back."
)
