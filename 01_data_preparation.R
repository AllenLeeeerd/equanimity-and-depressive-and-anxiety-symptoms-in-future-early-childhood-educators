source("00_config.R")

out <- file.path(OUTPUT_DIR, "01_data_preparation")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

# Reconstruct initially valid waves and pre-QC matched data
t1_all <- read_raw_wave(RAW_FILE, "第1次", "T1")
t2_all <- read_raw_wave(RAW_FILE, "第2次", "T2")
paired_pre <- read_paired_raw(PAIRED_FILE)

if (nrow(t1_all) != INITIAL_VALID_T1_N) stop("T1 initially valid N mismatch.")
if (nrow(t2_all) != INITIAL_VALID_T2_N) stop("T2 initially valid N mismatch.")
if (nrow(paired_pre) != MATCHED_PRE_QC_N) stop("Pre-QC matched N mismatch.")

# Confirm strict straightlining counts in all initially valid wave data
n_t1_strict_all <- sum(t1_all$T1_STRICT_STRAIGHTLINE)
n_t2_strict_all <- sum(t2_all$T2_STRICT_STRAIGHTLINE)
if (n_t1_strict_all != EXPECTED_T1_STRICT_ALL) {
  stop("Expected ", EXPECTED_T1_STRICT_ALL, " T1 strict straightliners; found ", n_t1_strict_all)
}
if (n_t2_strict_all != EXPECTED_T2_STRICT_ALL) {
  stop("Expected ", EXPECTED_T2_STRICT_ALL, " T2 strict straightliners; found ", n_t2_strict_all)
}

# Confirm strict straightlining in the 392 pre-QC matched sample
n_t1_strict_matched <- sum(paired_pre$T1_STRICT_STRAIGHTLINE)
n_t2_strict_matched <- sum(paired_pre$T2_STRICT_STRAIGHTLINE)
n_union_matched <- sum(paired_pre$T1_STRICT_STRAIGHTLINE | paired_pre$T2_STRICT_STRAIGHTLINE)
n_both_matched <- sum(paired_pre$T1_STRICT_STRAIGHTLINE & paired_pre$T2_STRICT_STRAIGHTLINE)

if (n_t1_strict_matched != EXPECTED_T1_STRICT_MATCHED) stop("Matched T1 strict count mismatch.")
if (n_t2_strict_matched != EXPECTED_T2_STRICT_MATCHED) stop("Matched T2 strict count mismatch.")
if (n_union_matched != EXPECTED_MATCHED_STRICT_UNION) stop("Matched strict union count mismatch.")

# Final analytic matched sample: exclude if strict at either wave
paired <- paired_pre %>%
  filter(!T1_STRICT_STRAIGHTLINE & !T2_STRICT_STRAIGHTLINE)

if (nrow(paired) != EXPECTED_MATCHED_ANALYTIC_N) {
  stop("Expected final matched analytic N=", EXPECTED_MATCHED_ANALYTIC_N,
       "; found ", nrow(paired))
}

exact <- paired %>% filter(MATCH_TYPE == "pair")
if (nrow(exact) != EXPECTED_EXACT_ANALYTIC_N) {
  stop("Expected exact-match analytic N=", EXPECTED_EXACT_ANALYTIC_N,
       "; found ", nrow(exact))
}

# Full baseline-eligible cohort for attrition/FIML: T1 strict cases removed.
# Matched T2 strict cases are treated as having no valid T2 data.
cohort <- build_full_t1_eligible_cohort(RAW_FILE, PAIRED_FILE)

# Save analysis-ready datasets
write.csv(paired_pre, file.path(out, "matched_preQC_392.csv"),
          row.names=FALSE, fileEncoding="UTF-8")
write.csv(paired, file.path(out, "matched_analytic_clean.csv"),
          row.names=FALSE, fileEncoding="UTF-8")
write.csv(exact, file.path(out, "exact_match_analytic_clean.csv"),
          row.names=FALSE, fileEncoding="UTF-8")
write.csv(cohort, file.path(out, "full_T1_eligible_with_T2_missing.csv"),
          row.names=FALSE, fileEncoding="UTF-8")

# Strict straightlining audit files
strict_t1_all <- t1_all %>%
  filter(T1_STRICT_STRAIGHTLINE) %>%
  transmute(
    wave="T1", SID=SID, PHONE=PHONE,
    raw_response_category=T1_STRICT_RESPONSE_CATEGORY
  )

strict_t2_all <- t2_all %>%
  filter(T2_STRICT_STRAIGHTLINE) %>%
  transmute(
    wave="T2", SID=SID, PHONE=PHONE,
    raw_response_category=T2_STRICT_RESPONSE_CATEGORY
  )

strict_matched <- paired_pre %>%
  filter(T1_STRICT_STRAIGHTLINE | T2_STRICT_STRAIGHTLINE) %>%
  transmute(
    ID, MATCH_TYPE, SID_T1, PHONE_T1, SID_T2, PHONE_T2,
    T1_STRICT_STRAIGHTLINE, T2_STRICT_STRAIGHTLINE,
    T1_STRICT_RESPONSE_CATEGORY, T2_STRICT_RESPONSE_CATEGORY
  )

write.csv(bind_rows(strict_t1_all, strict_t2_all),
          file.path(out, "strict_straightliners_all_initially_valid_waves.csv"),
          row.names=FALSE)
write.csv(strict_matched,
          file.path(out, "strict_straightliners_matched_excluded.csv"),
          row.names=FALSE)

qc_summary <- data.frame(
  metric=c(
    "T1 initially valid", "T1 strict straightliners", "T1 baseline eligible",
    "T2 initially valid", "T2 strict straightliners",
    "Matched pre-QC", "Matched T1 strict", "Matched T2 strict",
    "Matched strict at both waves", "Matched excluded union",
    "Final matched analytic", "Exact-match pre-QC", "Final exact-match analytic"
  ),
  n=c(
    nrow(t1_all), n_t1_strict_all, nrow(cohort),
    nrow(t2_all), n_t2_strict_all,
    nrow(paired_pre), n_t1_strict_matched, n_t2_strict_matched,
    n_both_matched, n_union_matched,
    nrow(paired), EXACT_MATCH_PRE_QC_N, nrow(exact)
  )
)
write.csv(qc_summary, file.path(out, "response_style_qc_summary.csv"), row.names=FALSE)

# Sample-flow table
sample_flow <- tibble::tribble(
  ~stage, ~starting_n, ~retained_n, ~excluded_n, ~note,
  "T1 initial validity screening", RAW_T1_N, INITIAL_VALID_T1_N,
  RAW_T1_N - INITIAL_VALID_T1_N,
  "Excluded if any of three instructed-response attention checks failed or completion time < 5 minutes.",
  "T1 post-hoc response-style QC", INITIAL_VALID_T1_N,
  INITIAL_VALID_T1_N - n_t1_strict_all, n_t1_strict_all,
  "Excluded only under the pre-specified revision rule: identical raw response category across all 40 focal analysis items.",
  "T2 initial validity screening", RAW_T2_N, INITIAL_VALID_T2_N,
  RAW_T2_N - INITIAL_VALID_T2_N,
  "Excluded if any of three instructed-response attention checks failed or completion time < 5 minutes.",
  "T2 post-hoc response-style QC", INITIAL_VALID_T2_N,
  INITIAL_VALID_T2_N - n_t2_strict_all, n_t2_strict_all,
  "Same strict cross-instrument straightlining rule.",
  "Longitudinal matching before response-style QC", INITIAL_VALID_T1_N,
  MATCHED_PRE_QC_N, INITIAL_VALID_T1_N - MATCHED_PRE_QC_N,
  "Pre-QC matched sample based on the original matching procedure.",
  "Final longitudinal analytic sample", MATCHED_PRE_QC_N,
  nrow(paired), n_union_matched,
  "Excluded matched participants who met the strict straightlining rule at either wave."
)
write.csv(sample_flow, file.path(out, "sample_flow.csv"), row.names=FALSE)

# Audit recomputed dimension scores against legacy paired file
audit <- list()
for (w in c("T1","T2")) {
  for (d in c("EMS","HIN","EAC","NRT")) {
    legacy <- paste0(w,"_",d,"_legacy")
    current <- paste0(w,"_",d)
    if (legacy %in% names(paired)) {
      diff <- paired[[current]] - paired[[legacy]]
      audit[[paste0(w,"_",d)]] <- data.frame(
        variable=paste0(w,"_",d),
        max_abs_difference=max(abs(diff), na.rm=TRUE),
        mean_difference=mean(diff, na.rm=TRUE)
      )
    }
  }
}
write.csv(bind_rows(audit), file.path(out, "score_audit_vs_legacy.csv"), row.names=FALSE)

# Range checks
range_checks <- tibble(
  variable_group=c(
    "PHQ symptoms T1","PHQ symptoms T2",
    "GAD symptoms T1","GAD symptoms T2",
    "Equanimity dimensions T1","Equanimity dimensions T2"
  ),
  min_value=c(
    min(as.matrix(paired[,paste0("T1_",c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI"))]),na.rm=TRUE),
    min(as.matrix(paired[,paste0("T2_",c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI"))]),na.rm=TRUE),
    min(as.matrix(paired[,paste0("T1_",c("ANX","CTL","WRY","RLX","RST","IRR","AFR"))]),na.rm=TRUE),
    min(as.matrix(paired[,paste0("T2_",c("ANX","CTL","WRY","RLX","RST","IRR","AFR"))]),na.rm=TRUE),
    min(as.matrix(paired[,paste0("T1_",EQ_CODES)]),na.rm=TRUE),
    min(as.matrix(paired[,paste0("T2_",EQ_CODES)]),na.rm=TRUE)
  ),
  max_value=c(
    max(as.matrix(paired[,paste0("T1_",c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI"))]),na.rm=TRUE),
    max(as.matrix(paired[,paste0("T2_",c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI"))]),na.rm=TRUE),
    max(as.matrix(paired[,paste0("T1_",c("ANX","CTL","WRY","RLX","RST","IRR","AFR"))]),na.rm=TRUE),
    max(as.matrix(paired[,paste0("T2_",c("ANX","CTL","WRY","RLX","RST","IRR","AFR"))]),na.rm=TRUE),
    max(as.matrix(paired[,paste0("T1_",EQ_CODES)]),na.rm=TRUE),
    max(as.matrix(paired[,paste0("T2_",EQ_CODES)]),na.rm=TRUE)
  )
)
write.csv(range_checks, file.path(out, "range_checks.csv"), row.names=FALSE)

# Reliability audit on the final analytic matched sample
reliability_one <- function(dat, vars, scale, wave) {
  X <- as.data.frame(dat[,vars,drop=FALSE])
  a <- tryCatch(psych::alpha(X, warnings=FALSE, check.keys=FALSE)$total$raw_alpha,
                error=function(e) NA_real_)
  o <- tryCatch(psych::omega(X, nfactors=1, plot=FALSE, warnings=FALSE)$omega.tot,
                error=function(e) NA_real_)
  data.frame(scale=scale, wave=wave, n_items=length(vars), alpha=a, omega_total=o)
}

rel <- list()
ctr <- 1L
for (w in c("T1","T2")) {
  scale_vars <- list(
    PHQ9=paste0(w,"_",c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI")),
    GAD7=paste0(w,"_",c("ANX","CTL","WRY","RLX","RST","IRR","AFR")),
    EMS=paste0(w,"_EQ",sprintf("%02d",c(1,2,5,8))),
    HIN=paste0(w,"_EQ",sprintf("%02d",9:14),"R"),
    EAC=paste0(w,"_ES",sprintf("%02d",1:8)),
    NRT=paste0(w,"_ES",sprintf("%02d",11:16),"R")
  )
  for (sc in names(scale_vars)) {
    rel[[ctr]] <- reliability_one(paired, scale_vars[[sc]], sc, w)
    ctr <- ctr + 1L
  }
}
write.csv(bind_rows(rel), file.path(out, "reliability_alpha_omega.csv"), row.names=FALSE)

message("01_data_preparation v2 complete. Review: ", normalizePath(out))
