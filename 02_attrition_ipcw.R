source("00_config.R")

infile <- file.path(OUTPUT_DIR, "01_data_preparation", "full_T1_eligible_with_T2_missing.csv")
if (!file.exists(infile)) stop("Run 01_data_preparation.R first.")
cohort <- read.csv(infile, check.names=FALSE)

if (nrow(cohort) != EXPECTED_BASELINE_ELIGIBLE_N) stop("Baseline-eligible N mismatch.")
if (sum(cohort$RETAINED) != EXPECTED_MATCHED_ANALYTIC_N) stop("Final retained N mismatch.")

out <- file.path(OUTPUT_DIR, "02_attrition_ipcw")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

vars_cont <- c("T1_AGE","T1_PHQ_TOTAL","T1_GAD_TOTAL","T1_EAC","T1_NRT","T1_EMS","T1_HIN")

# Group descriptives
desc <- cohort %>%
  mutate(RETAINED_LABEL=ifelse(RETAINED==1,"Retained","Not retained")) %>%
  group_by(RETAINED_LABEL) %>%
  group_modify(~make_descriptives(.x, vars_cont)) %>%
  ungroup()
write.csv(desc, file.path(out, "retention_group_descriptives.csv"), row.names=FALSE)

extract_es <- function(es_obj) {
  if (is.null(es_obj)) return(c(value=NA, low=NA, high=NA))
  z <- as.data.frame(es_obj)
  value_name <- grep("Hedges", names(z), value=TRUE)[1]
  if (is.na(value_name)) value_name <- names(z)[1]
  c(
    value=as.numeric(z[[value_name]][1]),
    low=as.numeric(z$CI_low[1]),
    high=as.numeric(z$CI_high[1])
  )
}

# Welch t tests + Hedges g. Direction is explicitly retained minus non-retained.
cont_tests <- lapply(vars_cont, function(v) {
  x1 <- cohort[[v]][cohort$RETAINED==1]
  x0 <- cohort[[v]][cohort$RETAINED==0]
  tt <- t.test(x1, x0, var.equal=FALSE)

  eg <- tryCatch(
    effectsize::hedges_g(x1, x0, paired=FALSE, ci=.95, alternative="two.sided"),
    error=function(e) NULL
  )
  es <- extract_es(eg)

  data.frame(
    variable=v,
    retained_n=sum(!is.na(x1)),
    dropout_n=sum(!is.na(x0)),
    retained_mean=mean(x1,na.rm=TRUE),
    retained_sd=sd(x1,na.rm=TRUE),
    dropout_mean=mean(x0,na.rm=TRUE),
    dropout_sd=sd(x0,na.rm=TRUE),
    mean_difference_retained_minus_dropout=mean(x1,na.rm=TRUE)-mean(x0,na.rm=TRUE),
    mean_diff_ci_low=tt$conf.int[1],
    mean_diff_ci_high=tt$conf.int[2],
    t=unname(tt$statistic),
    df=unname(tt$parameter),
    p=tt$p.value,
    hedges_g_retained_minus_dropout=es["value"],
    hedges_g_ci_low=es["low"],
    hedges_g_ci_high=es["high"]
  )
}) %>% bind_rows()
write.csv(cont_tests, file.path(out, "attrition_continuous_tests.csv"), row.names=FALSE)

# Gender attrition
gender_tab <- table(Gender=cohort$T1_GEN, Retained=cohort$RETAINED)
chi <- chisq.test(gender_tab, correct=FALSE)

cv <- tryCatch(
  effectsize::cramers_v(gender_tab, adjust=TRUE, ci=.95, alternative="two.sided"),
  error=function(e) NULL
)

extract_cv <- function(obj) {
  if (is.null(obj)) return(c(value=NA, low=NA, high=NA))
  z <- as.data.frame(obj)
  valname <- grep("^Cramers_v", names(z), value=TRUE)[1]
  if (is.na(valname)) valname <- names(z)[1]
  c(value=as.numeric(z[[valname]][1]),
    low=as.numeric(z$CI_low[1]),
    high=as.numeric(z$CI_high[1]))
}
cvv <- extract_cv(cv)

gender_test <- data.frame(
  chi_square=unname(chi$statistic),
  df=unname(chi$parameter),
  p=chi$p.value,
  cramers_v=cvv["value"],
  cramers_v_ci_low=cvv["low"],
  cramers_v_ci_high=cvv["high"]
)
write.csv(as.data.frame.matrix(gender_tab), file.path(out, "gender_by_retention_counts.csv"))
write.csv(gender_test, file.path(out, "gender_attrition_test.csv"), row.names=FALSE)

# Explicit gender-specific retention and attrition rates
gender_rates <- cohort %>%
  mutate(Gender=ifelse(T1_GEN==1,"Male",ifelse(T1_GEN==2,"Female",as.character(T1_GEN)))) %>%
  group_by(Gender) %>%
  summarise(
    n_total=n(),
    n_retained=sum(RETAINED==1),
    n_not_retained=sum(RETAINED==0),
    retention_rate=n_retained/n_total,
    attrition_rate=n_not_retained/n_total,
    .groups="drop"
  )
write.csv(gender_rates, file.path(out, "gender_specific_retention_attrition_rates.csv"), row.names=FALSE)

# Logistic retention model based only on observed baseline characteristics
dat <- cohort %>%
  mutate(
    AGE_Z=as.numeric(scale(T1_AGE)),
    PHQ_Z=as.numeric(scale(T1_PHQ_TOTAL)),
    GAD_Z=as.numeric(scale(T1_GAD_TOTAL)),
    EAC_Z=as.numeric(scale(T1_EAC)),
    NRT_Z=as.numeric(scale(T1_NRT)),
    EMS_Z=as.numeric(scale(T1_EMS)),
    HIN_Z=as.numeric(scale(T1_HIN)),
    GEN01=ifelse(T1_GEN==2,1,ifelse(T1_GEN==1,0,NA))
  )

ret_fit <- glm(
  RETAINED ~ AGE_Z + GEN01 + PHQ_Z + GAD_Z + EAC_Z + NRT_Z + EMS_Z + HIN_Z,
  data=dat, family=binomial()
)

coef_tab <- summary(ret_fit)$coefficients %>%
  as.data.frame() %>%
  tibble::rownames_to_column("predictor") %>%
  rename(estimate=Estimate, se=`Std. Error`, z=`z value`, p=`Pr(>|z|)`) %>%
  mutate(
    odds_ratio=exp(estimate),
    or_ci_low=exp(estimate - 1.96*se),
    or_ci_high=exp(estimate + 1.96*se)
  )
write.csv(coef_tab, file.path(out, "retention_logistic_model.csv"), row.names=FALSE)
saveRDS(ret_fit, file.path(out, "retention_logistic_model.rds"))

# Stabilized inverse probability of censoring weights (IPCW)
dat$P_RETAIN <- predict(ret_fit, type="response")
p_marginal <- mean(dat$RETAINED==1)
dat$IPCW_RAW <- ifelse(dat$RETAINED==1, 1/dat$P_RETAIN, NA)
dat$IPCW_STAB <- ifelse(dat$RETAINED==1, p_marginal/dat$P_RETAIN, NA)

q <- quantile(dat$IPCW_STAB[dat$RETAINED==1], IPCW_TRIM, na.rm=TRUE)
dat$IPCW_TRIM <- dat$IPCW_STAB
dat$IPCW_TRIM[dat$RETAINED==1] <-
  pmin(pmax(dat$IPCW_STAB[dat$RETAINED==1], q[1]), q[2])

w <- dat$IPCW_TRIM[dat$RETAINED==1]
ess <- sum(w)^2 / sum(w^2)

weight_diag <- data.frame(
  baseline_eligible_n=nrow(dat),
  n_retained=sum(dat$RETAINED==1),
  marginal_retention_probability=p_marginal,
  predicted_retention_min=min(dat$P_RETAIN,na.rm=TRUE),
  predicted_retention_p50=median(dat$P_RETAIN,na.rm=TRUE),
  predicted_retention_max=max(dat$P_RETAIN,na.rm=TRUE),
  stabilized_weight_min=min(dat$IPCW_STAB[dat$RETAINED==1]),
  stabilized_weight_p50=median(dat$IPCW_STAB[dat$RETAINED==1]),
  stabilized_weight_max=max(dat$IPCW_STAB[dat$RETAINED==1]),
  trim_low=q[1],
  trim_high=q[2],
  trimmed_weight_min=min(w),
  trimmed_weight_max=max(w),
  effective_sample_size=ess
)
write.csv(weight_diag, file.path(out, "ipcw_diagnostics.csv"), row.names=FALSE)

# Attach weights to final analytic paired sample
paired_file <- file.path(OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv")
paired <- read.csv(paired_file, check.names=FALSE)
paired$SID_T1 <- normalize_sid(paired$SID_T1)
paired$PHONE_T1 <- normalize_phone4(paired$PHONE_T1)

cohort_weights <- dat %>%
  filter(RETAINED==1) %>%
  select(T1_KEY, P_RETAIN, IPCW_STAB, IPCW_TRIM)

paired <- paired %>%
  mutate(T1_KEY=paste(SID_T1, PHONE_T1, sep="::")) %>%
  left_join(cohort_weights, by="T1_KEY")

if (any(is.na(paired$IPCW_TRIM))) stop("Missing IPCW after joining to final paired sample.")
write.csv(paired, file.path(out, "matched_analytic_with_ipcw.csv"), row.names=FALSE)

p <- ggplot(paired, aes(x=IPCW_TRIM)) +
  geom_histogram(bins=30) +
  theme_minimal() +
  labs(x="Stabilized, truncated IPCW", y="Count")
ggsave(file.path(out, "ipcw_distribution.pdf"), p, width=7, height=5)

message("02_attrition_ipcw v2 complete. Review: ", normalizePath(out))
