# ============================================================
# Data helpers - v2
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

normalize_sid <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.0+$", "", x)
  trimws(x)
}

normalize_phone4 <- function(x) {
  x <- as.character(x)
  x <- gsub("\\.0+$", "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  out <- ifelse(is.na(x), NA_character_, sprintf("%04s", x))
  gsub(" ", "0", out)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

recode_1to4_to_0to3 <- function(x, name = deparse(substitute(x))) {
  x <- safe_numeric(x)
  vals <- sort(unique(x[!is.na(x)]))
  if (!all(vals %in% 1:4)) {
    stop(name, " contains values outside 1-4: ", paste(vals, collapse = ", "))
  }
  x - 1
}

reverse_1to5 <- function(x) {
  x <- safe_numeric(x)
  vals <- sort(unique(x[!is.na(x)]))
  if (!all(vals %in% 1:5)) {
    stop("A 1-5 item contains values outside 1-5: ", paste(vals, collapse = ", "))
  }
  6 - x
}

# Strict response-style rule approved for the revision:
# all 40 analysis items at a wave have exactly the same RAW response category.
# 40 items = PHQ9 + GAD7 + 10 retained EQUA-S items + 14 retained ES items.
strict_straightline_flag <- function(dat, vars) {
  x <- dat[, vars, drop = FALSE]
  apply(x, 1, function(row) {
    z <- suppressWarnings(as.numeric(as.character(row)))
    if (any(is.na(z))) return(FALSE)
    length(unique(z)) == 1L
  })
}

strict_response_category <- function(dat, vars, flag) {
  out <- rep(NA_real_, nrow(dat))
  if (any(flag)) {
    out[flag] <- suppressWarnings(as.numeric(as.character(dat[flag, vars[1]])))
  }
  out
}

raw_wave_strict_vars <- function() {
  c(
    paste0("PHQ", sprintf("%02d", 1:9)),
    paste0("GAD", sprintf("%02d", 1:7)),
    paste0("EQ", sprintf("%02d", c(1,2,5,8,9:14))),
    paste0("ES", sprintf("%02d", c(1:8,11:16)))
  )
}

paired_wave_strict_vars <- function(wave) {
  c(
    paste0(wave, "_", c(
      "ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI",
      "ANX","CTL","WRY","RLX","RST","IRR","AFR"
    )),
    paste0(wave, "_EQ", sprintf("%02d", c(1,2,5,8,9:14))),
    paste0(wave, "_ES", sprintf("%02d", c(1:8,11:16)))
  )
}

read_raw_wave <- function(file, sheet, wave = c("T1", "T2")) {
  wave <- match.arg(wave)
  dat <- readxl::read_excel(file, sheet = sheet, col_types = "text")

  if (ncol(dat) < 56) {
    stop("Raw sheet ", sheet, " has fewer than 56 columns; layout is not as expected.")
  }
  dat <- dat[, 1:56]

  names(dat) <- c(
    "IP", "AGE", "PHONE", "SID", "GEN", "EDU", "ONL", "MEX", "MAJ",
    paste0("PHQ", sprintf("%02d", 1:9)),
    paste0("GAD", sprintf("%02d", 1:7)),
    paste0("EQ", sprintf("%02d", 1:14)),
    paste0("ES", sprintf("%02d", 1:16)),
    "TIMEPOINT"
  )

  # IMPORTANT: strict straightlining is assessed BEFORE any recoding/reversal.
  svars <- raw_wave_strict_vars()
  dat$STRICT_STRAIGHTLINE <- strict_straightline_flag(dat, svars)
  dat$STRICT_RESPONSE_CATEGORY <- strict_response_category(dat, svars, dat$STRICT_STRAIGHTLINE)

  dat <- dat %>%
    mutate(
      SID = normalize_sid(SID),
      PHONE = normalize_phone4(PHONE),
      AGE = safe_numeric(AGE),
      GEN = safe_numeric(GEN),
      MEX = safe_numeric(MEX)
    )

  for (v in c(paste0("PHQ", sprintf("%02d", 1:9)),
              paste0("GAD", sprintf("%02d", 1:7)))) {
    dat[[v]] <- recode_1to4_to_0to3(dat[[v]], paste0(sheet, ":", v))
  }

  for (v in c(paste0("EQ", sprintf("%02d", 1:14)),
              paste0("ES", sprintf("%02d", 1:16)))) {
    dat[[v]] <- safe_numeric(dat[[v]])
  }

  for (v in paste0("EQ", sprintf("%02d", 9:14))) {
    dat[[paste0(v, "R")]] <- reverse_1to5(dat[[v]])
  }
  for (v in paste0("ES", sprintf("%02d", 11:16))) {
    dat[[paste0(v, "R")]] <- reverse_1to5(dat[[v]])
  }

  dat <- dat %>%
    mutate(
      PHQ_TOTAL = rowSums(across(all_of(paste0("PHQ", sprintf("%02d", 1:9)))), na.rm = FALSE),
      GAD_TOTAL = rowSums(across(all_of(paste0("GAD", sprintf("%02d", 1:7)))), na.rm = FALSE),
      EMS = rowMeans(across(all_of(c("EQ01","EQ02","EQ05","EQ08"))), na.rm = FALSE),
      HIN = rowMeans(across(all_of(paste0("EQ", sprintf("%02d", 9:14), "R"))), na.rm = FALSE),
      EAC = rowMeans(across(all_of(paste0("ES", sprintf("%02d", 1:8)))), na.rm = FALSE),
      NRT = rowMeans(across(all_of(paste0("ES", sprintf("%02d", 11:16), "R"))), na.rm = FALSE)
    )

  symptom_map <- c(
    ANH="PHQ01", DEP="PHQ02", SLP="PHQ03", FAT="PHQ04", APP="PHQ05",
    GIL="PHQ06", CON="PHQ07", AGI="PHQ08", SUI="PHQ09",
    ANX="GAD01", CTL="GAD02", WRY="GAD03", RLX="GAD04",
    RST="GAD05", IRR="GAD06", AFR="GAD07"
  )
  for (nm in names(symptom_map)) dat[[nm]] <- dat[[symptom_map[[nm]]]]

  keep_unprefixed <- c("SID","PHONE")
  prefixable <- setdiff(names(dat), keep_unprefixed)
  names(dat)[match(prefixable, names(dat))] <- paste0(wave, "_", prefixable)

  dat
}

read_paired_raw <- function(file) {
  dat <- readxl::read_excel(file, sheet = "Raw_Wide_Selected", col_types = "text")

  # Strict straightlining BEFORE symptom recoding and item reversal.
  for (w in c("T1","T2")) {
    svars <- paired_wave_strict_vars(w)
    missing <- setdiff(svars, names(dat))
    if (length(missing) > 0) {
      stop("Missing raw variables needed for strict straightlining: ",
           paste(missing, collapse=", "))
    }
    flag <- strict_straightline_flag(dat, svars)
    dat[[paste0(w, "_STRICT_STRAIGHTLINE")]] <- flag
    dat[[paste0(w, "_STRICT_RESPONSE_CATEGORY")]] <-
      strict_response_category(dat, svars, flag)
  }

  for (v in c("SID_T1","SID_T2")) dat[[v]] <- normalize_sid(dat[[v]])
  for (v in c("PHONE_T1","PHONE_T2")) dat[[v]] <- normalize_phone4(dat[[v]])
  for (v in c("AGE_T1","AGE_T2","GEN_T1","GEN_T2")) dat[[v]] <- safe_numeric(dat[[v]])

  symptom_codes <- c(
    "ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI",
    "ANX","CTL","WRY","RLX","RST","IRR","AFR"
  )
  for (w in c("T1","T2")) {
    for (s in symptom_codes) {
      v <- paste0(w, "_", s)
      dat[[v]] <- recode_1to4_to_0to3(dat[[v]], v)
    }
  }

  for (w in c("T1","T2")) {
    for (i in 1:14) {
      v <- paste0(w, "_EQ", sprintf("%02d", i))
      dat[[v]] <- safe_numeric(dat[[v]])
    }
    for (i in 1:16) {
      v <- paste0(w, "_ES", sprintf("%02d", i))
      dat[[v]] <- safe_numeric(dat[[v]])
    }

    for (i in 9:14) {
      src <- paste0(w, "_EQ", sprintf("%02d", i))
      dat[[paste0(src, "R")]] <- reverse_1to5(dat[[src]])
    }
    for (i in 11:16) {
      src <- paste0(w, "_ES", sprintf("%02d", i))
      dat[[paste0(src, "R")]] <- reverse_1to5(dat[[src]])
    }

    dat[[paste0(w, "_EMS_calc")]] <-
      rowMeans(dat[, paste0(w, "_EQ", sprintf("%02d", c(1,2,5,8))), drop=FALSE])
    dat[[paste0(w, "_HIN_calc")]] <-
      rowMeans(dat[, paste0(w, "_EQ", sprintf("%02d", 9:14), "R"), drop=FALSE])
    dat[[paste0(w, "_EAC_calc")]] <-
      rowMeans(dat[, paste0(w, "_ES", sprintf("%02d", 1:8)), drop=FALSE])
    dat[[paste0(w, "_NRT_calc")]] <-
      rowMeans(dat[, paste0(w, "_ES", sprintf("%02d", 11:16), "R"), drop=FALSE])

    phq_nodes <- paste0(w, "_", c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI"))
    gad_nodes <- paste0(w, "_", c("ANX","CTL","WRY","RLX","RST","IRR","AFR"))
    dat[[paste0(w, "_PHQ_TOTAL_calc")]] <- rowSums(dat[, phq_nodes, drop=FALSE])
    dat[[paste0(w, "_GAD_TOTAL_calc")]] <- rowSums(dat[, gad_nodes, drop=FALSE])
  }

  for (w in c("T1","T2")) {
    for (d in c("EMS","HIN","EAC","NRT")) {
      legacy <- paste0(w, "_", d)
      if (legacy %in% names(dat)) {
        dat[[paste0(legacy, "_legacy")]] <- safe_numeric(dat[[legacy]])
      }
      dat[[legacy]] <- dat[[paste0(w, "_", d, "_calc")]]
    }
    dat[[paste0(w, "_PHQ_TOTAL")]] <- dat[[paste0(w, "_PHQ_TOTAL_calc")]]
    dat[[paste0(w, "_GAD_TOTAL")]] <- dat[[paste0(w, "_GAD_TOTAL_calc")]]
  }

  dat
}

build_full_t1_eligible_cohort <- function(raw_file, paired_file) {
  t1 <- read_raw_wave(raw_file, "第1次", "T1")
  paired <- read_paired_raw(paired_file)

  if (nrow(t1) != INITIAL_VALID_T1_N) {
    stop("Expected ", INITIAL_VALID_T1_N, " initially valid T1 rows, found ", nrow(t1))
  }
  if (nrow(paired) != MATCHED_PRE_QC_N) {
    stop("Expected ", MATCHED_PRE_QC_N, " pre-QC matched rows, found ", nrow(paired))
  }

  t1 <- t1 %>% mutate(T1_KEY = paste(SID, PHONE, sep="::"))
  paired <- paired %>% mutate(T1_KEY = paste(SID_T1, PHONE_T1, sep="::"))

  if (anyDuplicated(t1$T1_KEY)) stop("Duplicate T1 SID+PHONE keys in initially valid T1 data.")
  if (anyDuplicated(paired$T1_KEY)) stop("Duplicate T1 SID+PHONE keys in paired data.")

  missing_pair_keys <- setdiff(paired$T1_KEY, t1$T1_KEY)
  if (length(missing_pair_keys) > 0) {
    stop("Some paired T1 rows cannot be found in the initially valid T1 cohort.")
  }

  t1_base <- t1 %>% select(-starts_with("T1_IP"), -starts_with("T1_TIMEPOINT"))
  t2_cols <- names(paired)[startsWith(names(paired), "T2_")]

  t2_from_pairs <- paired %>%
    select(T1_KEY, MATCH_TYPE, T1_STRICT_STRAIGHTLINE,
           all_of(t2_cols))

  cohort_all <- t1_base %>%
    left_join(t2_from_pairs, by="T1_KEY", suffix=c("", "_PAIRED")) %>%
    mutate(
      RETAINED_PRE_QC = ifelse(!is.na(MATCH_TYPE), 1L, 0L),
      T2_STRICT_STRAIGHTLINE = dplyr::coalesce(T2_STRICT_STRAIGHTLINE, FALSE),
      T2_INVALID_STRAIGHTLINE = RETAINED_PRE_QC == 1L & T2_STRICT_STRAIGHTLINE,
      RETAINED = ifelse(RETAINED_PRE_QC == 1L & !T2_STRICT_STRAIGHTLINE, 1L, 0L)
    )

  # T1 strict straightliners are not eligible baseline observations.
  cohort <- cohort_all %>% filter(!T1_STRICT_STRAIGHTLINE)

  if (nrow(cohort) != EXPECTED_BASELINE_ELIGIBLE_N) {
    stop("Expected baseline-eligible N=", EXPECTED_BASELINE_ELIGIBLE_N,
         "; found ", nrow(cohort))
  }
  if (sum(cohort$RETAINED) != EXPECTED_MATCHED_ANALYTIC_N) {
    stop("Expected final retained analytic N=", EXPECTED_MATCHED_ANALYTIC_N,
         "; found ", sum(cohort$RETAINED))
  }

  cohort
}

make_descriptives <- function(dat, vars, group = NULL) {
  if (is.null(group)) {
    purrr::map_dfr(vars, function(v) {
      x <- dat[[v]]
      tibble(
        variable = v,
        n = sum(!is.na(x)),
        mean = mean(x, na.rm=TRUE),
        sd = sd(x, na.rm=TRUE),
        min = min(x, na.rm=TRUE),
        max = max(x, na.rm=TRUE)
      )
    })
  } else {
    dat %>%
      group_by(.data[[group]]) %>%
      group_modify(~make_descriptives(.x, vars)) %>%
      ungroup()
  }
}
