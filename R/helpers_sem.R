# ============================================================
# SEM / longitudinal measurement-invariance helpers - v2
# ============================================================

suppressPackageStartupMessages({
  library(lavaan)
  library(semTools)
  library(dplyr)
  library(tidyr)
})

safe_fit_measure <- function(fm, candidates) {
  for (nm in candidates) {
    if (nm %in% names(fm) && is.finite(fm[[nm]])) return(unname(fm[[nm]]))
  }
  NA_real_
}

fit_indices <- function(fit, model_name) {
  fm <- tryCatch(lavaan::fitMeasures(fit), error=function(e) numeric())

  data.frame(
    model=model_name,
    chisq=safe_fit_measure(fm, c("chisq.scaled","chisq")),
    df=safe_fit_measure(fm, c("df.scaled","df")),
    pvalue=safe_fit_measure(fm, c("pvalue.scaled","pvalue")),
    cfi=safe_fit_measure(fm, c("cfi.robust","cfi.scaled","cfi")),
    tli=safe_fit_measure(fm, c("tli.robust","tli.scaled","tli")),
    rmsea=safe_fit_measure(fm, c("rmsea.robust","rmsea.scaled","rmsea")),
    rmsea_ci_low=safe_fit_measure(fm, c("rmsea.ci.lower.robust","rmsea.ci.lower.scaled","rmsea.ci.lower")),
    rmsea_ci_high=safe_fit_measure(fm, c("rmsea.ci.upper.robust","rmsea.ci.upper.scaled","rmsea.ci.upper")),
    srmr=safe_fit_measure(fm, c("srmr")),
    stringsAsFactors=FALSE
  )
}

make_long_ind_names <- function(factor_items) {
  out <- list()
  for (fac in names(factor_items)) {
    items <- factor_items[[fac]]
    for (i in seq_along(items)) {
      nm <- paste0(fac, "_item", i)
      out[[nm]] <- c(
        setNames(paste0("T1_", items[i]), paste0("T1_", fac)),
        setNames(paste0("T2_", items[i]), paste0("T2_", fac))
      )
    }
  }
  out
}

build_configural_cfa <- function(factor_items) {
  lines <- character()
  for (fac in names(factor_items)) {
    items <- factor_items[[fac]]
    lines <- c(
      lines,
      paste0("T1_", fac, " =~ ", paste0("T1_", items, collapse=" + ")),
      paste0("T2_", fac, " =~ ", paste0("T2_", items, collapse=" + "))
    )
  }
  paste(lines, collapse="\n")
}

category_structure_table <- function(dat, base_items) {
  bind_rows(lapply(base_items, function(item) {
    x1 <- dat[[paste0("T1_",item)]]
    x2 <- dat[[paste0("T2_",item)]]
    data.frame(
      item=item,
      T1_categories=paste(sort(unique(x1[!is.na(x1)])), collapse="|"),
      T2_categories=paste(sort(unique(x2[!is.na(x2)])), collapse="|"),
      T1_n_categories=length(unique(x1[!is.na(x1)])),
      T2_n_categories=length(unique(x2[!is.na(x2)])),
      same_number_categories=
        length(unique(x1[!is.na(x1)])) == length(unique(x2[!is.na(x2)]))
    )
  }))
}

run_longitudinal_invariance <- function(dat, factor_items, scale_name, out_dir) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  config_model <- build_configural_cfa(factor_items)
  long_fac_names <- lapply(
    names(factor_items),
    function(fac) c(paste0("T1_",fac), paste0("T2_",fac))
  )
  names(long_fac_names) <- names(factor_items)
  long_ind_names <- make_long_ind_names(factor_items)

  base_items <- unique(unlist(factor_items))
  ordered_vars <- unique(unlist(lapply(base_items, function(item) {
    c(paste0("T1_",item), paste0("T2_",item))
  })))

  cat_tab <- category_structure_table(dat, base_items)
  write.csv(cat_tab,
            file.path(out_dir, paste0(scale_name, "_category_structure.csv")),
            row.names=FALSE)

  if (any(!cat_tab$same_number_categories)) {
    bad <- cat_tab$item[!cat_tab$same_number_categories]
    stop(
      scale_name,
      ": longitudinal ordered indicators do not have the same number of endorsed categories: ",
      paste(bad, collapse=", "),
      ". Review the category-harmonization step before fitting invariance."
    )
  }

  subdat <- dat[, ordered_vars, drop=FALSE]
  subdat[] <- lapply(subdat, function(x) ordered(as.integer(x)))

  make_syn <- function(long_equal=NULL) {
    args <- list(
      configural.model=config_model,
      longFacNames=long_fac_names,
      longIndNames=long_ind_names,
      ID.fac="std.lv",
      ID.cat="Wu.Estabrook.2016",
      auto=1L,
      data=subdat,
      ordered=ordered_vars,
      parameterization="theta"
    )
    if (!is.null(long_equal)) args$long.equal <- long_equal
    do.call(semTools::measEq.syntax, args)
  }

  fit_syn <- function(syn) {
    lavaan::cfa(
      as.character(syn),
      data=subdat,
      ordered=ordered_vars,
      estimator="WLSMV",
      parameterization="theta"
    )
  }

  # Recommended Wu-Estabrook sequence for ordered longitudinal indicators:
  # configural -> thresholds -> thresholds+loadings -> thresholds+loadings+intercepts
  syn_config <- make_syn(NULL)
  fit_config <- fit_syn(syn_config)

  syn_threshold <- make_syn(c("thresholds"))
  fit_threshold <- fit_syn(syn_threshold)

  syn_metric <- make_syn(c("thresholds","loadings"))
  fit_metric <- fit_syn(syn_metric)

  syn_scalar <- make_syn(c("thresholds","loadings","intercepts"))
  fit_scalar <- fit_syn(syn_scalar)

  fits <- bind_rows(
    fit_indices(fit_config, "configural"),
    fit_indices(fit_threshold, "threshold_invariance"),
    fit_indices(fit_metric, "threshold_plus_loading_invariance"),
    fit_indices(fit_scalar, "threshold_loading_intercept_invariance")
  )

  fits$delta_cfi_from_previous <- c(NA, diff(fits$cfi))
  fits$delta_rmsea_from_previous <- c(NA, diff(fits$rmsea))
  fits$delta_srmr_from_previous <- c(NA, diff(fits$srmr))

  write.csv(fits,
            file.path(out_dir, paste0(scale_name, "_invariance_fit.csv")),
            row.names=FALSE)

  safe_lrt <- function(a,b,label) {
    z <- tryCatch(
      as.data.frame(lavaan::lavTestLRT(a,b)),
      error=function(e) data.frame(error=e$message)
    )
    write.csv(z, file.path(out_dir, paste0(scale_name, "_", label, "_LRT.csv")),
              row.names=FALSE)
  }

  safe_lrt(fit_config, fit_threshold, "config_vs_threshold")
  safe_lrt(fit_threshold, fit_metric, "threshold_vs_metric")
  safe_lrt(fit_metric, fit_scalar, "metric_vs_scalar")

  saveRDS(
    list(
      configural=fit_config,
      threshold=fit_threshold,
      metric=fit_metric,
      scalar=fit_scalar,
      syntax=list(
        configural=syn_config,
        threshold=syn_threshold,
        metric=syn_metric,
        scalar=syn_scalar
      )
    ),
    file.path(out_dir, paste0(scale_name, "_invariance_models.rds"))
  )

  invisible(list(
    configural=fit_config,
    threshold=fit_threshold,
    metric=fit_metric,
    scalar=fit_scalar,
    fit_table=fits,
    category_table=cat_tab
  ))
}
