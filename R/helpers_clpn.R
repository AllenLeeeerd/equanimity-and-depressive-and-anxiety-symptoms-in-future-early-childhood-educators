# ============================================================
# CLPN helpers
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmnet)
  library(bootnet)
  library(qgraph)
  library(networktools)
  library(ggplot2)
})

SYMPTOM_CODES <- c(
  "ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI",
  "ANX","CTL","WRY","RLX","RST","IRR","AFR"
)

SYMPTOM_LABELS <- c(
  "Anhedonia","Depressed mood","Sleep problems","Fatigue","Appetite change",
  "Guilt / worthlessness","Concentration problems","Psychomotor agitation / retardation",
  "Suicidal ideation","Nervous / anxious","Uncontrollable worry","Excessive worry",
  "Trouble relaxing","Restlessness","Irritability","Fear something awful may happen"
)

EQ_CODES <- c("EAC","NRT","EMS","HIN")
EQ_LABELS <- c(
  EAC="Experiential acceptance",
  NRT="Non-reactivity",
  EMS="Even-minded state of mind",
  HIN="Hedonic independence"
)

z_safe <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) stop("Cannot standardize a constant variable.")
  (x - mean(x, na.rm = TRUE)) / s
}

make_foldid <- function(n, k = N_FOLDS, seed = SEED) {
  set.seed(seed + n)
  sample(rep(seq_len(k), length.out = n))
}

# Display domains remain Equanimity / Depression / Anxiety so readers can see
# the PHQ/GAD distinction in figures and edge tables.
make_display_communities <- function(node_codes) {
  ifelse(
    node_codes %in% EQ_CODES, "Equanimity",
    ifelse(node_codes %in% c("ANH","DEP","SLP","FAT","APP","GIL","CON","AGI","SUI","PHQ_TOTAL"),
           "Depression",
           ifelse(node_codes %in% c("ANX","CTL","WRY","RLX","RST","IRR","AFR","GAD_TOTAL"),
                  "Anxiety", "Other"))
  )
}

# Primary bridge metrics use two communities:
# Equanimity vs Distress.
# This choice is based on the 500-bootstrap EGA:
# - the four equanimity dimensions formed one empirical community at T2;
# - at T1 EAC/EMS/HIN clustered together and NRT was unassigned rather than
#   forming a second equanimity system;
# - the depression/anxiety boundary was less stable across waves.
make_bridge_communities <- function(node_codes) {
  ifelse(node_codes %in% EQ_CODES, "Equanimity", "Distress")
}

# Backward-compatible alias used by bootnet bridge statistics.
make_communities <- make_bridge_communities

make_group_index <- function(node_codes) {
  comm <- make_display_communities(node_codes)
  split(
    seq_along(node_codes),
    factor(comm, levels=c("Equanimity","Depression","Anxiety"))
  )
}

prepare_clpn_frame <- function(dat, node_codes, covariates = PRIMARY_COVARIATES,
                               weight_col = NULL) {
  t1 <- paste0("T1_", node_codes)
  t2 <- paste0("T2_", node_codes)
  needed <- c(t1, t2, covariates, weight_col)
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) stop("Missing CLPN variables: ", paste(missing, collapse=", "))

  out <- dat %>% select(all_of(needed))

  # Network variables are standardized explicitly before estimation so the
  # stored coefficients are standardized regression coefficients.
  for (v in c(t1, t2)) out[[v]] <- z_safe(out[[v]])

  # Age is standardized; gender remains a 0/1 adjustment variable.
  if ("AGE_T1" %in% names(out)) out$AGE_T1 <- z_safe(out$AGE_T1)

  if ("GEN_T1" %in% names(out)) {
    g <- as.numeric(out$GEN_T1)
    vals <- sort(unique(g[!is.na(g)]))
    if (all(vals %in% c(1,2))) {
      # 0 = male, 1 = female
      out$GEN_T1 <- ifelse(g == 2, 1, 0)
    } else if (all(vals %in% c(0,1))) {
      out$GEN_T1 <- g
    } else {
      stop("GEN_T1 coding is not recognized: ", paste(vals, collapse=", "))
    }
  }

  out
}

build_clpn_estimator <- function(node_codes, covariates = PRIMARY_COVARIATES,
                                 weight_col = NULL,
                                 lambda_rule = LAMBDA_RULE,
                                 nfolds = N_FOLDS,
                                 seed = SEED) {
  k <- length(node_codes)
  num_cov <- length(covariates)
  has_weight <- !is.null(weight_col)

  force(k); force(num_cov); force(node_codes); force(covariates)
  force(weight_col); force(lambda_rule); force(nfolds); force(seed)

  function(df) {
    df <- as.data.frame(df)

    # Input order: T1 nodes, T2 nodes, covariates, optional weight
    t1_idx <- seq_len(k)
    t2_idx <- k + seq_len(k)
    cov_idx <- if (num_cov > 0) (2*k + 1):(2*k + num_cov) else integer(0)
    wt_idx <- if (has_weight) ncol(df) else integer(0)

    X <- as.matrix(df[, c(t1_idx, cov_idx), drop=FALSE])
    # Network predictors are penalized; covariates are adjustment variables and
    # are kept unpenalized.
    penalty_factor <- c(rep(1, k), rep(0, num_cov))

    if (has_weight) {
      w <- as.numeric(df[[wt_idx]])
      if (any(!is.finite(w)) || any(w <= 0)) stop("Invalid IPCW values.")
    } else {
      w <- rep(1, nrow(df))
    }

    foldid <- make_foldid(nrow(df), nfolds, seed)
    graph <- matrix(0, nrow=k, ncol=k)
    lambda_tbl <- data.frame(
      outcome=node_codes, lambda_min=NA_real_, lambda_1se=NA_real_,
      selected_lambda=NA_real_, n_nonzero=NA_integer_
    )

    for (j in seq_len(k)) {
      y <- as.numeric(df[[t2_idx[j]]])

      fit <- glmnet::cv.glmnet(
        x = X,
        y = y,
        weights = w,
        family = "gaussian",
        alpha = 1,
        standardize = FALSE,
        intercept = TRUE,
        penalty.factor = penalty_factor,
        nfolds = nfolds,
        foldid = foldid,
        type.measure = "mse"
      )

      beta <- as.numeric(coef(fit, s=lambda_rule))[-1]
      graph[, j] <- beta[seq_len(k)]

      lambda_tbl$lambda_min[j] <- fit$lambda.min
      lambda_tbl$lambda_1se[j] <- fit$lambda.1se
      lambda_tbl$selected_lambda[j] <- if (lambda_rule == "lambda.min") fit$lambda.min else fit$lambda.1se
      lambda_tbl$n_nonzero[j] <- sum(abs(beta[seq_len(k)]) > 0)
    }

    rownames(graph) <- node_codes
    colnames(graph) <- node_codes

    list(
      graph = graph,
      results = list(lambda = lambda_tbl)
    )
  }
}

estimate_clpn <- function(dat, node_codes, out_dir,
                          covariates = PRIMARY_COVARIATES,
                          weight_col = NULL,
                          suffix = "main") {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  frame <- prepare_clpn_frame(dat, node_codes, covariates, weight_col)
  label_map <- c(
    EQ_LABELS,
    setNames(SYMPTOM_LABELS, SYMPTOM_CODES),
    PHQ_TOTAL="PHQ-9 total",
    GAD_TOTAL="GAD-7 total"
  )
  labels <- unname(label_map[node_codes])
  communities <- make_bridge_communities(node_codes)
  display_communities <- make_display_communities(node_codes)

  est_fun <- build_clpn_estimator(
    node_codes=node_codes,
    covariates=covariates,
    weight_col=weight_col
  )

  net <- bootnet::estimateNetwork(
    data = as.matrix(frame),
    fun = est_fun,
    labels = node_codes,
    directed = TRUE
  )

  graph <- net$graph
  rownames(graph) <- node_codes
  colnames(graph) <- node_codes

  write.csv(graph, file.path(out_dir, paste0("edge_weights_", suffix, ".csv")))
  saveRDS(net, file.path(out_dir, paste0("network_object_", suffix, ".rds")))

  # Export node-level indices
  idx <- compute_directed_indices(graph, node_codes)
  write.csv(idx, file.path(out_dir, paste0("node_indices_", suffix, ".csv")), row.names=FALSE)

  # Export all edges and cross-community edges
  et <- graph_to_edge_table(graph, node_codes)
  write.csv(et, file.path(out_dir, paste0("all_edges_", suffix, ".csv")), row.names=FALSE)
  write.csv(
    et %>% filter(from_community != to_community),
    file.path(out_dir, paste0("cross_community_edges_", suffix, ".csv")),
    row.names=FALSE
  )
  write.csv(
    et %>% filter(
      (from_community == "Equanimity" & to_community != "Equanimity") |
      (from_community != "Equanimity" & to_community == "Equanimity")
    ),
    file.path(out_dir, paste0("equanimity_distress_edges_", suffix, ".csv")),
    row.names=FALSE
  )

  plot_clpn(graph, node_codes, out_dir, suffix=suffix, reduced=FALSE)
  plot_clpn(graph, node_codes, out_dir, suffix=suffix, reduced=TRUE)
  plot_indices(idx, out_dir, suffix=suffix)

  list(
    network=net, graph=graph, frame=frame, indices=idx,
    labels=labels, communities=communities, display_communities=display_communities, node_codes=node_codes
  )
}

graph_to_edge_table <- function(graph, node_codes) {
  disp <- setNames(make_display_communities(node_codes), node_codes)
  bridge <- setNames(make_bridge_communities(node_codes), node_codes)

  expand.grid(from=node_codes, to=node_codes, stringsAsFactors=FALSE) %>%
    mutate(
      weight = mapply(function(a,b) graph[a,b], from, to),
      abs_weight = abs(weight),
      from_community = disp[from],
      to_community = disp[to],
      from_bridge_community = bridge[from],
      to_bridge_community = bridge[to],
      autoregressive = from == to
    ) %>%
    arrange(desc(abs_weight))
}

compute_directed_indices <- function(graph, node_codes) {
  bridge_comm <- setNames(make_bridge_communities(node_codes), node_codes)
  display_comm <- setNames(make_display_communities(node_codes), node_codes)

  out_ei <- rowSums(graph)
  in_ei <- colSums(graph)
  out_strength <- rowSums(abs(graph))
  in_strength <- colSums(abs(graph))

  bridge_ei <- bridge_strength <- numeric(length(node_codes))
  names(bridge_ei) <- names(bridge_strength) <- node_codes

  for (i in seq_along(node_codes)) {
    node <- node_codes[i]
    other <- node_codes[bridge_comm[node_codes] != bridge_comm[node]]
    bridge_ei[i] <- sum(graph[node, other, drop=TRUE])
    bridge_strength[i] <- sum(abs(graph[node, other, drop=TRUE]))
  }

  data.frame(
    node=node_codes,
    community=unname(display_comm[node_codes]),
    bridge_community=unname(bridge_comm[node_codes]),
    out_expected_influence=unname(out_ei[node_codes]),
    in_expected_influence=unname(in_ei[node_codes]),
    out_strength=unname(out_strength[node_codes]),
    in_strength=unname(in_strength[node_codes]),
    bridge_expected_influence=unname(bridge_ei[node_codes]),
    bridge_strength=unname(bridge_strength[node_codes])
  )
}

plot_clpn <- function(graph, node_codes, out_dir, suffix="main", reduced=FALSE) {
  g <- graph
  comm <- make_display_communities(node_codes)
  groups <- make_group_index(node_codes)
  labels <- node_codes

  if (reduced) {
    # Hide within-equanimity, depression-depression, anxiety-anxiety, and
    # depression-anxiety edges. The full model is still estimated unchanged.
    for (i in seq_along(node_codes)) {
      for (j in seq_along(node_codes)) {
        if (!((comm[i] == "Equanimity" && comm[j] != "Equanimity") ||
              (comm[i] != "Equanimity" && comm[j] == "Equanimity"))) {
          g[i,j] <- 0
        }
      }
    }
  }

  pdf(
    file.path(out_dir, paste0(ifelse(reduced,"network_reduced_","network_full_"), suffix, ".pdf")),
    width=11, height=9
  )
  qgraph::qgraph(
    g,
    groups=groups,
    labels=labels,
    directed=TRUE,
    layout="spring",
    minimum=PLOT_MIN_EDGE,
    legend=TRUE,
    details=TRUE,
    posCol="#4A90E2",
    negCol="#D95F5F",
    color=c("#CFE8D6","#F6D1C1","#D7E7F5")
  )
  dev.off()
}

plot_indices <- function(idx, out_dir, suffix="main") {
  long <- idx %>%
    select(node, community,
           out_expected_influence, in_expected_influence,
           out_strength, in_strength,
           bridge_expected_influence, bridge_strength) %>%
    pivot_longer(-c(node,community), names_to="metric", values_to="value")

  p <- ggplot(long, aes(x=value, y=reorder(node, value))) +
    geom_col() +
    facet_wrap(~metric, scales="free_x", ncol=2) +
    theme_minimal(base_size=10) +
    labs(x=NULL, y=NULL)

  ggsave(file.path(out_dir, paste0("node_indices_", suffix, ".pdf")),
         p, width=11, height=12)
}

run_network_bootstrap <- function(result, out_dir, suffix="main",
                                  n_boot=N_BOOT, n_cores=N_CORES) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  stats <- c(
    "edge",
    "outExpectedInfluence","inExpectedInfluence",
    "outStrength","inStrength",
    "bridgeExpectedInfluence","bridgeStrength"
  )

  set.seed(SEED)
  b_nonparam <- bootnet::bootnet(
    result$network,
    nBoots=n_boot,
    nCores=n_cores,
    type="nonparametric",
    directed=TRUE,
    includeDiagonal=TRUE,
    communities=result$communities,
    statistics=stats,
    memorysaver=FALSE
  )

  set.seed(SEED + 1L)
  b_case <- bootnet::bootnet(
    result$network,
    nBoots=n_boot,
    nCores=n_cores,
    type="case",
    directed=TRUE,
    includeDiagonal=TRUE,
    communities=result$communities,
    statistics=c(
      "outExpectedInfluence","inExpectedInfluence",
      "outStrength","inStrength",
      "bridgeExpectedInfluence","bridgeStrength"
    ),
    memorysaver=FALSE
  )

  saveRDS(b_nonparam, file.path(out_dir, paste0("bootstrap_nonparametric_", suffix, ".rds")))
  saveRDS(b_case, file.path(out_dir, paste0("bootstrap_case_", suffix, ".rds")))

  cs <- bootnet::corStability(b_case)
  write.csv(
    data.frame(metric=names(cs), cs_coefficient=as.numeric(cs)),
    file.path(out_dir, paste0("cs_coefficients_", suffix, ".csv")),
    row.names=FALSE
  )

  pdf(file.path(out_dir, paste0("edge_accuracy_", suffix, ".pdf")), width=11, height=9)
  print(plot(b_nonparam, "edge", labels=FALSE, order="sample"))
  dev.off()

  pdf(file.path(out_dir, paste0("centrality_case_stability_", suffix, ".pdf")), width=12, height=9)
  print(plot(
    b_case,
    c("outExpectedInfluence","inExpectedInfluence","outStrength","inStrength",
      "bridgeExpectedInfluence","bridgeStrength"),
    facet=TRUE
  ))
  dev.off()

  # Full bootstrap tables
  write.csv(b_nonparam$sampleTable,
            file.path(out_dir, paste0("bootstrap_sample_table_", suffix, ".csv")), row.names=FALSE)
  write.csv(b_nonparam$bootTable,
            file.path(out_dir, paste0("bootstrap_full_table_", suffix, ".csv")), row.names=FALSE)

  export_bootstrap_ci(b_nonparam, out_dir, suffix)
  export_selection_frequency(b_nonparam, out_dir, suffix)
  export_equanimity_node_difference_tests(b_nonparam, out_dir, suffix)

  invisible(list(nonparametric=b_nonparam, case=b_case, cs=cs))
}

export_bootstrap_ci <- function(bootobj, out_dir, suffix) {
  bt <- bootobj$bootTable
  ci <- bt %>%
    group_by(type, id) %>%
    summarise(
      boot_mean=mean(value, na.rm=TRUE),
      boot_median=median(value, na.rm=TRUE),
      ci_low=quantile(value, .025, na.rm=TRUE, type=6),
      ci_high=quantile(value, .975, na.rm=TRUE, type=6),
      .groups="drop"
    ) %>%
    left_join(
      bootobj$sampleTable %>% select(type, id, sample_value=value),
      by=c("type","id")
    )
  write.csv(ci, file.path(out_dir, paste0("bootstrap_95CI_", suffix, ".csv")), row.names=FALSE)
}

export_selection_frequency <- function(bootobj, out_dir, suffix) {
  bt <- bootobj$bootTable %>% filter(type == "edge")
  sel <- bt %>%
    group_by(id) %>%
    summarise(
      proportion_nonzero=mean(abs(value) > 1e-12, na.rm=TRUE),
      proportion_positive=mean(value > 0, na.rm=TRUE),
      proportion_negative=mean(value < 0, na.rm=TRUE),
      .groups="drop"
    ) %>%
    left_join(
      bootobj$sampleTable %>% filter(type=="edge") %>% select(id, sample_weight=value),
      by="id"
    ) %>%
    arrange(desc(abs(sample_weight)))
  write.csv(sel, file.path(out_dir, paste0("edge_selection_frequency_", suffix, ".csv")), row.names=FALSE)
}

generic_boot_difference <- function(bootobj, id1, id2, metric, alpha=.05) {
  bt <- bootobj$bootTable %>% filter(type == metric, id %in% c(id1,id2))
  wide <- bt %>%
    select(name, id, value) %>%
    tidyr::pivot_wider(names_from=id, values_from=value)

  if (!all(c(id1,id2) %in% names(wide))) {
    return(data.frame(id1=id1,id2=id2,metric=metric,lower=NA,upper=NA,significant=NA))
  }
  d <- wide[[id2]] - wide[[id1]]
  lo <- quantile(d, alpha/2, na.rm=TRUE, type=6)
  hi <- quantile(d, 1-alpha/2, na.rm=TRUE, type=6)
  data.frame(
    id1=id1, id2=id2, metric=metric,
    lower=unname(lo), upper=unname(hi),
    significant=!(lo <= 0 && hi >= 0)
  )
}

export_equanimity_node_difference_tests <- function(bootobj, out_dir, suffix) {
  metrics <- c(
    "outExpectedInfluence","inExpectedInfluence",
    "outStrength","inStrength",
    "bridgeExpectedInfluence","bridgeStrength"
  )
  eq_ids <- EQ_CODES

  out <- list()
  ctr <- 1L
  for (m in metrics) {
    valid <- intersect(eq_ids, unique(bootobj$bootTable$id[bootobj$bootTable$type == m]))
    if (length(valid) >= 2) {
      cmb <- combn(valid, 2, simplify=FALSE)
      for (z in cmb) {
        out[[ctr]] <- generic_boot_difference(bootobj, z[1], z[2], m)
        ctr <- ctr + 1L
      }
    }
  }
  if (length(out) > 0) {
    tab <- bind_rows(out)
    write.csv(tab, file.path(out_dir, paste0("equanimity_node_difference_tests_", suffix, ".csv")), row.names=FALSE)
  }
}

compare_networks <- function(graph_a, graph_b, label_a="A", label_b="B") {
  stopifnot(all(dim(graph_a) == dim(graph_b)))
  va <- as.vector(graph_a)
  vb <- as.vector(graph_b)
  nz_a <- abs(va) > 1e-12
  nz_b <- abs(vb) > 1e-12
  both <- nz_a | nz_b

  data.frame(
    model_a=label_a,
    model_b=label_b,
    edge_weight_correlation=cor(va, vb, use="pairwise.complete.obs"),
    edge_weight_correlation_nonzero_union=if (sum(both)>2) cor(va[both], vb[both]) else NA,
    selection_jaccard=sum(nz_a & nz_b) / max(1, sum(nz_a | nz_b)),
    sign_agreement_selected_union=mean(sign(va[both]) == sign(vb[both])),
    n_nonzero_a=sum(nz_a),
    n_nonzero_b=sum(nz_b)
  )
}
