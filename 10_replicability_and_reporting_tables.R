source("00_config.R")

out <- file.path(OUTPUT_DIR, "10_replicability_and_reporting_tables")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

primary_dir <- file.path(OUTPUT_DIR, "05_combined_clpn_primary")
boot_file <- file.path(primary_dir, "bootstrap_nonparametric_combined_primary_unweighted.rds")
result_file <- file.path(primary_dir, "combined_primary_result.rds")

if (!file.exists(boot_file) || !file.exists(result_file)) {
  stop("Run 05_combined_clpn_primary.R first.")
}

bootobj <- readRDS(boot_file)
res <- readRDS(result_file)

sel_file <- file.path(primary_dir, "edge_selection_frequency_combined_primary_unweighted.csv")
sel <- read.csv(sel_file, check.names=FALSE)

rep_summary <- data.frame(
  threshold=c(.50,.70,.80,.90),
  n_edges_with_selection_frequency_at_least=
    sapply(c(.50,.70,.80,.90), function(t) sum(sel$proportion_nonzero >= t))
)
write.csv(rep_summary, file.path(out, "bootstrap_edge_replicability_summary.csv"), row.names=FALSE)

matched <- read.csv(
  file.path(OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv"),
  check.names=FALSE
)

desc_vars <- c(
  "T1_EAC","T1_NRT","T1_EMS","T1_HIN","T1_PHQ_TOTAL","T1_GAD_TOTAL",
  "T2_EAC","T2_NRT","T2_EMS","T2_HIN","T2_PHQ_TOTAL","T2_GAD_TOTAL"
)
write.csv(
  make_descriptives(matched, desc_vars),
  file.path(out, "descriptive_statistics_core_variables.csv"),
  row.names=FALSE
)

cors <- cor(matched[,desc_vars], use="pairwise.complete.obs", method="pearson")
write.csv(cors, file.path(out, "zero_order_correlations_core_variables.csv"))

edges <- graph_to_edge_table(res$graph, c(EQ_CODES,SYMPTOM_CODES)) %>%
  filter(abs(weight) > 1e-12)

write.csv(edges, file.path(out, "all_nonzero_cross_lagged_edges.csv"), row.names=FALSE)
write.csv(edges %>% filter(autoregressive),
          file.path(out, "autoregressive_edges.csv"), row.names=FALSE)
write.csv(
  edges %>% filter(
    (from_community=="Equanimity" & to_community!="Equanimity") |
    (from_community!="Equanimity" & to_community=="Equanimity")
  ),
  file.path(out, "nonzero_equanimity_distress_edges.csv"),
  row.names=FALSE
)

sink(file.path(out, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("10_replicability_and_reporting_tables v2 complete: ", normalizePath(out))
