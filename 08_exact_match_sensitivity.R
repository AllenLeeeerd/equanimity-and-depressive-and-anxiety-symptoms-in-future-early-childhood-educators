source("00_config.R")

main_file <- file.path(OUTPUT_DIR, "01_data_preparation", "matched_analytic_clean.csv")
exact_file <- file.path(OUTPUT_DIR, "01_data_preparation", "exact_match_analytic_clean.csv")
if (!file.exists(main_file) || !file.exists(exact_file)) stop("Run 01_data_preparation.R first.")

main_dat <- read.csv(main_file, check.names=FALSE)
exact_dat <- read.csv(exact_file, check.names=FALSE)

out <- file.path(OUTPUT_DIR, "08_exact_match_sensitivity")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

nodes <- c(EQ_CODES, SYMPTOM_CODES)

main_res <- estimate_clpn(
  main_dat, nodes, out, PRIMARY_COVARIATES, NULL,
  suffix=paste0("analytic_reference_N",nrow(main_dat))
)
exact_res <- estimate_clpn(
  exact_dat, nodes, out, PRIMARY_COVARIATES, NULL,
  suffix=paste0("exact_match_sensitivity_N",nrow(exact_dat))
)

comp <- compare_networks(
  main_res$graph, exact_res$graph,
  paste0("Analytic N=", nrow(main_dat)),
  paste0("Exact match N=", nrow(exact_dat))
)
write.csv(comp, file.path(out, "analytic_vs_exact_match_comparison.csv"), row.names=FALSE)

edge_main <- graph_to_edge_table(main_res$graph, nodes) %>%
  filter(
    (from_community=="Equanimity" & to_community!="Equanimity") |
    (from_community!="Equanimity" & to_community=="Equanimity")
  ) %>%
  select(from,to,weight_main=weight)

edge_exact <- graph_to_edge_table(exact_res$graph, nodes) %>%
  filter(
    (from_community=="Equanimity" & to_community!="Equanimity") |
    (from_community!="Equanimity" & to_community=="Equanimity")
  ) %>%
  select(from,to,weight_exact=weight)

edge_compare <- full_join(edge_main, edge_exact, by=c("from","to")) %>%
  mutate(
    retained_main=abs(weight_main)>1e-12,
    retained_exact=abs(weight_exact)>1e-12,
    same_sign=sign(weight_main)==sign(weight_exact)
  )
write.csv(edge_compare, file.path(out, "equanimity_distress_edge_comparison.csv"), row.names=FALSE)

message("08_exact_match_sensitivity v2 complete: ", normalizePath(out))
