plot_differential_enrichment_heatmap <- function(differential_enrichment_hits, enrichment_results, group_names) {
  
  # We just want the ids from the differential enrichment
  differential_enrichment_hits$ID -> ids_to_plot

  enrichment_plot_data <- list()
  
  for (i in 1:length(group_names)) {
    enrichment_results[[i]] |>
      as_tibble() |>
      filter(ID %in% ids_to_plot) |>
      add_column(gene_set=group_names[i]) -> enrichment_plot_data[[i]]
  }
  
  do.call(bind_rows,enrichment_plot_data) -> enrichment_plot_data
  
  tidyheatmaps::tidy_heatmap(
    enrichment_plot_data,
    rows=Description,
    columns = gene_set,
    values=FoldEnrichment,
    cluster_rows=TRUE,
    cluster_cols = TRUE,
    scale="none",
    colors = c("white","red3")
  )
  
}