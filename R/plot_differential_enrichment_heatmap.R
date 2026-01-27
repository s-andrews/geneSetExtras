#' Plot the results of differential enrichment as a heatmap
#'
#' @param differential_enrichment_hits The output of compare_enrichment_results
#' @param scale A TRUE/FALSE value to say whether to scale the heatmap per row
#' @param log_scale Whether to log tranform the enrichment values
#'
#' @return
#' @export
#'
#' @examples
plot_differential_enrichment_heatmap <- function(differential_enrichment_hits,scale=TRUE, log_scale=TRUE) {
  
  differential_enrichment_hits |>
    pivot_longer(
      cols=starts_with("FoldEnrichment_"),
      names_to="group",
      values_to="Enrichment"
    ) |>
    mutate(group=str_replace(group,"^FoldEnrichment_","")) -> plot_data

  if (log_scale) {
    plot_data |>
      mutate(Enrichment = replace(Enrichment, Enrichment==0, 0.1)) |>
      mutate(Enrichment=log2(Enrichment)) -> plot_data
  }  

  scale_value <- "none"

  if(scale) {
    scale_value="row"
  }
  

  tidyheatmaps::tidy_heatmap(
    plot_data,
    rows=Description,
    columns = group,
    values=Enrichment,
    cluster_rows=TRUE,
    cluster_cols = TRUE,
    scale=scale_value
  )
  
}