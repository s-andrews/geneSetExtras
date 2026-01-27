#' Create a plot of enrichment vs significance for clustered results
#'
#' @param cluster_result An output from cluster_enrich_result
#' @param interactive Whether to return a traditional ggplot (FALSE) or interactive plotly (TRUE) plot
#'
#' @return A ggplot or plotly plot
#' @export
#'
#' @examples
cluster_volcano_plot <- function(cluster_result, interactive=FALSE) {
  cluster_result |>
    mutate(FDR_Phred = -10*log10(p.adjust)) -> plot_data
  
  plot_data |>
    ggplot2::ggplot(aes(x=FoldEnrichment, y=FDR_Phred, size=term_cluster_size, ID=ID, Description=Description)) +
    ggplot2::geom_point(pch=21) +
    ggplot2::coord_cartesian(ylim=c(0,max(plot_data$FDR_Phred)))-> volcano_plot
  
  
  if(interactive) {
    return(plotly::ggplotly(volcano_plot, width=800, height=800))
  }
  else{
    return(volcano_plot)
  }
}