cluster_volcano_plot <- function(cluster_result, interactive=FALSE) {
  cluster_result |>
    mutate(FDR_Phred = -10*log10(p.adjust)) |>
    ggplot2::ggplot(aes(x=FoldEnrichment, y=FDR_Phred, size=term_cluster_size, ID=ID, Description=Description)) +
    ggplot2::geom_point(pch=21) -> volcano_plot
  
  
  if(interactive) {
    return(plotly::ggplotly(volcano_plot, width=800, height=800))
  }
  else{
    return(volcano_plot)
  }
}