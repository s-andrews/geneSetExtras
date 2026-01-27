#' Plot a volcano plot of categorical enrichment results
#'
#' @param enrich_result A clusterProfiler categorical enrichment results from enrichGO or similar
#' @param interactive Whether to return a traditional ggplot (FALSE) or interactive plotly (TRUE) plot
#'
#' @return A ggplot or plotly plot
#' @export
#'
#' @examples
volcano_plot <- function(enrich_result, interactive=FALSE, title=NA) {
  
  enrich_result |>
    dplyr::as_tibble() |>
    dplyr::mutate(setSize=as.integer(stringr::str_replace(BgRatio,"/.*$",""))) |>
    dplyr::mutate(significant = dplyr::if_else(p.adjust<0.05,"Significant","Non-Significant")) |>
    dplyr::mutate(FDR_Phred = -10*log10(p.adjust)) |>
    arrange(desc(p.adjust)) -> plot_data
  
  
  plot_data |>
    ggplot2::ggplot(aes(x=FoldEnrichment, y=FDR_Phred, size=setSize, ID=ID, Description=Description, colour=significant)) +
    ggplot2::geom_point(pch=21) +
    scale_colour_manual(values=c("grey","red2")) -> volcano_plot
  
  if (!is.na(title)) {
    volcano_plot + ggtitle(title) -> volcano_plot
  }
  
  
  if (interactive) {
      return(plotly::ggplotly(volcano_plot, width=800, height=500))
  }
  else {
    return(volcano_plot)
  }
  
}