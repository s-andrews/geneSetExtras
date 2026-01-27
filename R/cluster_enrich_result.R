#' Heirarchical clustering of enrichment results based on the similarity of terms within them
#'
#' @param enrich_result A clusterProfiler categorical enrichment results from enrichGO or similar
#' @param cluster_threshold The tree height at which to create clusters. Higher values give more clusters
#' @param min_fold_enrichment Fold Enrichment filter for which terms to cluster
#' @param max_p_adjust Adjusted P value filter for which terms to cluster
#' @param clust_method How cluster agglomeration method to use - see hclust for details
#' @param min_cluster_size The smallest cluster term size to report
#'
#' @return
#' @export
#'
#' @examples
cluster_enrich_result <- function(enrich_result, cluster_threshold=0.9, min_fold_enrichment=1, max_p_adjust=0.05, clust_method="average", min_cluster_size=2) {
  enrich_result |>
    dplyr::filter(FoldEnrichment>=min_fold_enrichment) |>
    dplyr::filter(p.adjust<=max_p_adjust) -> data_for_similarity
  
  if (nrow(data_for_similarity) < 2) {
    return(NULL)
  }
  data_for_similarity |>
    enrichplot::pairwise_termsim() -> termsim_results
  
  
  stats::as.dist(1-termsim_results@termsim) |>
    stats::hclust(method = clust_method) |>
    stats::cutree(h=cluster_threshold) -> clusters
  
  tibble::tibble(
    Description=names(clusters),
    term_cluster=clusters
  ) |>
    dplyr::group_by(term_cluster) |>
    dplyr::mutate(term_cluster_size=n()) |>
    dplyr::filter(term_cluster_size >= min_cluster_size) -> clusters
  
  
  enrich_result |>
    tibble::as_tibble() |>
    dplyr::inner_join(clusters) -> clusters
  
  clusters |>
    dplyr::arrange(term_cluster_size,pvalue) |>
    dplyr::group_by(term_cluster) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(ID,Description,term_cluster,term_cluster_size,everything()) -> final_data
  
  return(final_data)
    
}