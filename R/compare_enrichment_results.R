#' Test for differentially enriched categories between gene sets
#'
#' @param enrichment_results A list of results from enrichGO or similar
#'
#' @return A tibble of statistical results for differentially enriched categories
#' @export
#'
#' @examples
compare_enrichment_results <- function(enrichment_results) {
  
  # We'll get variable numbers of enrichment results from enrichGO or similar
  # For each one we need to extract out the IDs, Description and Gene and Background
  # Ratios and put them into a common tibble

  list() -> data_for_testing
  
  for (i in 1:length(enrichment_results)) {
    enrichment_results[[i]] |>
      dplyr::as_tibble() |>
      dplyr::select(ID,Description, GeneRatio, BgRatio) |>
      tidyr::separate(GeneRatio, into=c("GeneSuccess","GeneTotal"), sep="/", convert=TRUE) |>
      tidyr::separate(BgRatio, into=c("BgSuccess","BgTotal"), sep="/", convert=TRUE) |>
      tibble::add_column(result_id=i) -> data_for_testing[[i]]
  }
  
  do.call(dplyr::bind_rows, data_for_testing) -> data_for_testing
  
  # We have a problem here which is that if enrichGO doesn't have any hits to a category
  # in its hit list (not the background) then it doesn't report the result, not matter
  # what the filters are set to.  We can therefore only test categories which appear
  # at least once in all of the results unless we kludge the results.
  # 
  # This *may* be a fixed problem - the devel version of enrichGO uses enrichit as the
  # back end and this looks like it might report 0 overlap hits.  Need to wait for that
  # to hit production and check again.
  
  data_for_testing |>
    dplyr::group_by(ID) |>
    dplyr::filter(
      n() == length(enrichment_results)
    ) |>
    dplyr::ungroup() -> data_for_testing
  
  data_for_testing |>
    distinct(ID,Description) -> annotation_data
  
  # A bit more stuff to do
  
  data_for_testing |>
    dplyr::arrange(ID) |>
    dplyr::mutate(Gene_Fail = GeneTotal-GeneSuccess, Bg_Fail=BgTotal-BgSuccess) |>
    dplyr::select(ID,result_id, GeneSuccess,Gene_Fail, BgSuccess,Bg_Fail)  |>
    dplyr::rename(Gene_Success=GeneSuccess, Bg_Success=BgSuccess) |>
    tidyr::pivot_longer(
      cols=c(Gene_Success,Gene_Fail,Bg_Success,Bg_Fail),
      names_to="count_type",
      values_to="count"
    ) |>
    tidyr::separate(count_type, into=c("group","type"),sep="_") |>
    tidyr::pivot_wider(
      names_from=type,
      values_from=count
    ) |>
    dplyr::mutate(result_id=factor(result_id)) -> data_for_testing
  
  data_for_testing |>
    dplyr::group_by(ID) |>
    dplyr::group_modify(perform_test) -> results_data

  results_data |>
    dplyr::ungroup() |>
    dplyr::left_join(annotation_data) |>
    dplyr::arrange(pvalue) |>
    dplyr::mutate(p.adjust=p.adjust(pvalue,method="fdr")) |>
    dplyr::select(ID,Description,p.adjust,everything()) -> results_data
  
  return(results_data)
  
}


perform_test <- function(data,id) {

  stats::glm (
    cbind(Success,Fail) ~ result_id * group,
    data = data,
    family=binomial()
  ) -> fit
  
  stats::glm (
    cbind(Success,Fail) ~ result_id + group,
    data = data,
    family=binomial()
  ) -> fit_no_interaction
  
  stats::anova(fit_no_interaction,fit, test="Chisq") -> test_result
  
  test_result$`Pr(>Chi)`[-1] -> p_value
  
  return(tibble::tibble(pvalue=p_value))
  
}


