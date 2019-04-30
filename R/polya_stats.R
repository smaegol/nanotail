#' Calculate basic statistics
#'
#' @param polya_data
#' @param min_reads
#'
#' @return tibble with statistics
#' @export
#'
#' @examples
calculate_polya_stats <- function(polya_data, min_reads = 0, grouping_factor = "group") {



  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assert_that(is.number(min_reads),msg = "Non-numeric parameter provided (min_reads)")

  ### TBD
  # batch - make conditional (testing as formula??)
  # remove batch for now as it is raising some Evaluation error: contrasts can be applied only to factors with 2 or more levels. errors

  ###

  #test_formula = reformulate(grouping_factor,polya_length)

  # leave only those tanscripts which were identified in all conditions
  polya_data_complete_cases <-
    polya_data %>% dplyr::group_by(.dots = c("transcript", grouping_factor)) %>% dplyr::add_count() %>% dplyr::filter(n > min_reads) %>% dplyr::slice(1) %>%
    dplyr::ungroup() %>% dplyr::group_by(transcript) %>% dplyr::add_count() %>% dplyr::filter(nn > 1)
  # All 0 values in polya_length transformed to 1 (log2(1)=0) to avoid errors in glm call
  polya_data_complete <-
    polya_data %>% dplyr::semi_join(polya_data_complete_cases %>% dplyr::select(transcript)) %>% dplyr::mutate(polya_length=ifelse(polya_length==0,1,polya_length))
  # calculate statistics using Wilcoxon test (non-parametric)
  polya_data_stat_wilcoxon <-
    polya_data_complete %>% dplyr::ungroup() %>% dplyr::group_by(transcript) %>% dplyr::mutate(
      stats = suppressWarnings(wilcox.test(as.formula(paste0("polya_length ~",grouping_factor))))$p.value,
      ks_stats = suppressWarnings(FSA::ksTest(as.formula(paste0("polya_length ~",grouping_factor))))$p.value)


  # summarise statistics
  polyA_data_stat_wilcoxon_summary <-
    polya_data_stat_wilcoxon %>% dplyr::group_by(.dots = c("transcript", grouping_factor)) %>% dplyr::summarise(
      p.value = max(stats),
      p.value_ks = max(ks_stats),
      counts = dplyr::n(),
      polya_mean = mean(polya_length),
      polya_sd = sd(polya_length),
      polya_median = median(polya_length),
      polya_gm_mean = gm_mean(polya_length)
    ) %>% tidyr::gather(
      key = variable,
      value = wart,
      counts,
      polya_mean,
      polya_gm_mean,
      polya_median,
      polya_sd,
      p.value,
      p.value_ks
    ) %>% dplyr::mutate(group_var = paste(group, variable, sep = "_")) %>% dplyr::select(-c(group,
                                                                                     variable)) %>% tidyr::spread(group_var, wart) %>% dplyr::rename(p.value = wt_p.value,
                                                                                                                                       p.value_ks = wt_p.value_ks) %>% dplyr::select(-mut_p.value) %>%
    dplyr::mutate(fold_change = wt_polya_mean / mut_polya_mean) %>% dplyr::select(
      transcript,
      wt_counts,
      mut_counts,
      wt_polya_mean,
      mut_polya_mean,
      wt_polya_gm_mean,
      mut_polya_gm_mean,
      wt_polya_median,
      mut_polya_median,
      wt_polya_sd,
      mut_polya_sd,
      fold_change_means = fold_change,
      p.value,
      p.value_ks
    )
  pvals_BH_summary <-
    p.adjust(polyA_data_stat_wilcoxon_summary$p.value, method = "BH")
  polyA_data_stat_wilcoxon_summary$p.corr <- pvals_BH_summary
  pvals_BH_summary_KS <-
    p.adjust(polyA_data_stat_wilcoxon_summary$p.value_ks, method = "BH")
  polyA_data_stat_wilcoxon_summary$p.corr_ks <-
    pvals_BH_summary_KS
 # pvals_BH_summary_glm <-
   # p.adjust(polyA_data_stat_wilcoxon_summary$p.value_glm, method = "BH")
 # polyA_data_stat_wilcoxon_summary$p.corr_glm <-
    #pvals_BH_summary_glm


  return(
    list(stats = polya_data_stat_wilcoxon, summary = polyA_data_stat_wilcoxon_summary)
  )
}

#' Calculate basic statistics
#'
#' @param polya_data
#' @param min_reads
#'
#' @return tibble with statistics
#' @export
#'
#' @examples
calculate_polya_stats_glm <- function(polya_data, min_reads = 0, grouping_factors = c("group")) {



  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assert_that(is.number(min_reads),msg = "Non-numeric parameter provided (min_reads)")

  ### TBD
  # batch - make conditional (testing as formula??)
  # remove batch for now as it is raising some Evaluation error: contrasts can be applied only to factors with 2 or more levels. errors

  ###

  # leave only those tanscripts which were identified in all conditions
  polya_data_complete_cases <-
    polya_data %>% dplyr::group_by(transcript, group) %>% dplyr::add_count() %>% dplyr::filter(n > min_reads) %>% dplyr::slice(1) %>%
    dplyr::ungroup() %>% dplyr::group_by(transcript) %>% dplyr::add_count() %>% dplyr::filter(nn > 1)
  # All 0 values in polya_length transformed to 1 (log2(1)=0) to avoid errors in glm call
  polya_data_complete <-
    polya_data %>% dplyr::semi_join(polya_data_complete_cases %>% dplyr::select(transcript)) %>% dplyr::mutate(polya_length=ifelse(polya_length==0,1,polya_length))
  # calculate statistics using Wilcoxon test (non-parametric)
  polya_data_stat_wilcoxon <-
    polya_data_complete %>% dplyr::ungroup() %>% dplyr::group_by(transcript) %>% dplyr::mutate(
      stats = suppressWarnings(wilcox.test(polya_length ~
                                             as.factor(group)))$p.value,
      ks_stats = suppressWarnings(FSA::ksTest(polya_length ~ group))$p.value,
      glm_stats = coef(summary(glm(
        formula = log2(polya_length) ~ group)))[2, 4])


  # summarise statistics
  polyA_data_stat_wilcoxon_summary <-
    polya_data_stat_wilcoxon %>% dplyr::group_by(transcript, group) %>% dplyr::summarise(
      p.value = max(stats),
      p.value_ks = max(ks_stats),
      p.value_glm = max(glm_stats),
      counts = dplyr::n(),
      polya_mean = mean(polya_length),
      polya_sd = sd(polya_length),
      polya_median = median(polya_length),
      polya_gm_mean = gm_mean(polya_length)
    ) %>% tidyr::gather(
      key = variable,
      value = wart,
      counts,
      polya_mean,
      polya_gm_mean,
      polya_median,
      polya_sd,
      p.value,
      p.value_ks,
      p.value_glm
    ) %>% dplyr::mutate(group_var = paste(group, variable, sep = "_")) %>% dplyr::select(-c(group,
                                                                                            variable)) %>% tidyr::spread(group_var, wart) %>% dplyr::rename(p.value = wt_p.value,
                                                                                                                                                            p.value_ks = wt_p.value_ks,
                                                                                                                                                            p.value_glm = wt_p.value_glm) %>% dplyr::select(-mut_p.value) %>%
    dplyr::mutate(fold_change = wt_polya_mean / mut_polya_mean) %>% dplyr::select(
      transcript,
      wt_counts,
      mut_counts,
      wt_polya_mean,
      mut_polya_mean,
      wt_polya_gm_mean,
      mut_polya_gm_mean,
      wt_polya_median,
      mut_polya_median,
      wt_polya_sd,
      mut_polya_sd,
      fold_change_means = fold_change,
      p.value,
      p.value_ks,
      p.value_glm
    )
  pvals_BH_summary <-
    p.adjust(polyA_data_stat_wilcoxon_summary$p.value, method = "BH")
  polyA_data_stat_wilcoxon_summary$p.corr <- pvals_BH_summary
  pvals_BH_summary_KS <-
    p.adjust(polyA_data_stat_wilcoxon_summary$p.value_ks, method = "BH")
  polyA_data_stat_wilcoxon_summary$p.corr_ks <-
    pvals_BH_summary_KS
  pvals_BH_summary_glm <-
    p.adjust(polyA_data_stat_wilcoxon_summary$p.value_glm, method = "BH")
  polyA_data_stat_wilcoxon_summary$p.corr_glm <-
    pvals_BH_summary_glm


  return(
    list(stats = polya_data_stat_wilcoxon, summary = polyA_data_stat_wilcoxon_summary)
  )
}



#' Title
#'
#' @param polya_data
#'
#' @return long-format tibble with per-transcript statistics for each sample
#' @export
#'
#' @examples
summarize_polya <- function(polya_data,summary_factors = c("group")) {

  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assert_that(has_rows(polya_data),msg = "Empty data.frame provided as an input")
  assert_that(is_character(summary_factors),msg = "Non-character argument is not alowed for `summary factors`. Please provide either string or vector of strings")
  assert_that(all(summary_factors %in% colnames(polya_data)),msg="Non-existent column name provided as the argument (summary_factors)")

  polya_data_summarized <-
    polya_data %>% dplyr::ungroup() %>% dplyr::group_by(.dots = c("transcript",summary_factors)) %>% dplyr::summarise(
      counts = dplyr::n(),
      polya_mean = mean(polya_length),
      polya_sd = sd(polya_length),
      polya_median = median(polya_length),
      polya_gm_mean = gm_mean(polya_length)
    )
  return(polya_data_summarized)
}
#' title
#'
#' @param polya_data
#'
#' @return
#' @export
#'
#' @examples
polya_summary_pivot_to_wide <- function(polya_data) {
  polya_data_summarized <-
    polya_data %>% dplyr::ungroup() %>% dplyr::group_by(sample_name, transcript, ensembl_transcript_id_short) %>% dplyr::summarise(
      counts = dplyr::n(),
      polya_mean = mean(polya_length),
      polya_sd = sd(polya_length),
      polya_median = median(polya_length),
      polya_gm_mean = gm_mean(polya_length)
    ) %>% tidyr::gather(
      key = variable,
      value = wart,
      counts,
      polya_mean,
      polya_gm_mean,
      polya_median,
      polya_sd) %>% dplyr::mutate(group_var = paste(sample_name, variable, sep = "_")) %>% dplyr::ungroup() %>% dplyr::select(-c(sample_name, variable)) %>% tidyr::spread(group_var, wart)
      return(polya_data_summarized)
}


#' Title
#'
#' @param polya_data_summarized
#'
#' @return
#' @export
#'
#' @examples
calculate_pca <- function(polya_data_summarized) {


}

#' Title
#'
#' @param polya_data
#'
#' @return
#' @export
#'
#' @examples
calculate_processing_statistics <- function(polya_data) {

  return_list <- list()
  return_list$number_all_reads <- nrow(polya_data)
  return_list$number_adapter_reads <- polya_data %>% dplyr::filter(qc_tag=='ADAPTER') %>% nrow()
  return_list$number_noregion_reads <- polya_data %>% dplyr::filter(qc_tag=='NOREGION') %>% nrow()
  return_list$number_read_load_failed_reads <- polya_data %>% dplyr::filter(qc_tag=='READ_FAILED_LOAD') %>% nrow()
  return_list$number_suffclip_reads <- polya_data %>% dplyr::filter(qc_tag=='SUFFCLIP') %>% nrow()
  return_list$number_pass_reads <- polya_data %>% dplyr::filter(qc_tag=='PASS') %>% nrow()
  return (return_list)
}
