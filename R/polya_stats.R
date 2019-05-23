#' Calculates basic statistics for polya lengths
#'
#' Takes polyA predictions table as input and checks if there is significant difference in polyA lengths between chosen conditions for each transcript.
#' By default, Wilcoxon Rank Sum (\link{wilcox.test}) test is used.
#'
#' @param polya_data input table with polyA predictions
#' @param min_reads minimum number of reads to include transcript in the analysis
#' @param grouping_factor which column defines groups (default: sample_name)
#' @param condition1 if `grouping_factor` has more than 2 levels, which level use for comparison
#' @param condition2 if `grouping_factor` has more than 2 levels, which level use for comparison
#' @param stat_test what statistical test to use for testing, currently supports "Wilcoxon" (for \link{wilcox.test}), "KS" (for \link[FSA]{ksTest} from FSA package) or "glm" (for \link{glm})
#' @param use_dwell_time if TRUE, will use dwell time instead of estimated polya length for statistics
#' @param custom_glm_formula provides custom formula to be used with glm statistical test
#'
#' @return A list with elements:\itemize{
#' \item summary - summary table with pvalues and median/mean values associated to each transcript
#' \item stats - input table with p.value column added}
#'
#' @export
#'
calculate_polya_stats <- function(polya_data, min_reads = 0, grouping_factor = "sample_name",condition1=NA,condition2=NA,stat_test="Wilcoxon",use_dwell_time=FALSE,custom_glm_formula,alpha=0.05) {



  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }


  available_statistical_tests = c("Wilcoxon","KS","glm")

  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data.frame provided as an input")

  assertthat::assert_that(stat_test %in% available_statistical_tests,msg = "Please provide one of available statistical tests (Wilcoxon, KS or glm)")
  assertthat::assert_that(assertthat::is.number(min_reads),msg = "Non-numeric parameter provided (min_reads)")
  assertthat::assert_that(assertive::is_a_bool(use_dwell_time),msg="Non-boolean value provided for option use_dwell_time.")
  assertthat::assert_that(grouping_factor %in% colnames(polya_data),msg=paste0(grouping_factor," is not a column of input dataset"))

  # if grouping factor has more than two levels
  if (length(levels(polya_data[[grouping_factor]]))>2) {
    if(is.na(condition1) && is.na(condition2)) {
      #throw error when no conditions for comparison are specified
      stop(paste0("grouping_factor ",grouping_factor," has more than 2 levels. Please specify condtion1 and condition2 to select comparison pairs"))
    }
    else {
      # filter input data leaving only specified conditions, dropping other factor levels
      assertthat::assert_that(condition1 %in% levels(polya_data[[grouping_factor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[grouping_factor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (grouping_factor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(grouping_factor) %in% c(condition1,condition2)) %>% dplyr::mutate() %>% droplevels()
    }
  }
  else if (length(levels(polya_data[[grouping_factor]]))==1) {
    stop("Only 1 level present for grouping factor. Choose another groping factor for comparison")
  }
  else {
    condition1 = levels(polya_data[[grouping_factor]])[1]
    condition2 = levels(polya_data[[grouping_factor]])[2]
  }

  ### TBD
  # batch - make conditional (testing as formula??)
  # remove batch for now as it is raising some Evaluation error: contrasts can be applied only to factors with 2 or more levels. errors

  ###

  #test_formula = reformulate(grouping_factor,polya_length)

  if (!missing(custom_glm_formula)) {
    glm_groups<-all.vars(as.formula(custom_glm_formula))[-1]
    polya_data_complete_cases <-
      polya_data %>% dplyr::group_by(.dots = c("transcript", glm_groups)) %>% dplyr::add_count() %>% dplyr::filter(n > min_reads) %>% dplyr::slice(1) %>%
      dplyr::ungroup() %>% dplyr::group_by(transcript) %>% dplyr::summarize(nn = n()) %>% dplyr::filter(nn > 2^(length(glm_groups)-1))
  }


  else {
  # leave only those tanscripts which were identified in all conditions
  polya_data_complete_cases <-
    polya_data %>% dplyr::group_by(.dots = c("transcript", grouping_factor)) %>% dplyr::add_count() %>% dplyr::filter(n > min_reads) %>% dplyr::slice(1) %>%
    dplyr::ungroup() %>% dplyr::group_by(transcript) %>% dplyr::summarize(nn = n()) %>% dplyr::filter(nn > 1)
  }
  # All 0 values in polya_length transformed to 1 (log2(1)=0) to avoid errors in glm call
  polya_data_complete <-
    polya_data %>% dplyr::semi_join(polya_data_complete_cases %>% dplyr::select(transcript)) %>% dplyr::mutate(polya_length=ifelse(polya_length==0,1,polya_length))


  polya_data_stat <-
    polya_data_complete %>% dplyr::ungroup() %>% dplyr::group_by(transcript)


  if (use_dwell_time) {
    statistics_formula <- paste0("dwell_time ~",grouping_factor)
  }
  else {
    statistics_formula <- paste0("polya_length ~",grouping_factor)
  }

  if (stat_test=="Wilcoxon") {
    polya_data_stat <- polya_data_stat  %>% dplyr::mutate(stats = suppressWarnings(wilcox.test(as.formula(statistics_formula))$p.value))
  }
  else if (stat_test=="KS") {
    polya_data_stat <- polya_data_stat  %>% dplyr::mutate(stats = suppressWarnings(FSA::ksTest(as.formula(statistics_formula))$p.value))
  }
  else if (stat_test=="glm") {
    if(!missing(custom_glm_formula)) {
      statistics_formula = custom_glm_formula
    }
    mcp_call <- paste0("multcomp::mcp(",grouping_factor,' = "Tukey")')
    polya_data_stat <- polya_data_stat  %>% dplyr::mutate(stats = suppressWarnings(summary(multcomp::glht(glm(formula = as.formula(statistics_formula)),eval(parse(text = mcp_call))))$test$pvalues[1]))
    #polya_data_stat <- polya_data_stat  %>% dplyr::mutate(stats = suppressWarnings(coef(summary(glm(formula = as.formula(statistics_formula))))[2,4]))
  }
  else {
    stop("wrong stat_test parameter provided")
  }

  # summarise statistics
  polyA_data_stat_summary <-
    polya_data_stat %>% dplyr::group_by(.dots = c("transcript", grouping_factor)) %>% dplyr::summarise(
      p.value = max(stats),
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
      p.value
    ) %>% dplyr::mutate(group_var = paste(!!rlang::sym(grouping_factor), variable, sep = "_")) %>% dplyr::select(-c(!!rlang::sym(grouping_factor),
                                                                                                                    variable)) %>% tidyr::spread(group_var, wart) %>% dplyr::rename(p.value = !!rlang::sym(paste0(condition1,"_p.value"))) %>% dplyr::select(- !!rlang::sym(paste0(condition2,"_p.value"))) %>%
    dplyr::mutate(fold_change = !!rlang::sym(paste0(condition1,"_polya_median")) / !!rlang::sym(paste0(condition2,"_polya_median"))) %>% dplyr::select(
      transcript,
      !!rlang::sym(paste0(condition1,"_counts")),
      !!rlang::sym(paste0(condition2,"_counts")),
      !!rlang::sym(paste0(condition1,"_polya_mean")),
      !!rlang::sym(paste0(condition2,"_polya_mean")),
      !!rlang::sym(paste0(condition1,"_polya_gm_mean")),
      !!rlang::sym(paste0(condition2,"_polya_gm_mean")),
      !!rlang::sym(paste0(condition1,"_polya_median")),
      !!rlang::sym(paste0(condition2,"_polya_median")),
      !!rlang::sym(paste0(condition1,"_polya_sd")),
      !!rlang::sym(paste0(condition2,"_polya_sd")),
      fold_change,
      p.value
    )
  pvals_BH_summary <-
    p.adjust(polyA_data_stat_summary$p.value, method = "BH")
  polyA_data_stat_summary$padj <- pvals_BH_summary

  polyA_data_stat_summary <- polyA_data_stat_summary %>% dplyr::mutate(significance = ifelse(padj < alpha, paste0("FDR<", alpha), "NotSig"))

  polyA_data_stat_summary_short <- polyA_data_stat_summary %>% dplyr::select(transcript,dplyr::ends_with("counts"),dplyr::ends_with("gm_mean"),p.value,padj)


  return(
    list(stats = polya_data_stat, summary = polyA_data_stat_summary,summary_short = polyA_data_stat_summary_short)
  )
}



#' Summarizes input polya table
#'
#' Summarizes input table with polyA predictions, calculating medians, mean, geometric means and standard deviation values for each transcript (default).
#'
#' @param polya_data input table with polyA predictions
#' @param summary_factors specifies column used for grouping (default: group)
#' @param transcript_id_column specifies which column use as transcript identifier (default: transcript)
#'
#' @return long-format \link[tibble]{tibble} with per-transcript statistics for each sample
#' @export
#'
summarize_polya <- function(polya_data,summary_factors = c("group"),transcript_id_column = c("transcript")) {

  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data.frame provided as an input")
  assertthat::assert_that(assertive::is_character(summary_factors),msg = "Non-character argument is not alowed for `summary factors`. Please provide either string or vector of strings")
  assertthat::assert_that(all(summary_factors %in% colnames(polya_data)),msg="Non-existent column name provided as the argument (summary_factors)")

  polya_data_summarized <-
    polya_data %>% dplyr::ungroup() %>% dplyr::group_by(.dots = c(transcript_id_column,summary_factors)) %>% dplyr::summarise(
      counts = dplyr::n(),
      polya_mean = mean(polya_length),
      polya_sd = sd(polya_length),
      polya_median = median(polya_length),
      polya_gm_mean = gm_mean(polya_length)
    )
  return(polya_data_summarized)
}




#' Calculates PCA using polya predictions or counts
#'
#' @param polya_data_summarized summarized polyA predictions. Generate use \link{summarize_polya}
#' @param parameter - parameter used for PCA calculation. One of: polya_median,polya_mean,polya_gm_mean,counts
#'
#' @return pca object
#' @export
#'
calculate_pca <- function(polya_data_summarized,parameter="polya_median") {


  assertthat::assert_that(parameter %in% colnames(polya_data_summarized),msg=paste0(parameter," is not a column of input dataset"))

  polya_data_summarized <- polya_data_summarized %>% dplyr::select(transcript,sample_name,!!rlang::sym(parameter)) %>% tidyr::spread(sample_name,!!rlang::sym(parameter)) %>% as.data.frame()
  polya_data_summarized[is.na(polya_data_summarized)] <- 0
  sample_names <- colnames(polya_data_summarized[,-1])
  transcript_names <- polya_data_summarized[,1]
  polya_data_summarized_t<-t(polya_data_summarized[,-1])
  colnames(polya_data_summarized_t) <- transcript_names
  pca.test <- prcomp(polya_data_summarized_t,center=T,scale=T)

  return_list <- list(pca = pca.test,sample_names = sample_names)
  return(return_list)
}





#' Get information about nanopolish processing
#'
#' Process the information returned by \code{nanopolish polya} in the \code{qc_tag} column
#'
#' @param polya_data A data.frame or tibble containig unfiltered polya output from Nanopolish,
#' @param grouping_factor How to group results (e.g. by sample_name)
#' best read with \link[nanotail]{read_polya_single} or \link[nanotail]{read_polya_multiple}
#' @return A \link[table]{tibble} with counts for each processing state

#' @export
#'
get_nanopolish_processing_info <- function(polya_data,grouping_factor=NA) {


  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data.frame provided as an input")


  if(!is.na(grouping_factor)) {
    assertthat::assert_that(grouping_factor %in% colnames(polya_data),msg=paste0(grouping_factor," is not a column of input dataset"))
    processing_info <- polya_data %>% dplyr::mutate(qc_tag=forcats::fct_relevel(qc_tag,"PASS", after = Inf)) %>% dplyr::group_by(!!rlang::sym(grouping_factor),qc_tag) %>% dplyr::count()
  }
  else {
    processing_info <- polya_data %>% dplyr::mutate(qc_tag=forcats::fct_relevel(qc_tag,"PASS", after = Inf)) %>% dplyr::group_by(qc_tag) %>% dplyr::count()
  }

  return (processing_info)
}



#' Performs differential expression analysis
#'
#' Uses counts for each identified transcript to calculate differential expression between specified groups.
#' This function is a wrapper for \code{\link[edgeR]{binomTest}} from \code{edgeR} package
#'
#'
#' @param polya_data polya_data tibble
#' @param grouping_factor name of column containing factor with groups for comparison
#' @param condition1 first condition to compare
#' @param condition2 second condition to compare
#' @param alpha threshold for a pvalue, to treat the result as significant (default = 0.05)
#' @param summarized_input is input table already summarized?
#'
#' @return a tibble with differential expression results
#' @export
#'
#' @seealso \link[edgeR]{binomTest}
#'
calculate_diff_exp_binom <- function(polya_data,grouping_factor=NA,condition1=NA,condition2=NA,alpha=0.05,summarized_input=FALSE) {



  if (missing(polya_data)) {
    stop("PolyA data are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  if (missing(grouping_factor)) {
    stop("Grouping factor is missing. Please specify one",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::is_a_bool(summarized_input),msg="Non-boolean value provided for option summarized_input")
  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data.frame provided as an input")
  if (!is.na(grouping_factor)) {
    assertthat::assert_that(grouping_factor %in% colnames(polya_data),msg=paste0(grouping_factor," is not a column of input dataset"))
    if (!is.na(condition1)) {
      if(!is.na(condition2)) {
        assertthat::assert_that(condition1 %in% levels(polya_data[[grouping_factor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (grouping_factor)"))
        assertthat::assert_that(condition2 %in% levels(polya_data[[grouping_factor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (grouping_factor)"))
        assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      }
      else {
        stop("Please provide both condition1 and condition2 for comparison")
      }
    }
    else {
      stop("Please provide both condition1 and condition2 for comparison")
    }
  }
  else {
    stop("Please provide valid grouping_factor")
  }

  assertthat::assert_that(is.numeric(alpha),msg = "Non-numeric parameter provided (alpha)")

  if (summarized_input){
    polya_data_summarized <- polya_data
  }
  else{
    polya_data_summarized <- summarize_polya(polya_data,summary_factors = grouping_factor)
  }

  libsizes <- polya_data_summarized %>% dplyr::ungroup() %>% dplyr::group_by(!!rlang::sym(grouping_factor)) %>% dplyr::summarize(lib_size=sum(counts)) %>% dplyr::mutate(sizeFactor=lib_size/mean(lib_size))

  #sum counts for all samples in given group
  polyA_data_counts_summarized<-polya_data_summarized %>% dplyr::group_by(transcript,!!rlang::sym(grouping_factor)) %>% dplyr::summarize(counts_sum=sum(counts)) %>% tidyr::spread(!!rlang::sym(grouping_factor),counts_sum) %>% as.data.frame()

  polyA_data_counts_summarized[is.na(polyA_data_counts_summarized)] <- 0

  polyA_data_counts_summarized <- cbind(polyA_data_counts_summarized$transcript,polyA_data_counts_summarized[,-1]/libsizes$sizeFactor)
  colnames(polyA_data_counts_summarized)[1] <- "transcript"


  binom_test_pvalues<-edgeR::binomTest(y1=polyA_data_counts_summarized[[condition1]],y2=polyA_data_counts_summarized[[condition2]],n1=sum(polyA_data_counts_summarized[[condition1]]),n2=sum(polyA_data_counts_summarized[[condition2]]))
  binom_test_adjusted_pvalues <- p.adjust(binom_test_pvalues,method="BH")

  binom_test_results <-
    data.frame(transcript = polyA_data_counts_summarized$transcript,
               pvalue = binom_test_pvalues,
               padj = binom_test_adjusted_pvalues) %>%
    dplyr::left_join(polyA_data_counts_summarized) %>%
    dplyr::mutate(fold_change = (!! rlang::sym(condition2))/(!! rlang::sym(condition1))) %>%
    #dplyr::mutate_(fold_change = ifelse((condition2 > 0 & condition1 > 0), fold_change, 0)) %>%
    dplyr::mutate(significance = ifelse(padj < alpha, paste0("FDR<", alpha), "NotSig")) %>%
    dplyr::mutate(mean_expr = ((!!rlang::sym(condition1)) + (!!rlang::sym(condition2)))/2) %>%
    dplyr::arrange(padj)

  return(binom_test_results)

}


