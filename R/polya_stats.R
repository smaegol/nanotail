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
#' @param transcript_id_column - name of the column with transcript ids (default = "transcript")
#' @param alpha - alpha value to consider a hit significant (default - 0.05)
#' @param add_summary - add summary (mean polya lengths, counts) to statistics results?
#' @param length_summary_to_show - which length summary to show ("median"/"mean"/"gm_mean")
#' @param ... - additional parameters to pass to .polya_stats (custom_glm_formula,use_dwell_time)
#' @param stat_test what statistical test to use for testing, currently supports "Wilcoxon" (for \link{wilcox.test}), "KS" (for \link[FSA]{ksTest} from FSA package) or "glm" (for \link{glm})
#'
#' @return summary table with pvalues and median/mean values associated to each transcript
#'
#' @export
#'
calculate_polya_stats <- function(polya_data, transcript_id_column = "transcript", min_reads = 0, grouping_factor = "sample_name",condition1=NA,condition2=NA,stat_test="Wilcoxon",alpha=0.05,add_summary=TRUE,length_summary_to_show = "gm_mean",...)
{


  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }



  assertthat::assert_that(assertive::has_rows(polya_data),msg = "Empty data.frame provided as an input")





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

  polya_data_stat <-
    polya_data %>% dplyr::mutate(transcript2=transcript) %>% dplyr::group_by(.dots = c(transcript_id_column)) %>% tidyr::nest()

  #future::plan(future::multiprocess())
  polya_data_stat <- polya_data_stat %>% dplyr::mutate(stats=purrr::map(data,.polya_stats,grouping_factor=grouping_factor,stat_test=stat_test,min_reads=min_reads,...)) %>% dplyr::select(-data) %>% tidyr::unnest()
  message("calculating statistics")
  #polya_data_stat <- polya_data_stat %>% dplyr::mutate(stats=furrr::future_map(data,.polya_stats,grouping_factor=grouping_factor,stat_test=stat_test,min_reads=min_reads,...)) %>% dplyr::select(-data) %>% tidyr::unnest()
  message("Finished")
  if (add_summary) {
    polyA_data_stat_summary <- summarize_polya(polya_data,summary_factors = grouping_factor,transcript_id_column = transcript_id_column) %>% dplyr::select(!!rlang::sym(transcript_id_column),counts,!!rlang::sym(grouping_factor),!!rlang::sym(paste0("polya_",length_summary_to_show))) %>% spread_multiple(!!rlang::sym(grouping_factor),c(counts,!!rlang::sym(paste0("polya_",length_summary_to_show))))
    polya_data_stat <- polya_data_stat %>% dplyr::full_join(polyA_data_stat_summary,by=transcript_id_column)
    polya_data_stat$length_diff <- polya_data_stat[[paste0(condition2,"_polya_",length_summary_to_show)]] - polya_data_stat[[paste0(condition1,"_polya_",length_summary_to_show)]]
    polya_data_stat$fold_change <- polya_data_stat[[paste0(condition2,"_polya_",length_summary_to_show)]] / polya_data_stat[[paste0(condition1,"_polya_",length_summary_to_show)]]
  }
  message("Adjusting p.value")
  polya_data_stat$padj <- p.adjust(polya_data_stat$p.value, method = "BH")
  polya_data_stat<- polya_data_stat %>% dplyr::mutate(effect_size=dplyr::case_when((abs(cohen_d))<0.2 ~ "negligible",(abs(cohen_d)<0.5) ~ "small", (abs(cohen_d)<0.8) ~ "medium",(abs(cohen_d)>=0.8) ~ "large",TRUE ~ "NA"))
  
  # create significance factor
  polya_data_stat <-
    polya_data_stat %>% dplyr::mutate(significance = dplyr::case_when(is.na(padj)  ~ "NotSig",
                                                                      (padj < alpha) ~ paste0("FDR<", alpha),
                                                                      TRUE ~ "NotSig"))

  polya_data_stat$stats_code <- sapply(polya_data_stat$stats_code,FUN = function(x) {stat_codes_list[[x]]},simplify = "vector",USE.NAMES = FALSE) %>% unlist()



  polya_data_stat <- polya_data_stat %>% dplyr::arrange(padj)
  #polya_data_stat_short <- polya_data_stat %>% dplyr::select(!!rlang::sym(transcript_id_column),dplyr::ends_with("counts"),dplyr::ends_with("gm_mean"),p.value,padj,stats_code)


  #return(
  # list(summary = polya_data_stat,summary_short = polya_data_stat_short))
  return(polya_data_stat)
}



stat_codes_list = list(OK = "OK",
                   G1_NA = "GROUP1_NA",
                   G2_NA = "GROUP2_NA",
                   G1_LC = "G1_LOW_COUNT",
                   G2_LC = "G2_LOW_COUNT",
                   B_NA = "DATA FOR BOTH GROUPS NOT AVAILABLE",
                   B_LC = "LOW COUNTS FOR BOTH GROUPS",
                   G_LC = "LOW COUNT FOR ONE GROUP",
                   G_NA = "DATA FOR ONE GROUP NOT AVAILABLE",
                   ERR = "OTHER ERROR",
                   GLM_GROUP = "MISSING FACTOR LEVEL IN GLM CALL")

#' Calculates polyA statistics for single group of reads (for single transcript)
#'
#' @param polya_data - input data frame with polyA predictions
#' @param stat_test - statistical test to use. One of : Wilcoxon, KS (Kolmogorov-Smirnov) or glm (Generalized Linear Model). All tests use log2(polya_length) as a response variable
#' @param grouping_factor - factor defining groups (Need to have 2 levels)
#' @param min_reads - minimum reads per group to include in the statistics calculation
#' @param use_dwell_time - use dwell time instead of calculated polya length for statistics
#' @param custom_glm_formula - custom glm formula (when using glm for statistics)
#'
#' @return data frame
#'
.polya_stats <- function(polya_data,stat_test,grouping_factor,min_reads=0,use_dwell_time = FALSE,custom_glm_formula = NA) {


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
      polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(grouping_factor) %in% c(condition1,condition2)) %>% droplevels()
    }
  }
  else if (length(levels(polya_data[[grouping_factor]]))==1) {
    stop("Only 1 level present for grouping factor. Choose another groping factor for comparison")
  }
  else {
    condition1 = levels(polya_data[[grouping_factor]])[1]
    condition2 = levels(polya_data[[grouping_factor]])[2]
  }


  # initial status code
  stats_code = codes_stats = "OK"
  # calculate group counts
  group_counts = polya_data  %>% dplyr::group_by(.dots = c(grouping_factor)) %>% dplyr::count()

  stats <- NA
  cohend <- NA
  if (use_dwell_time) {
    statistics_formula <- paste0("dwell_time ~",grouping_factor)
  }
  else {
    statistics_formula <- paste0("log2(polya_length) ~",grouping_factor)
  }

  if ((!missing(custom_glm_formula)) & (stat_test!='glm')) {
    warning("custom_glm_formula specified but glm is not used for statistics calculation. Formula will be ignored")
  }

 
  if (nrow(group_counts)==2) {
    if (group_counts[1,]$n < min_reads) {
      if (group_counts[2,]$n < min_reads) {
        #message("Not enough counts for both groups")
        stats_code = "B_LC"
      }
      else {
        #message(paste0("Not enough counts for group ",group_counts[1,]$group))
        stats_code = "G_LC"
      }
    }
    else if (group_counts[2,]$n < min_reads) {
      #message(paste0("Not enough counts for group ",group_counts[2,]$group))
      stats_code = "G_LC"
    }
    else {
      #calculate cohen's d parameter
      cohend<-effsize::cohen.d(data=polya_data,as.formula(statistics_formula))$estimate
      if (stat_test=="Wilcoxon") {
        #print(polya_data$transcript2)
        stats <- suppressWarnings(wilcox.test(as.formula(statistics_formula),polya_data))$p.value
      }
      else if (stat_test=="KS") {
        stats <- FSA::ksTest(as.formula(statistics_formula),data=polya_data)$p.value
      }
      else if (stat_test=="glm") {
        valid_glm_groups = TRUE
        if(!missing(custom_glm_formula)) {
          custom_glm_formula <- substitute(custom_glm_formula)
          glm_groups<-all.vars(as.formula(custom_glm_formula))[-1]
          assertthat::assert_that(length(glm_groups)>0,msg="Please provide valid custom_glm_formula (with the correct categorical explanatory variable)")
          assertthat::assert_that(all(glm_groups %in% colnames(polya_data)),msg = "wrong custom glm formula. Please provide valid column names as explanatory variable")
          # check that at least 2 levels present for each explanatory variable
          for (glm_group in glm_groups) {
            if (length(unique(polya_data[[glm_group]]))<2) {
              valid_glm_groups = FALSE
            }
          }
          statistics_formula = custom_glm_formula
        }

        if (valid_glm_groups) {
          # required, as log2(0) produces inf, throwing error in glm
          polya_data <- polya_data %>% dplyr::mutate(polya_length = ifelse(polya_length==0,1,polya_length))
          mcp_call <- paste0("multcomp::mcp(",grouping_factor,' = "Tukey")')
          stats <- summary(multcomp::glht(glm(formula = as.formula(statistics_formula),data=polya_data),eval(parse(text = mcp_call))))$test$pvalues[1]
          
          #polya_data_stat <- polya_data_stat  %>% dplyr::mutate(stats = suppressWarnings(coef(summary(glm(formula = as.formula(statistics_formula))))[2,4]))
        }
        else{
          stats <- NA
          stats_code <- "GLM_GROUP"
        }
      }
      else {
        stop("wrong stat_test parameter provided")
      }
    }
  }
  else if (nrow(group_counts)==1) {
    stats_code = "G_NA"
  }
  else if (nrow(group_counts)==0) {
    stats_code = "B_NA"
  }
  else {
    stats_code = "ERR"
  }


  stats<-tibble::tibble(p.value=stats,stats_code=as.character(stats_code),cohen_d=cohend)

  return(stats)

}


#' Summarizes input polya table
#'
#' Summarizes input table with polyA predictions, calculating medians, mean, geometric means and standard deviation values for each transcript (default).
#' To get overall summary for each sample or group, specify `transcript_id_column=NULL`
#'
#' @param polya_data input table with polyA predictions
#' @param summary_factors specifies column used for grouping (default: group)
#' @param transcript_id_column specifies which column use as transcript identifier (default: transcript). Set to `NULL` to omit per-transcript stats
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
      polya_gm_mean = gm_mean(polya_length),
      polya_sem = polya_sd/sqrt(counts)
    )
  return(polya_data_summarized)
}




#' Calculates PCA using polya predictions or counts
#'
#' Needs polyA predictions table summarized by \link{summarize_polya} function, using "sample_name" as summary_factors
#'
#' @param polya_data_summarized summarized polyA predictions. Generate use \link{summarize_polya}
#' @param parameter - parameter used for PCA calculation. One of: polya_median,polya_mean,polya_gm_mean,counts
#' @param transcript_id_column column which respresnrt transcript id
#'
#' @return pca object
#' @export
#'
calculate_pca <- function(polya_data_summarized,parameter="polya_median",transcript_id_column = "transcript") {


  if (missing(polya_data_summarized)) {
    stop("Summarized PolyA predictions are missing. Please provide a valid polya_data_summarized argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(polya_data_summarized),msg = "Empty data.frame provided as an input")
  assertthat::assert_that(assertive::is_character(parameter),msg = "Non-character argument is not alowed for `parameter`.")
  assertthat::assert_that(transcript_id_column %in% colnames(polya_data_summarized),msg=paste0("Required `transcript`` column is missing from input dataset."))
  assertthat::assert_that("sample_name" %in% colnames(polya_data_summarized),msg=paste0("Required `sample_name` column is missing from input dataset."))
  assertthat::assert_that(parameter %in% colnames(polya_data_summarized),msg=paste0(parameter," is not a column of input dataset."))

  polya_data_summarized <- polya_data_summarized %>% dplyr::select(!!rlang::sym(transcript_id_column),sample_name,!!rlang::sym(parameter)) %>% tidyr::spread(sample_name,!!rlang::sym(parameter)) %>% as.data.frame()
  polya_data_summarized[is.na(polya_data_summarized)] <- 0
  sample_names <- colnames(polya_data_summarized[,-1])
  transcript_names <- polya_data_summarized[,1]
  polya_data_summarized_t<-t(polya_data_summarized[,-1])
  # step required to remove zero-variance columns (based on https://stackoverflow.com/questions/40315227/how-to-solve-prcomp-default-cannot-rescale-a-constant-zero-column-to-unit-var)
  zero_variance_transcripts<-which(apply(polya_data_summarized_t, 2, var)==0)
  if (length(zero_variance_transcripts>0)) {
    polya_data_summarized_t <- polya_data_summarized_t[ , -zero_variance_transcripts]
    colnames(polya_data_summarized_t) <- transcript_names[-zero_variance_transcripts]
  }
  else {
    colnames(polya_data_summarized_t) <- transcript_names
  }
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
#' @return A \link[tibble]{tibble} with counts for each processing state

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


