#' PCA biplot
#'
#' Plots PCA biplot using \link[ggbiplot]{ggbiplot}
#'
#' @param pca_object pca object
#' @param samples_names names of samples to be shown on the plot
#'
#' @return \link[ggplot2]{ggplot} object

#' @export
#'
plot_polyA_PCA <- function(pca_object,samples_names) {

  ggbiplot::ggbiplot(pca_object,var.axes = F,labels = samples_names)

}


#' Plots polya distribution as ndensity plot
#'
#' @param polya_data input data with polyA predictions
#' @param groupingFactor how to group
#' @param parameter_to_plot what to plot on x scale (defaults to polya_length)
#' @param condition1 First condition to include on the plot
#' @param condition2 Second condition to include on the plot
#' @param show_center_values Show center values as vertical line. Possible values: "none","median","mean"
#' @param subsample Subsample input table, provide either absolute number or fraction
#' @param ndensity Should ndensity (scaled density) be plotted instead of normal denisty (default = TRUE)
#' @param mode_method method used for the mode calculation (argument to modeest::mlv method parameter)
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_polya_distribution <- function(polya_data, groupingFactor=NA, parameter_to_plot = "polya_length", condition1=NA,condition2=NA,show_center_values="none",subsample=NA,ndensity=TRUE,mode_method="density",...) {


  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }


  assertthat::assert_that(show_center_values %in% c("none","median","mean","gm_mean","mode"))



  if (!is.na(groupingFactor)) {
    assertthat::assert_that(groupingFactor %in% colnames(polya_data),msg=paste0(groupingFactor," is not a column of input dataset"))
    if (!is.na(condition1)) {
      if(!is.na(condition2)) {
        assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",groupingFactor," (groupingFactor)"))
        assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",groupingFactor," (groupingFactor)"))
        assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
        polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
      }
    }
    if(!is.na(subsample)) {
      polya_data <- subsample_table(polya_data,groupingFactor = groupingFactor,subsample=subsample)
    }
    distribution_plot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=parameter_to_plot,color=groupingFactor))
    if (ndensity) {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::ylab("normalized density")
    }
    else {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..density..)) + ggplot2::ylab("density")
    }
    distribution_plot <- distribution_plot + ggplot2::theme_bw()

    if(show_center_values!="none") {
      center_values = polya_data %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::summarize(median_value = median(!!rlang::sym(parameter_to_plot),na.rm = TRUE),mean_value=mean(!!rlang::sym(parameter_to_plot),na.rm=TRUE),gm_mean_value=gm_mean(!!rlang::sym(parameter_to_plot),na.rm=TRUE),mode_value=modeest::mlv(!!rlang::sym(parameter_to_plot),method=mode_method,na.rm=TRUE))
      if(show_center_values=='median') {
            distribution_plot <- distribution_plot + ggplot2::geom_vline(data=center_values,ggplot2::aes_string(xintercept="median_value",color=groupingFactor),linetype="longdash")
      }
      else if(show_center_values=='mean') {
        distribution_plot <- distribution_plot + ggplot2::geom_vline(data=center_values,ggplot2::aes_string(xintercept="mean_value",color=groupingFactor),linetype="longdash")
      }
      else if(show_center_values=='gm_mean') {
        distribution_plot <- distribution_plot + ggplot2::geom_vline(data=center_values,ggplot2::aes_string(xintercept="gm_mean_value",color=groupingFactor),linetype="longdash")
      }
      else if(show_center_values=='mode') {
        distribution_plot <- distribution_plot + ggplot2::geom_vline(data=center_values,ggplot2::aes_string(xintercept="mode_value",color=groupingFactor),linetype="longdash")
      }
    }
  }
  else {
    if(!is.na(subsample)) {
      polya_data <- subsample_table(polya_data,subsample=subsample)
    }
    distribution_plot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=parameter_to_plot))
    if (ndensity) {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::ylab("normalized density")
    }
    else {
      distribution_plot <- distribution_plot + ggplot2::geom_line(stat="density",size=1,ggplot2::aes(y=..density..)) + ggplot2::ylab("density")
    }

    distribution_plot <- distribution_plot  + ggplot2::theme_bw()
    if(show_center_values=="median") {
    }
    else if (show_center_values=="mean") {
      distribution_plot <- distribution_plot + ggplot2::geom_vline(aes(xintercept=mean(polya_length)),linetype="longdash")
    }
  }

  distribution_plot <- distribution_plot + ggplot2::xlab("poly(A) length")

  distribution_plot <- .basic_aesthetics(distribution_plot,...)




  distribution_plot <- distribution_plot + nanotail_ggplot2_theme

  return(distribution_plot)

}




#' Plots boxplot of estimated polya lengths
#'
#' @param polya_data input table with polyA predictions
#' @param groupingFactor which factor to use for grouping
#' @param condition1 First condition to include on the plot
#' @param condition2 Second condition to include on the plot
#' @param violin Should violin plot be plotted instead of boxplot?
#' @param additional_grouping_factor additional coloring grouping factor
#' @param add_points should individual points be plotted (only if less than max_points). Represented as \link[ggforce]{geom_sina}
#' @param max_points maximum number of points to be plotted if add_points is specified
#'
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)

#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_polya_boxplot <- function(polya_data, groupingFactor,additional_grouping_factor=NA,condition1=NA,condition2=NA,violin=FALSE,add_points=FALSE,max_points=500,...) {


  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }

  assertthat::assert_that(groupingFactor %in% colnames(polya_data),msg=paste0(groupingFactor," is not a column of input dataset"))
  if(!is.na(additional_grouping_factor)) {
    assertthat::assert_that(additional_grouping_factor %in% colnames(polya_data),msg=paste0(additional_grouping_factor," is not a column of input dataset"))
  }

  if (!is.na(condition1)) {
    if(!is.na(condition2)) {
      assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",groupingFactor," (groupingFactor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",groupingFactor," (groupingFactor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
    }
  }

  if (!is.na(additional_grouping_factor)) {
    transcripts_boxplot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=groupingFactor,y="polya_length",color=additional_grouping_factor))
  }
  else {
    transcripts_boxplot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=groupingFactor,y="polya_length"))
  }
  if(violin) {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_violin(position=ggplot2::position_dodge())
  }
  else{
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_boxplot(position=ggplot2::position_dodge())
  }
  if(add_points) {
    points_counts<-polya_data %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::count()
    if (!any(points_counts$n>max_points)) {
      transcripts_boxplot <- transcripts_boxplot + ggforce::geom_sina(position=ggplot2::position_dodge())
    }
  }

  transcripts_boxplot <- .basic_aesthetics(transcripts_boxplot,...)

  return(transcripts_boxplot)
}





#' Plots violin plot of estimated polya lengths, other version
#'
#' @param polya_data input table with polyA predictions
#' @param groupingFactor which factor to use for grouping
#' @param condition1 First condition to include on the plot
#' @param condition2 Second condition to include on the plot
#' @param violin Should violin plot be plotted instead of boxplot?
#' @param additional_grouping_factor additional coloring grouping factor
#' @param add_points should individual points be plotted (only if less than max_points). Represented as \link[ggforce]{geom_sina}
#' @param max_points maximum number of points to be plotted if add_points is specified
#' @param add_boxplot Add boxplot inside violin?
#'
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)

#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_polya_violin <- function(polya_data, groupingFactor,additional_grouping_factor=NA,condition1=NA,condition2=NA,violin=FALSE,add_points=FALSE,max_points=500,add_boxplot=TRUE,fill_by=NA,...) {
  
  
  if (missing(polya_data)) {
    stop("PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }
  
  assertthat::assert_that(groupingFactor %in% colnames(polya_data),msg=paste0(groupingFactor," is not a column of input dataset"))
  if(!is.na(additional_grouping_factor)) {
    assertthat::assert_that(additional_grouping_factor %in% colnames(polya_data),msg=paste0(additional_grouping_factor," is not a column of input dataset"))
  }
  
  if(!is.na(fill_by)) {
    assertthat::assert_that(fill_by %in% colnames(polya_data),msg=paste0(additional_grouping_factor," is not a column of input dataset"))
  }
  else{
    fill_by = groupingFactor
  }
  
  if (!is.na(condition1)) {
    if(!is.na(condition2)) {
      assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",groupingFactor," (groupingFactor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",groupingFactor," (groupingFactor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      polya_data <- polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
    }
  }
  
  polya_data <- polya_data %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::add_count() %>% dplyr::ungroup() %>% dplyr::mutate(label = paste0("(n=", n, ")"))
  
 
  transcripts_boxplot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x="label",y="polya_length",fill=fill_by)) + geom_violin(trim = FALSE)
  
  if(add_points) {
    points_counts<-polya_data %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::count()
    if (!any(points_counts$n>max_points)) {
      transcripts_boxplot <- transcripts_boxplot + ggforce::geom_sina(position=ggplot2::position_dodge())
    }
  }
  if(add_boxplot) {
  transcripts_boxplot <- transcripts_boxplot + geom_boxplot(
    width = 0.2,
    fill = "white",
    notch = TRUE,
    outlier.shape = NA
  )  
  }
  
  transcripts_boxplot <- transcripts_boxplot + 
    facet_grid(. ~ group, drop = T, scales = "free_x") + 
    nanotail::stat_median_line(color = "red", linetype = "dashed") + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0)) + 
    scale_x_discrete() + 
    xlab("") 
  
  transcripts_boxplot <- .basic_aesthetics(transcripts_boxplot,...)
  
  return(transcripts_boxplot)
}

#' Title
#'
#' @param polya_data_summarized polyA predictions table, summarized using \link{summarize_polya}
#' @param groupingFactor name of column used for grouping
#' @param condition1 first condition to use for plotting
#' @param condition2 second condition to use for plotting
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#' @param min_counts minimum number of counts to be shown
#' @param max_counts maximum number of counts to be shown
#' @param points_coloring_factor factor specifying how to color points
#' @param repel_elements TBD
#' @param repel_group TBD
#' @param transcript_id_column TBD
#' #'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_counts_scatter <- function(polya_data_summarized, groupingFactor = NA, condition1 = NA, condition2 = NA,min_counts = 0, max_counts = 0,points_coloring_factor =NA, repel_elements=NA,repel_group=NA,transcript_id_column="transcript",...) {


  if (missing(polya_data_summarized)) {
    stop("Summarized PolyA predictions are missing. Please provide a valid polya_data argument",
         call. = FALSE)
  }


  assertthat::assert_that(!is.na(condition1),msg = "Please specify conditions for comparison")
  assertthat::assert_that(!is.na(condition2),msg = "Please specify conditions for comparison")
  assertthat::assert_that(!is.na(groupingFactor),msg = "Please specify groupingFactor")
  assertthat::assert_that(condition1 %in% levels(polya_data_summarized[[groupingFactor]]),msg=paste0(condition1," is not a level of ",groupingFactor," (groupingFactor)"))
  assertthat::assert_that(condition2 %in% levels(polya_data_summarized[[groupingFactor]]),msg=paste0(condition2," is not a level of ",groupingFactor," (groupingFactor)"))
  assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
  assertthat::assert_that("counts" %in% colnames(polya_data_summarized),msg = "Please provide summarized polya table (using summarize_polya()) as an input")



  if(!is.na(points_coloring_factor)){
    assertthat::assert_that(points_coloring_factor %in% colnames(polya_data_summarized),msg = "Please provide valid points_coloring_factor as an input")
    polya_data_summarized_counts_xy<-polya_data_summarized %>% dplyr::group_by(!!rlang::sym(transcript_id_column),!!rlang::sym(groupingFactor),!!rlang::sym(points_coloring_factor)) %>% dplyr::summarize(counts_sum=sum(counts)) %>% tidyr::spread_(groupingFactor,"counts_sum")
  }
  else {
    polya_data_summarized_counts_xy<-polya_data_summarized %>% dplyr::group_by(!!rlang::sym(transcript_id_column),!!rlang::sym(groupingFactor)) %>% dplyr::summarize(counts_sum=sum(counts)) %>% tidyr::spread_(groupingFactor,"counts_sum")
  }
  polya_data_summarized_counts_xy[is.na(polya_data_summarized_counts_xy)] <- 0

  polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)>=min_counts,!!rlang::sym(condition2)>=min_counts,!!rlang::sym(condition1)<=max_counts,!!rlang::sym(condition2)<=max_counts)

  if (!is.na(max_counts)) {
    assertthat::assert_that(assertive::is_numeric(max_counts),msg="Please provide numeric value for max_counts")
    polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)<=max_counts,!!rlang::sym(condition2)<=max_counts)
  }

  if(!is.na(points_coloring_factor)){
    counts_scatter_plot<-ggplot2::ggplot(polya_data_summarized_counts_xy,ggplot2::aes(x=!!rlang::sym(condition1),y=!!rlang::sym(condition2),colour=!!rlang::sym(points_coloring_factor))) + ggplot2::geom_point(ggplot2::aes_string(text=transcript_id_column),alpha=0.7)
  }
  else
    {
    counts_scatter_plot<-ggplot2::ggplot(polya_data_summarized_counts_xy,ggplot2::aes(x=!!rlang::sym(condition1),y=!!rlang::sym(condition2))) + ggplot2::geom_point(ggplot2::aes_string(text=transcript_id_column),alpha=0.7)
  }

  if (!is.na(repel_elements)) {
    assertthat::assert_that(assertive::is_numeric(repel_elements),msg="Please provide numeric paraemter for repel_elements")
    if(!is.na(repel_group)) {
    assertthat::assert_that(assertive::is_numeric(repel_elements),msg="Please provide numeric paraemter for repel_elements")
    counts_scatter_plot <- counts_scatter_plot + ggrepel::geom_text_repel(data=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::filter(!!rlang::sym(points_coloring_factor) == repel_group) %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1)))[1:repel_elements,], ggplot2::aes(label=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::filter(!!rlang::sym(points_coloring_factor) == repel_group) %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1))) %>% dplyr::select(!!rlang::sym(transcript_id_column)) %>% as.vector()[1:repel_elements]))
    }
    else {
      counts_scatter_plot <- counts_scatter_plot + ggrepel::geom_text_repel(data=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1))) %>% dplyr::slice(1:repel_elements), ggplot2::aes(label=polya_data_summarized_counts_xy %>% dplyr::ungroup() %>% dplyr::arrange(dplyr::desc(!!rlang::sym(condition1))) %>% dplyr::slice(1:repel_elements) %>% dplyr::select(!!rlang::sym(transcript_id_column)) %>% as.vector()))
    }
  }
  counts_scatter_plot <- .basic_aesthetics(counts_scatter_plot,...)


  return(counts_scatter_plot)
}



#' Plot Nanopolish polya QC
#'
#' @param nanopolish_processing_info output of \link{get_nanopolish_processing_info}
#' @param frequency show frequency plot instead of counts plot
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_nanopolish_qc <- function(nanopolish_processing_info, frequency=TRUE,...) {


  if (missing(nanopolish_processing_info)) {
    stop("nanopolish processing info is missing. Please provide a valid nanopolish_processing_info argument",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(nanopolish_processing_info),msg = "Empty data.frame provided as an input")
  basic_colnames = c("qc_tag","n")
  assertthat::assert_that(basic_colnames[1] %in% colnames(nanopolish_processing_info),msg="qc_tag column is missing in the input. Is that valid output of get_nanopolish_processing_info()?")
  assertthat::assert_that(basic_colnames[2] %in% colnames(nanopolish_processing_info),msg="n column is missing in the input. Is that valid output of get_nanopolish_processing_info()?")

  assertthat::assert_that(assertive::is_a_bool(frequency),msg="Non-boolean value provided for option frequency")


  # if there were multiple samples compared
  if (ncol(nanopolish_processing_info)>2) {
    grouping_colname = setdiff(colnames(nanopolish_processing_info),basic_colnames)
    nanopolish_qc_plot <- ggplot2::ggplot(nanopolish_processing_info,ggplot2::aes(x=!!rlang::sym(grouping_colname),fill=qc_tag,y=n))
    if(frequency) {
      nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::geom_bar(stat="identity",position="fill") + ggplot2::ylab("frequency")
    }
    else{
      nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::geom_bar(stat="identity",position="stack") + ggplot2::ylab("count")
    }
    nanopolish_qc_plot <- nanopolish_qc_plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,vjust = 0.7))
  }
  else {
    nanopolish_qc_plot <- ggplot2::ggplot(nanopolish_processing_info,ggplot2::aes(x=qc_tag,y=n)) + ggplot2::geom_bar(stat="identity")
  }

  nanopolish_qc_plot <- .basic_aesthetics(nanopolish_qc_plot,color_mode = "fill",...)
  return(nanopolish_qc_plot)
}



#' Plots volcano plot of differential expression analysis
#'

#' @param input_data a table with output from \link{calculate_diff_exp_binom} or \link{calculate_polya_stats}
#'
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#' @param transcript_id_column column used for transcript id
#' @param labels show point labels using ggrepel
#' @param nlabels number of labels to show
#' #'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_volcano <- function(input_data,transcript_id_column,labels=FALSE,nlabels=10,...) {


  if (missing(input_data)) {
    stop("input data missing. Please provide a valid input",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(input_data),msg = "Empty data.frame provided as an input")
  assertthat::assert_that("fold_change" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). fold_change column is missing.")
  assertthat::assert_that("padj" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). padj column is missing.")
  assertthat::assert_that("significance" %in% colnames(input_data),msg = "Input table is not a valid input for plot_volcano(). significance column is missing.")


  input_data <- input_data %>% dplyr::mutate(padj=as.numeric(padj)) %>% dplyr::filter(!is.na(padj))

  volcano_plot <- ggplot2::ggplot(input_data,ggplot2::aes(x=log2(as.numeric(fold_change)),y=-log10(as.numeric(padj)),col=significance))
  if(!missing(transcript_id_column)) {
    volcano_plot <- volcano_plot + ggplot2::geom_point(ggplot2::aes(text=!!rlang::sym(transcript_id_column)))
  }
  else {
    volcano_plot <- volcano_plot + ggplot2::geom_point()
  }

  if (labels) {
    if(!missing(transcript_id_column)) {
      labels_df <- input_data %>% dplyr::arrange(padj) %>% dplyr::filter(grepl("^FDR",significance))
      labels_df <- head(labels_df,nlabels)
      volcano_plot <- volcano_plot + ggrepel::geom_text_repel(data=labels_df,ggplot2::aes(label=!!rlang::sym(transcript_id_column)),colour="black")
    }
    else {
      warning("transcript_id_column is missing. Labels will not be produced")
    }
  }

  volcano_plot <- .basic_aesthetics(volcano_plot,...)

  return(volcano_plot)

}


#' Plots MA plot of differential expression analysis
#'
#' Crates simple MA plot, with log10(mean expression) on the X-axis and log2(fold_change) on the Y-axis
#'
#'
#' @param input_data a table with output from \link{calculate_diff_exp_binom}
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#' @param transcript_id_column column used for transcript id
#' @param labels show point labels using ggrepel
#' @param nlabels number of labels to show
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'

plot_MA <- function(input_data,transcript_id_column,labels=FALSE,nlabels=10,...) {


  if (missing(input_data)) {
    stop("input data missing. Please provide a valid input",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(input_data),msg = "Empty data.frame provided as an input")
  assertthat::assert_that("fold_change" %in% colnames(input_data),msg = "Input table is not a valid input for plot_MA(). fold_change column is missing.")
  assertthat::assert_that("mean_expr" %in% colnames(input_data),msg = "Input table is not a valid input for plot_MA(). mean_expr column is missing.")
  assertthat::assert_that("significance" %in% colnames(input_data),msg = "Input table is not a valid input for plot_MA(). significance column is missing.")

  input_data <- input_data %>% dplyr::filter(!is.na(fold_change))

  MA_plot <- ggplot2::ggplot(input_data,ggplot2::aes(x=log10(mean_expr),y=log2(fold_change),col=significance))
  if(!missing(transcript_id_column)) {
    MA_plot <- MA_plot + ggplot2::geom_point(ggplot2::aes(text=!!rlang::sym(transcript_id_column)))
  }
  else {
    MA_plot <- MA_plot + ggplot2::geom_point()
  }

  if (labels) {
    if(!missing(transcript_id_column)) {
      labels_df <- input_data %>% dplyr::arrange(desc(mean_expr)) %>% dplyr::filter(grepl("^FDR",significance))
      labels_df <- head(labels_df,nlabels)
      MA_plot <- MA_plot + ggrepel::geom_text_repel(data=labels_df,ggplot2::aes(label=!!rlang::sym(transcript_id_column)),colour="black")
    }
    else {
      warning("transcript_id_column is missing. Labels will not be produced")
    }
  }


  MA_plot <- .basic_aesthetics(MA_plot,...)

  return(MA_plot)

}

#' Title
#'
#' @param annotated_polya_data data frame(or tibble) with polyA predictions and associated annotations
#' @param grouping_factor column in polya_data_table specifing factor grouping samples
#' @param annotation_factor column specifying factor grouping transcripts by annotation
#' @param condition1 if only 2 conditions to show, choose which one is first
#' @param condition2 if only 2 conditions to show, choose which one is second
#' @param annotation_levels vector specifying selected annotation levels from annotation_factor
#' @param ... parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)
#' @param violin plot violin instead of boxplot?
#'
#' @return \link[ggplot2]{ggplot} object
#' @export
#'
plot_annotations_comparison_boxplot <- function(annotated_polya_data,annotation_factor = NA,grouping_factor = NA,condition1=NA,condition2=NA,annotation_levels=c(),violin=FALSE,...) {

  if (missing(annotated_polya_data)) {
    stop("Annotated polya data table is missing. Please provide a valid input",
         call. = FALSE)
  }

  assertthat::assert_that(assertive::has_rows(annotated_polya_data),msg = "Empty data.frame provided as an input")


  if(!is.na(annotation_factor)) {
    assertthat::assert_that(annotation_factor %in% colnames(annotated_polya_data),msg=paste0(annotation_factor," is not a column of input data.frame"))
  }

  if (length(annotation_levels)>0) {
    #assertthat::assert_that(all(annotation_levels %in% levels(annotated_polya_data[[annotation_factor]])),msg="non-existing factor levels specified for annotation_levels parameter")
    annotated_polya_data <- annotated_polya_data %>% dplyr::filter(!!rlang::sym(annotation_factor) %in% annotation_levels)
    print(annotation_levels)
  }


  if (!is.na(condition1)) {
    if(!is.na(condition2)) {
      assertthat::assert_that(condition1 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition1," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 %in% levels(polya_data[[groupingFactor]]),msg=paste0(condition2," is not a level of ",grouping_factor," (groupingFactor)"))
      assertthat::assert_that(condition2 != condition1,msg="condition2 should be different than condition1")
      annotated_polya_data <- annotated_polya_data %>% dplyr::filter(!!rlang::sym(groupingFactor) %in% c(condition1,condition2))
    }
  }

  if (!is.na(grouping_factor)) {
    assertthat::assert_that(grouping_factor %in% colnames(annotated_polya_data))
    transcripts_boxplot <- ggplot2::ggplot(annotated_polya_data,ggplot2::aes_string(x=annotation_factor,y="polya_length",color=grouping_factor))
  }
  else {
    transcripts_boxplot <- ggplot2::ggplot(annotated_polya_data,ggplot2::aes_string(x=annotation_factor,y="polya_length"))
  }
  if(violin) {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_violin(position = ggplot2::position_dodge())
  }
  else{
    transcripts_boxplot <- transcripts_boxplot + ggplot2::geom_boxplot(position=ggplot2::position_dodge())
  }

  transcripts_boxplot <- .basic_aesthetics(transcripts_boxplot,...)
  transcripts_boxplot <- transcripts_boxplot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust = 1, hjust = 1))
  return(transcripts_boxplot)

}


#' Title
#'
#' @param ggplot_object ggplot2 object to manipulate asesthetics
#' @param scale_x_limit_low lower limit of x continuous scale
#' @param scale_x_limit_high upper limit of x continuous scale
#' @param scale_y_limit_low lower limit of y continuous scale
#' @param scale_y_limit_high upper limit of y continuous scale
#' @param color_palette color palette (one from RColorBrewer of ggsci packages)
#' @param plot_title Title of the plot
#' @param color_mode if using color, fill or both, when specifying color_palette
#' @param axis_titles_size size of axis titles
#'
#' @return \link[ggplot2]{ggplot} object
#'
.basic_aesthetics <- function(ggplot_object,scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA,color_mode="color",axis_titles_size=16)
{


  if(missing(ggplot_object)) {
    stop("ggplot object required as the input")
  }

  if (!is.na(scale_x_limit_low)) {
    if (is.na(scale_x_limit_high)) {
      stop("Please provide both limits for x scale")
    }
    assertthat::assert_that(assertive::is_numeric(scale_x_limit_low),msg="Please provide numeric value for scale_x_limit_low")
    assertthat::assert_that(assertive::is_numeric(scale_x_limit_high),msg="Please provide numeric value for scale_x_limit_high")
    ggplot_object <- ggplot_object + ggplot2::scale_x_continuous(limits=c(scale_x_limit_low,scale_x_limit_high))
  }


  if (!is.na(scale_y_limit_low)) {
    if (is.na(scale_y_limit_high)) {
      stop("Please provide both limits for y scale")
    }
    assertthat::assert_that(assertive::is_numeric(scale_y_limit_low),msg="Please provide numeric value for scale_y_limit_low")
    assertthat::assert_that(assertive::is_numeric(scale_y_limit_high),msg="Please provide numeric value for scale_y_limit_high")
    ggplot_object <- ggplot_object + ggplot2::scale_y_continuous(limits=c(scale_y_limit_low,scale_y_limit_high))
  }


  valid_color_palettes_ggsci = c("npg", "aaas", "nejm", "lancet", "jama",
                          "jco", "ucscgb", "d3", "locuszoom",
                          "igv", "uchicago", "startrek", "tron",
                          "futurama", "rickandmorty", "simpsons")
  valid_color_palettes_RColorBrewer <- rownames(RColorBrewer::brewer.pal.info)

  if(!is.na(color_palette)) {
    if(!color_palette %in% c(valid_color_palettes_RColorBrewer,valid_color_palettes_ggsci)) {
      warning("Please provide valid color palette from RColorBrewer or ggsci packages")
    }
  }

  if(color_mode %in% c("color","both")) {
    if(color_palette %in% c(valid_color_palettes_ggsci)) {
      ggplot_object <- ggplot_object + eval(parse(text = paste0("ggsci::scale_color_",color_palette,"()")))
    }
    else if (color_palette %in% valid_color_palettes_RColorBrewer) {
      ggplot_object <- ggplot_object + ggplot2::scale_colour_brewer(palette = color_palette)
    }
  }
  else if (color_mode %in% c("fill","both")) {
    if(color_palette %in% c(valid_color_palettes_ggsci)) {
      ggplot_object <- ggplot_object + eval(parse(text = paste0("ggsci::scale_fill_",color_palette,"()")))
    }
    else if (color_palette %in% valid_color_palettes_RColorBrewer) {
      ggplot_object <- ggplot_object + ggplot2::scale_fill_brewer(palette = color_palette)
    }
  }

  if(!is.na(plot_title)){
    ggplot_object <- ggplot_object + ggplot2::ggtitle(plot_title)
  }

  return(ggplot_object)

}


#' Title
#'
#' @param input_data data frame (or tibble) with input values
#' @param groupingFactor factor to group input table by
#' @param valuesColumn column with numeric values to calculate density on
#' @param density_bw bandwidth of density (bw argument to \link{density} function, default=1 )
#' @param kernel kenrel to use
#' @param shift_bands shift bands by specified value
#' @param kernel_from from which value start density estimation
#' @param kernel_to to which value perform density estimation
#' @param scale_by_size should densities be scaled by sample size
#' @param scaling_vector vector containing scaling factors
#'
#' @return list containing plot and density estimates (data.frame)
#' @export
#'
plot_virtual_gel <- function(input_data, groupingFactor, valuesColumn, density_bw = 0.6, kernel="gaussian",shift_bands=0,kernel_from=0,kernel_to=NA,scale_by_size=T,scaling_vector=NA) {

  data_for_gel <- input_data  %>% dplyr::select(!!rlang::sym(groupingFactor),!!rlang::sym(valuesColumn)) %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::mutate(no_sequences=n()) %>% tidyr::nest()



  if (!is.na(scaling_vector)) {
    #data_for_gel$scale <- scaling_vector
    data_counts <- input_data  %>% dplyr::select(!!rlang::sym(groupingFactor),!!rlang::sym(valuesColumn)) %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::summarise(no_sequences=n())
    data_counts <- tibble::enframe(scaling_vector) %>% dplyr::inner_join(data_counts,by=c("name" = groupingFactor))
    data_counts <- data_counts %>% dplyr::ungroup() %>% dplyr::mutate(data_ratio=no_sequences/value) %>% dplyr::mutate(scale = data_ratio/max(data_ratio)) %>% dplyr::select(name,scale)
    #data_counts$scale <- data_counts$data_ratio/max(data_counts$data_ratio)
    data_for_gel <- data_for_gel %>% dplyr::inner_join(data_counts,by=setNames("name",groupingFactor))
  }



  if(is.na(kernel_from)) {
    kernel_from=0
  }
  if(is.na(kernel_to)) {
    kernel_to=max(input_data[[valuesColumn]])
  }

  #### TODO: improve scale_by_size approach
  if (scale_by_size) {
    if(!is.na(scaling_vector)) {
      data_for_gel <- data_for_gel %>% dplyr::group_by(!!rlang::sym(groupingFactor)) %>% dplyr::mutate(dens_x = purrr::map(data,~density(.x[[valuesColumn]],bw=density_bw,kernel=kernel,from=kernel_from,to=kernel_to)$x),dens_y = purrr::map(data,~density(.x[[valuesColumn]],bw=density_bw,kernel=kernel,from=kernel_from,to=kernel_to)$y*scale)) %>% dplyr::select(-data) %>% tidyr::unnest()
    }
    else{
      data_for_gel <- data_for_gel %>% dplyr::mutate(dens_x = purrr::map(data,~density(.x[[valuesColumn]],bw=density_bw,kernel=kernel,from=kernel_from,to=kernel_to)$x),dens_y = purrr::map(data,~density(.x[[valuesColumn]],bw=density_bw,kernel=kernel,from=kernel_from,to=kernel_to)$y*.x$no_sequences[1])) %>% dplyr::select(-data) %>% tidyr::unnest()
    }
  }
  else {
    data_for_gel <- data_for_gel %>% dplyr::mutate(dens_x = purrr::map(data,~density(.x[[valuesColumn]],bw=density_bw,kernel=kernel,from=kernel_from,to=kernel_to)$x),dens_y = purrr::map(data,~density(.x[[valuesColumn]],bw=density_bw,kernel=kernel,from=kernel_from,to=kernel_to)$y)) %>% dplyr::select(-data) %>% tidyr::unnest()
  }
  virtual_gel_plot<-ggplot2::ggplot(data_for_gel,ggplot2::aes(x=!!rlang::sym(groupingFactor),y=dens_x+shift_bands,fill=dens_y)) + ggplot2::geom_tile(width=0.9,show.legend = F) + ggplot2::theme_minimal() + ggplot2::scale_fill_gradient(low="white",high="black") + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank()) + ggplot2::ylab(valuesColumn)

  return(list(data=data_for_gel,plot=virtual_gel_plot))


}
