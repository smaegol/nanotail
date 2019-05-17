#' Title
#'
#' @param pca_object
#' @param samples_names
#'
#' @return
#' @export
#'
#' @examples
plot_polyA_PCA <- function(pca_object,samples_names) {

  ggbiplot::ggbiplot(pca_object,var.axes = F,labels = samples_names)

}


#' Title
#'
#' @param polya_data
#' @param groupingFactor
#' @param scale_x_limit_low
#' @param scale_x_limit_high
#' @param color_palette
#' @param reverse_palette
#'
#' @return
#' @export
#'
#' @examples
plot_polya_distribution <- function(polya_data, groupingFactor, scale_x_limit_low=NA, scale_x_limit_high=NA, color_palette = "Set1", reverse_palette = 0) {

  distribution_plot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x="polya_length",color=groupingFactor)) + ggplot2::geom_density(size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::theme_bw()

  if (!is.na(scale_x_limit_low)) {
    if (is.na(scale_x_limit_high)) {
      stop("Please provide upper limit for x scale")
    }
    distribution_plot <- distribution_plot + ggplot2::scale_x_continuous(limits=c(scale_x_limit_low,scale_x_limit_high))
  }
  if (reverse_palette) {
    distribution_plot <- distribution_plot + ggplot2::scale_colour_brewer(palette = color_palette,direction=-1)
  }
  else {
    distribution_plot <- distribution_plot + ggplot2::scale_colour_brewer(palette = color_palette)
  }

  return(distribution_plot)

}


#' Title
#'
#' @param polya_data
#' @param groupingFactor
#' @param scale_y_limit_low
#' @param scale_y_limit_high
#' @param color_palette
#' @param reverse_palette
#' @param plot_title
#'
#' @return
#' @export
#'
#' @examples
plot_polya_boxplot <- function(polya_data, groupingFactor, scale_y_limit_low=NA, scale_y_limit_high=NA, color_palette = "Set1", reverse_palette = 0, plot_title = "") {

  transcripts_boxplot <- ggplot2::ggplot(polya_data,ggplot2::aes_string(x=groupingFactor,y="polya_length")) + ggplot2::geom_boxplot()  + ggplot2::ggtitle(plot_title)
  if (!is.na(scale_y_limit_low)) {
    if (is.na(scale_y_limit_high)) {
      stop("Please provide both limits for y scale")
    }
    transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_y_continuous(limits=c(scale_y_limit_low,scale_y_limit_high))
  }
  if (reverse_palette) {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_colour_brewer(palette = color_palette,direction=-1)
  }
  else {
    transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_colour_brewer(palette = color_palette)
  }
}


#' Title
#'
#' @param polya_data_summarized
#' @param groupingFactor
#' @param color_palette
#' @param reverse_palette
#' @param condition1
#' @param condition2
#' @param min_counts
#' @param max_counts
#'
#' @return
#' @export
#'
#' @examples
plot_counts_scatter <- function(polya_data_summarized, groupingFactor, color_palette = "Set1", reverse_palette = 0,condition1 = NA, condition2 = NA,min_counts=0,max_counts=NA) {


  polya_data_summarized_counts_xy<-polya_data_summarized %>% dplyr::group_by(transcript,!!rlang::sym(groupingFactor)) %>% dplyr::summarize(counts_sum=sum(counts)) %>% tidyr::spread(group,counts_sum)

  polya_data_summarized_counts_xy[is.na(polya_data_summarized_counts_xy)] <- 0

  polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)>=min_counts,!!rlang::sym(condition2)>=min_counts)

  if (!is.na(max_counts)) {
    polya_data_summarized_counts_xy <- polya_data_summarized_counts_xy %>% dplyr::filter(!!rlang::sym(condition1)<=max_counts,!!rlang::sym(condition2)<=max_counts)
  }

  counts_scatter_plot<-ggplot2::ggplot(polya_data_summarized_counts_xy,ggplot2::aes(x=!!rlang::sym(condition1),y=!!rlang::sym(condition2))) + ggplot2::geom_point(ggplot2::aes(text=transcript),alpha=0.7)

  if (reverse_palette) {
    counts_scatter_plot <- counts_scatter_plot + ggplot2::scale_colour_brewer(palette = color_palette,direction=-1)
  }
  else {
    counts_scatter_plot <- counts_scatter_plot + ggplot2::scale_colour_brewer(palette = color_palette)
  }

}
