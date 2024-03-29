% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polya_plots.R
\name{plot_MA}
\alias{plot_MA}
\title{Plots MA plot of differential expression analysis}
\usage{
plot_MA(input_data, transcript_id_column, labels = FALSE, nlabels = 10, ...)
}
\arguments{
\item{input_data}{a table with output from \link{calculate_diff_exp_binom}}

\item{transcript_id_column}{column used for transcript id}

\item{labels}{show point labels using ggrepel}

\item{nlabels}{number of labels to show}

\item{...}{parameters passed to .basic_aesthetics function (scale_x_limit_low = NA, scale_x_limit_high = NA, scale_y_limit_low = NA, scale_y_limit_high = NA, color_palette = "Set1",plot_title=NA)}
}
\value{
\link[ggplot2]{ggplot} object
}
\description{
Crates simple MA plot, with log10(mean expression) on the X-axis and log2(fold_change) on the Y-axis
}
