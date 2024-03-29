% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polya_shiny_app.R
\name{nanoTailApp}
\alias{nanoTailApp}
\title{wrapper for NanoTail Shiny interface}
\usage{
nanoTailApp(
  polya_table,
  precomputed_polya_statistics = NA,
  precomputed_annotations = NA
)
}
\arguments{
\item{polya_table}{polyA predictions table. Can be obtained with \link{read_polya_multiple}}

\item{precomputed_polya_statistics}{precomputed differential adenylation table (obtained with \link{calculate_polya_stats})}

\item{precomputed_annotations}{precomputed annotations}
}
\description{
wrapper for NanoTail Shiny interface
}
