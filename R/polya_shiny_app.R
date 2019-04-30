# Define the NanoTail shiny interface UI
nanotail_shiny_ui <- fluidPage(


  # Application title
  titlePanel("NanoTail"),

  fluidRow(
    column(6,
           DT::dataTableOutput('transcripts_table')),
    column(4,plotly::plotlyOutput('polya_boxplot') %>% shinycssloaders::withSpinner(type = 4))
  )
)




# wrapper for NanoTail Shiny interface
#' Title
#'
#' @param summary_table - data frame with summarized polya lengths and statistics per transcript
#' @param polya_table - data.table object with polyA data for all reads
#'
#' @return
#' @export
#'
#' @examples
nanoTailApp <- function(summary_table,polya_table) {

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("NanoTail requires 'shiny'. Please install it using
         install.packages('shiny')")
  }
  if ( !requireNamespace('shinycssloaders',quietly = TRUE) ) {
    stop("pcaExplorer requires 'shinycssloaders'. Please install it using
         install.packages('shinycssloaders')")
  }
  if ( !requireNamespace('plotly',quietly = TRUE) ) {
    stop("pcaExplorer requires 'plotly'. Please install it using
         install.packages('plotly')")
  }

  #check if parameters are provided
  if (missing(summary_table)) {
    stop("The statistics summary table  (argument summary_table) is missing",
         call. = FALSE)
  }

  if (missing(polya_table)) {
    stop("The polyA lengths table  (argument polya_table) is missing",
         call. = FALSE)
  }

  assert_that(has_rows(summary_table),msg = "Empty data frame provided as an input (summary_table). Please provide the correct summary_table")
  assert_that(has_rows(summary_table),msg = "Empty data frame provided as an input (polya_table). Please provide the correct polya_table")

  shinyApp(ui = nanotail_shiny_ui, server = function(input, output) {


    output$transcripts_table = DT::renderDataTable(summary_table, server = TRUE, selection='single')


    output$polya_boxplot = plotly::renderPlotly({
      #selected_row <- input$table_rows_selected
      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript = data_transcript()
      isoforms_plot <- ggplot2::ggplot(data_transcript,ggplot2::aes(x=group,y=polya_length)) + ggplot2::geom_boxplot()  + ggplot2::ggtitle(selected_transcript)
      #if (input$split_replicates == TRUE) {
      #  isoforms_plot <- isoforms_plot  + ggplot2::facet_wrap(. ~ replicate)
      #}
      plotly::ggplotly(isoforms_plot)

    })

    data_transcript <- shiny::reactive({
      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      message(selected_transcript)
      data_transcript = subset(polya_table, transcript==selected_transcript)
      data_transcript %<>% dplyr::mutate(replicate = gsub(".*(.)$","\\1", sample_name))

      data_transcript
    })

  })
}
