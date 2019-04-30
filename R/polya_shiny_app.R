# Define the NanoTail shiny interface UI
nanotail_shiny_ui_old <- fluidPage(


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

  # Convert polya_table to data.table (if required)
  # TODO check if is already datatable or not
  polya_table <- data.table::data.table(polya_table)
  data.table::setkey(polya_table,"transcript")

  # Calculate processing statistics
  processing_statistics = calculate_processing_statistics(polya_table)

  #remove failed reads from polya_table
  polya_table <- remove_failed_reads(polya_table)

  number_of_samples = length(levels(polya_table$sample_name))
  number_of_transcripts = length(levels(as.factor(polya_table$transcript)))
  ## colorpicker based on Zbyszek Pietras Bar Plot app (https://github.com/zetp/Bar_plot_NC/blob/master/app.R)
  ### here is code to make color scale picker (palette picker) UI input
  ### it is taken from shinyWidgets website: https://dreamrs.github.io/shinyWidgets/articles/palette_picker.html

  # List of palettes
  colors_pal <- lapply(
    X = split(
      x = brewer.pal.info,
      f = factor(brewer.pal.info$category, labels = c("Diverging", "Qualitative", "Sequential"))
    ),
    FUN = rownames
  )

  # Get all colors given a palette name(s)
  get_brewer_name <- function(name) {
    pals <- brewer.pal.info[rownames(brewer.pal.info) %in% name, ]
    res <- lapply(
      X = seq_len(nrow(pals)),
      FUN = function(i) {
        brewer.pal(n = pals$maxcolors[i], name = rownames(pals)[i])
      }
    )
    unlist(res)
  }

  background_pals <- sapply(unlist(colors_pal, use.names = FALSE), get_brewer_name)

  # Calc linear gradient for CSS
  linear_gradient <- function(cols) {
    x <- round(seq(from = 0, to = 100, length.out = length(cols)+1))
    ind <- c(1, rep(seq_along(x)[-c(1, length(x))], each = 2), length(x))
    m <- matrix(data = paste0(x[ind], "%"), ncol = 2, byrow = TRUE)
    res <- lapply(
      X = seq_len(nrow(m)),
      FUN = function(i) {
        paste(paste(cols[i], m[i, 1]), paste(cols[i], m[i, 2]), sep = ", ")
      }
    )
    res <- unlist(res)
    res <- paste(res, collapse = ", ")
    paste0("linear-gradient(to right, ", res, ");")
  }

  background_pals <- unlist(lapply(X = background_pals, FUN = linear_gradient))
  head(background_pals, 3)

  colortext_pals <- rep(c("white", "black", "black"), times = sapply(colors_pal, length))


  nanotail_shiny_ui <- shinydashboard::dashboardPage(
    skin="blue",
    # header definition -----------------------------------------------------------
    shinydashboard::dashboardHeader(
      title = paste0("NanoTail - Interactive exploration of polyA lengths estimations ",
                     "done using Nanopore sequencing - version ",
                     packageVersion("nanotail")),
      titleWidth = 900),
    # sidebar definition -----------------------------------------------------------
    shinydashboard::dashboardSidebar(
      width = 180,
      shinydashboard::sidebarMenu(id="tabs",
        shinydashboard::menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
        shinydashboard::menuItem("Basic Info", icon = icon("info-square"), tabName = "basicInfo"),
        shinydashboard::menuItem("Global distribution", icon = icon("chart-line-down"), tabName = "global_distr"),
        shinydashboard::menuItem("Per transcript plots", icon = icon("dna"), tabName = "boxplots")
      ),
      shinyWidgets::pickerInput(inputId = "col_palette", label = "Colour scale",
                  choices = colors_pal, selected = "Set1", width = "100%",
                  choicesOpt = list(
                    content = sprintf(
                      "<div style='width:100%%;padding:5px;border-radius:4px;background:%s;color:%s'>%s</div>",
                      unname(background_pals), colortext_pals, names(background_pals)
                    )
                  )
      ),
      uiOutput("Colorblind"),
      checkboxInput("reverse", "reverse color scale", value = FALSE),
      selectInput("groupingFactor","Group by",choices=polya_table %>% dplyr::select_if(is.factor) %>% colnames)
    ),

    # body definition ---------------------------------------------------------
    shinydashboard::dashboardBody(

      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "dashboard",
                                h2("Nanotail")
        ),
        shinydashboard::tabItem(tabName = "basicInfo",
                                infoBoxOutput("numberOfSamples"),
                                infoBoxOutput("numberOfTranscripts"),
                                infoBoxOutput("numberOfAllReads"),
                                infoBoxOutput("numberOfPassReads"),
                                infoBoxOutput("numberOfAdapterReads"),
                                infoBoxOutput("numberOfReadLoadFailedReads"),
                                infoBoxOutput("numberOfNoregionReads")

        ),
        shinydashboard::tabItem(tabName = "global_distr",
                                h2("Global distribution of tails"),
                                fluidRow(
                                  column(4),
                                  column(6,plotly::plotlyOutput('polya_global') %>% shinycssloaders::withSpinner(type = 4))
                                )),
        shinydashboard::tabItem(tabName = "boxplots",
                                fluidRow(
                                  box(DT::dataTableOutput('transcripts_table',width = "40%")),
                                  box(plotly::plotlyOutput('polya_boxplot') %>% shinycssloaders::withSpinner(type = 4))
                                ),
                                fluidRow(
                                  box(),
                                  box(plotly::plotlyOutput('polya_distribution') %>% shinycssloaders::withSpinner(type = 4))
                                )
        ))))



  shinyApp(ui = nanotail_shiny_ui, server = function(input, output) {


    output$transcripts_table = DT::renderDataTable(summary_table, server = TRUE, selection='single')


    output$polya_boxplot = plotly::renderPlotly({
      #selected_row <- input$table_rows_selected
      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript = data_transcript()
      transcripts_boxplot <- ggplot2::ggplot(data_transcript,ggplot2::aes_string(x=input$groupingFactor,y="polya_length")) + ggplot2::geom_boxplot()  + ggplot2::ggtitle(selected_transcript)
      if (input$reverse) {
        transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_colour_brewer(palette = input$col_palette,direction=-1)
      }
      else {
        transcripts_boxplot <- transcripts_boxplot + ggplot2::scale_colour_brewer(palette = input$col_palette)
      }
      #if (input$split_replicates == TRUE) {
      #  transcripts_boxplot <- transcripts_boxplot  + ggplot2::facet_wrap(. ~ replicate)
      #}
      plotly::ggplotly(transcripts_boxplot)

    })

    output$polya_global = plotly::renderPlotly({
      #selected_row <- input$table_rows_selected
      #ggplot(polyA_all,aes(x=polya_length,color=group)) + geom_density(size=1,aes(y=..ndensity..)) + scale_x_continuous(limits=c(0,128)) + theme_bw()
      global_distribution_plot <- ggplot2::ggplot(polya_table,ggplot2::aes_string(x="polya_length",color=input$groupingFactor)) + ggplot2::geom_density(size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::scale_x_continuous(limits=c(0,128)) + ggplot2::theme_bw()
      if (input$reverse) {
        global_distribution_plot <- global_distribution_plot + ggplot2::scale_colour_brewer(palette = input$col_palette,direction=-1)
      }
      else {
        global_distribution_plot <- global_distribution_plot + ggplot2::scale_colour_brewer(palette = input$col_palette)
      }
      #if (input$split_replicates == TRUE) {
      #  transcripts_boxplot <- transcripts_boxplot  + ggplot2::facet_wrap(. ~ replicate)
      #}
      plotly::ggplotly(global_distribution_plot)

    })

    output$polya_distribution = plotly::renderPlotly({
      #selected_row <- input$table_rows_selected
      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript = data_transcript()
      #ggplot(polyA_all,aes(x=polya_length,color=group)) + geom_density(size=1,aes(y=..ndensity..)) + scale_x_continuous(limits=c(0,128)) + theme_bw()
      transcript_distribution_plot <- ggplot2::ggplot(data_transcript,ggplot2::aes_string(x="polya_length",color=input$groupingFactor)) + ggplot2::geom_density(size=1,ggplot2::aes(y=..ndensity..)) + ggplot2::scale_x_continuous(limits=c(0,128)) + ggplot2::theme_bw()
      if (input$reverse) {
        transcript_distribution_plot <- transcript_distribution_plot + ggplot2::scale_colour_brewer(palette = input$col_palette,direction=-1)
      }
      else {
        transcript_distribution_plot <- transcript_distribution_plot + ggplot2::scale_colour_brewer(palette = input$col_palette)
      }
      #if (input$split_replicates == TRUE) {
      #  transcripts_boxplot <- transcripts_boxplot  + ggplot2::facet_wrap(. ~ replicate)
      #}
      plotly::ggplotly(transcript_distribution_plot)

    })

    data_transcript <- shiny::reactive({
      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      message(selected_transcript)
      data_transcript = subset(polya_table, transcript==selected_transcript)
      data_transcript %<>% dplyr::mutate(replicate = gsub(".*(.)$","\\1", sample_name))

      data_transcript
    })


    output$Colorblind <- renderUI({
      if(brewer.pal.info %>% subset(rownames(.) == input$col_palette) %>% .$colorblind) {
        txt_ <- "Yes"
      } else {
        txt_ <- "No"
      }
      HTML(paste("<p>Is colour scale colour blind friendly?<b>", txt_, "</b></p>"))
    })


    output$numberOfSamples <- renderInfoBox({
      infoBox(
        "Samples to Compare ", number_of_samples, icon = icon("vials"),
        color = "yellow"
      )
    })

    output$numberOfTranscripts <- renderInfoBox({
      infoBox(
        "Transcripts analyzed ", number_of_transcripts, icon = icon("dna"),
        color = "blue"
      )
    })

    output$numberOfAllReads <- renderInfoBox({
      infoBox(
        "Reads analyzed ", processing_statistics$number_all_reads, icon = icon("dna"),
        color = "black"
      )
    })

    output$numberOfPassReads <- renderInfoBox({
      infoBox(
        "Reads passing filters ", processing_statistics$number_pass_reads, icon = icon("dna"),
        color = "green"
      )
    })

    output$numberOfAdapterReads <- renderInfoBox({
      infoBox(
        "Adapter reads ", processing_statistics$number_adapter_reads, icon = icon("dna"),
        color = "red"
      )
    })

    output$numberOfReadLoadFailedReads <- renderInfoBox({
      infoBox(
        "Reads failed to load ", processing_statistics$number_read_load_failed_reads, icon = icon("dna"),
        color = "red"
      )
    })

    output$numberOfNoregionReads <- renderInfoBox({
      infoBox(
        "NoRegion reads ", processing_statistics$number_noregion_reads, icon = icon("dna"),
        color = "red"
      )
    })


  })
}
