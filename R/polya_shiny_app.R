# TODO ShinyApp: Allow to choose columns shown in the table
# TODO ShinyApp: Allow launching on raw polya data.frame, calculate statistics using selected factor in the app
# TODO ShinyApp: Export report
# TODO allow statistics on more than 2 levels of factor (anova?)
# TODO Add withProgress

# wrapper for NanoTail Shiny interface
#' Title
#'
#' @param summary_table - data frame with summarized polya lengths and statistics per transcript
#' @param polya_table - data.table object with polyA data for all reads
#'
#' @return
#' @export
#'
nanoTailApp <- function(polya_table,summary_table) {

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
  #if (missing(summary_table)) {
  #  stop("The statistics summary table  (argument summary_table) is missing",
  #       call. = FALSE)
  #}

  if (missing(polya_table)) {
    stop("The polyA lengths table  (argument polya_table) is missing",
         call. = FALSE)
  }

  #assert_that(has_rows(summary_table),msg = "Empty data frame provided as an input (summary_table). Please provide the correct summary_table")
  assert_that(has_rows(polya_table),msg = "Empty data frame provided as an input (polya_table). Please provide the correct polya_table")

  # Calculate processing statistics
  # TODO Move to reactive values
  print("processing_stats")
  processing_statistics = get_nanopolish_processing_info(polya_table)
  print("filter")
  #remove failed reads from polya_table before further analysis
  polya_table_passed <- remove_failed_reads(polya_table)


  # Convert polya_table to data.table (if required)
  # TODO check if is already datatable or not
  polya_table_passed <- data.table::data.table(polya_table_passed)
  data.table::setkey(polya_table_passed,"transcript")

  number_of_samples = length(levels(polya_table_passed$sample_name))
  number_of_transcripts = length(levels(as.factor(polya_table_passed$transcript)))
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

        shinydashboard::menuItem("Basic Info", icon = icon("info-square"), tabName = "basicInfo"),
        shinydashboard::menuItem("Global distribution", icon = icon("chart-line-down"), tabName = "global_distr"),
        shinydashboard::menuItem("Per transcript plots", icon = icon("dna"), tabName = "boxplots"),
        shinydashboard::menuItem("Differential expression", icon = icon("dna"), tabName = "diff_exp"),
        shinydashboard::menuItem("Plot settings", icon = icon("settings"), tabName = "plot_settings"),
        shinydashboard::menuItem("About", tabName = "dashboard", icon = icon("dashboard"))
      ),selectInput("groupingFactor","Group by",choices=polya_table_passed %>% dplyr::select_if(is.factor) %>% colnames),
      uiOutput("condition1UI"),
      uiOutput("condition2UI")
),

    # body definition ---------------------------------------------------------
    shinydashboard::dashboardBody(

      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "dashboard",
                                fillPage(padding = 0, title = NULL, bootstrap = F, theme = NULL,
                                         wellPanel(style = "background-color: #ffffff; overflow-y:scroll; max-height: 750px;",
                                                   p("NanoTail - Interactive exploration of polyA lengths estimations")))
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
        shinydashboard::tabItem(tabName = "plot_settings",
                                box(
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
                                sliderInput("scale_limit_low","Plot scale limit low",0,1024,0,10),
                                sliderInput("scale_limit_high","Plot scale limit high",50,1024,200,10),
                                sliderInput("counts_scatter_low_value","Counts scatter lower counts limit",0,10024,0,50),
                                sliderInput("counts_scatter_high_value","Counts scatter high counts limit",100,10024,1000,50)
                                )),
        shinydashboard::tabItem(tabName = "global_distr",
                                h2("Global distribution of tails"),
                                fluidRow(
                                  box(),
                                  box(plotly::plotlyOutput('polya_global') %>% shinycssloaders::withSpinner(type = 4),collapsible=TRUE)
                                )),
        shinydashboard::tabItem(tabName = "boxplots",
                                fluidRow(
                                  box(actionButton("compute_diff_polya",
                                                     HTML("Compute Polya statistics"),
                                                     icon = icon("spinner")),
                                      DT::dataTableOutput('transcripts_table'),collapsible = TRUE,height = '100%'),
                                  tabBox(title = "per transcript plots",id="tabset1", height="500px",
                                         tabPanel("boxplot",plotly::plotlyOutput('polya_boxplot') %>% shinycssloaders::withSpinner(type = 4)),
                                         tabPanel("distribution plot",plotly::plotlyOutput('polya_distribution') %>% shinycssloaders::withSpinner(type = 4))
                                )
                                )),
        shinydashboard::tabItem(tabName = "diff_exp",
                                fluidRow(
                                  box(actionButton("compute_diff_exp",
                                                   HTML("Compute Differential Expression using Binomial Test"),
                                                   icon = icon("spinner")),
                                  DT::dataTableOutput('diff_exp_table')),
                                  tabBox(title = "plots",id="tabset2",height="500px",
                                    tabPanel("pca_biplot",plotly::plotlyOutput('pca_biplot') %>% shinycssloaders::withSpinner(type = 4)),
                                    tabPanel("volcano_plot",plotly::plotlyOutput('volcano') %>% shinycssloaders::withSpinner(type = 4)),
                                    tabPanel("counts_plot",plotly::plotlyOutput('counts_plot') %>% shinycssloaders::withSpinner(type = 4))
                                  )
        )
        ))))



  shinyApp(ui = nanotail_shiny_ui, server = function(input, output) {


    values <- reactiveValues()
    values$polya_table <- polya_table_passed




    output$transcripts_table = DT::renderDataTable(values$summary_table, server = TRUE, selection='single',options = list(dom = 'ftip'))
    output$diff_exp_table = DT::renderDataTable(values$diffexp, server = TRUE, selection='single')


    output$polya_boxplot = plotly::renderPlotly({
      summary_table = values$summary_table
      #selected_row <- input$table_rows_selected
      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript = data_transcript()
      polya_boxplot <- plot_polya_boxplot(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_y_limit_low = input$scale_limit_low,scale_y_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse,plot_title = selected_transcript)
      plotly::ggplotly(polya_boxplot)

    })

    output$polya_global = plotly::renderPlotly({
      #selected_row <- input$table_rows_selected
      #ggplot(polyA_all,aes(x=polya_length,color=group)) + geom_density(size=1,aes(y=..ndensity..)) + scale_x_continuous(limits=c(0,128)) + theme_bw()

      global_distribution_plot <- plot_polya_distribution(polya_data = polya_table_passed,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse)

      plotly::ggplotly(global_distribution_plot)

    })

    output$polya_distribution = plotly::renderPlotly({
      summary_table = values$summary_table

      #selected_row <- input$table_rows_selected
      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript = data_transcript()
      #ggplot(polyA_all,aes(x=polya_length,color=group)) + geom_density(size=1,aes(y=..ndensity..)) + scale_x_continuous(limits=c(0,128)) + theme_bw()
      transcript_distribution_plot <- plot_polya_distribution(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse)
      plotly::ggplotly(transcript_distribution_plot)

    })

    output$condition1UI <- renderUI({
      group_factor = input$groupingFactor
      selectizeInput("condition1_diff_exp","Condition 1",choices=levels(polya_table_passed[[group_factor]]),selected = levels(polya_table_passed[[group_factor]])[1])
      #selectInput("condition2_diff_exp","Condition 2",choices=levels(polya_table[[group_factor]]),selected = levels(levels(polya_table[[group_factor]])[2]))

    })
    output$condition2UI <- renderUI({
      group_factor = input$groupingFactor
      #selectInput("condition1_diff_exp","Condition 1",choices=levels(polya_table[[group_factor]]),selected = levels(polya_table[[group_factor]])[1])
      selectizeInput("condition2_diff_exp","Condition 2",choices=levels(polya_table_passed[[group_factor]]),selected = levels(polya_table_passed[[group_factor]])[2])

    })

    data_transcript <- shiny::reactive({
      summary_table = values$summary_table

      selected_row <- input$transcripts_table_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      message(selected_transcript)
      data_transcript = subset(polya_table_passed, transcript==selected_transcript)
      data_transcript %<>% dplyr::mutate(replicate = gsub(".*(.)$","\\1", sample_name))

      data_transcript
    })


    output$counts_plot <- plotly::renderPlotly({

      # minimal value of counts for a transcript to show, defined as a 25% quantile of all counts
      #min_counts_quantile =  quantile(values$polya_table_summarized[["counts"]],na.rm=TRUE,probs=c(25, 99)/100)[1]
      # maximal value of counts for a transcript to show, defined as a 95% quantile of all counts, as usually most abundant transcripts locate further from the main body
      #max_counts_quantile =  quantile(values$polya_table_summarized[["counts"]],na.rm=TRUE,probs=c(25, 99)/100)[2]
      counts_scatter_plot <- plot_counts_scatter(polya_data = values$polya_table_summarized,groupingFactor = input$groupingFactor,color_palette = input$col_palette,reverse_palette = input$reverse,condition1 = input$condition1_diff_exp, condition2 = input$condition2_diff_exp,min_counts = input$counts_scatter_low_value, max_counts = input$counts_scatter_high_value)

      counts_scatter_plot


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
        "Reads analyzed ", processing_statistics$number_all_reads, icon = icon("stream"),
        color = "black"
      )
    })

    output$numberOfPassReads <- renderInfoBox({
      infoBox(
        "Reads passing filters ", processing_statistics$number_pass_reads, icon = icon("check"),
        color = "green"
      )
    })

    output$numberOfAdapterReads <- renderInfoBox({
      infoBox(
        "Adapter reads ", processing_statistics$number_adapter_reads, icon = icon("times"),
        color = "red"
      )
    })

    output$numberOfReadLoadFailedReads <- renderInfoBox({
      infoBox(
        "Reads failed to load ", processing_statistics$number_read_load_failed_reads, icon = icon("times"),
        color = "red"
      )
    })

    output$numberOfNoregionReads <- renderInfoBox({
      infoBox(
        "NoRegion reads ", processing_statistics$number_noregion_reads, icon = icon("times"),
        color = "red"
      )
    })


    observeEvent(input$compute_diff_exp,
                 {
                   withProgress(message="Computing Differential epxression using binomial test...",
                                detail = "This step can take a little while",
                                value = 0,{
                                  print(input$condition1_diff_exp)
                                  print(input$condition2_diff_exp)
                                  print(input$groupingFactor)
                                  print("summarizing")
                                  values$polya_table_summarized <- summarize_polya(values$polya_table,summary_factors = c(input$groupingFactor))
                                  print("diffexp")

                                  values$diffexp <- calculate_diff_exp_binom(values$polya_table_summarized,grouping_factor = input$groupingFactor,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,summarized_input=TRUE)
                                  values$condition1 <- input$condition1_diff_exp
                                  values$condition2 <- input$condition2_diff_exp
                                })
                 })


    observeEvent(input$compute_diff_polya,
                 {
                   withProgress(message="Computing polya statistics...",
                                detail = "This step can take a little while",
                                value = 0,{
                                  print(input$condition1_diff_exp)
                                  print(input$condition2_diff_exp)
                                  print(input$groupingFactor)
                                  print("calculating stats")
                                  condition1 = input$condition1_diff_exp
                                  condition2 = input$condition2_diff_exp

                                  polya_stats <- calculate_polya_stats(values$polya_table,grouping_factor = input$groupingFactor,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp)
                                  values$summary_table <- polya_stats$summary %>% dplyr::select(transcript,!!rlang::sym(paste0(condition1,"_counts")),!!rlang::sym(paste0(condition2,"_counts")),!!rlang::sym(paste0(condition1,"_polya_median")),!!rlang::sym(paste0(condition2,"_polya_median")),p.value,p.corr)
                                  values$condition1 <- input$condition1_diff_exp
                                  values$condition2 <- input$condition2_diff_exp
                                })
                 })


  })
}
