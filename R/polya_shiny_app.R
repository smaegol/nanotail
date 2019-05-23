# TODO ShinyApp: Allow to choose columns shown in the table
# TODO ShinyApp: Allow launching on raw polya data.frame, calculate statistics using selected factor in the app
# TODO ShinyApp: Export report
# TODO allow statistics on more than 2 levels of factor (anova?)
# TODO Add withProgress

#' wrapper for NanoTail Shiny interface
#'
#' @param polya_table polyA predictions table. Can be obtained with \link{read_polya_multiple}
#'
#' @export
#'

nanoTailApp <- function(polya_table) {

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("NanoTail requires 'shiny'. Please install it using
         install.packages('shiny')")
  }
  if ( !requireNamespace('shinycssloaders',quietly = TRUE) ) {
    stop("NanoTail requires 'shinycssloaders'. Please install it using
         install.packages('shinycssloaders')")
  }
  if ( !requireNamespace('plotly',quietly = TRUE) ) {
    stop("NanoTail requires 'plotly'. Please install it using
         install.packages('plotly')")
  }

  #check if parameters are provided

  if (missing(polya_table)) {
    stop("The polyA lengths table  (argument polya_table) is missing",
         call. = FALSE)
  }


  assertthat::assert_that(assertive::has_rows(polya_table),msg = "Empty data frame provided as an input (polya_table). Please provide the correct polya_table")
  assertthat::assert_that("polya_length" %in% colnames(polya_table),msg = "Input polya_table should contain polya length predictions. Did you provide valid input?")
  assertthat::assert_that("transcript" %in% colnames(polya_table),msg = "Input polya_table should contain 'transcript' column. Did you provide valid input?")
  assertthat::assert_that("sample_name" %in% colnames(polya_table),msg = "Input polya_table should contain at least `sample_name` column for grouping. If data were read using read_polya_single() function please use read_polya_single(...,sample_name= SAMPLE_NAME ) or read_polya_multiple() to provide metadata")
  assertthat::assert_that(is.factor(polya_table$sample_name),msg = "Sample_name column should be a factor")



  # Calculate processing statistics
  # TODO Move to reactive values
  #remove failed reads from polya_table before further analysis
  polya_table_passed <- remove_failed_reads(polya_table)

  # make sure transcript names are not factor
  polya_table_passed$transcript <- as.character(polya_table_passed$transcript)

  # Convert polya_table to data.table (if required)
  # TODO check if is already datatable or not
  polya_table_passed <- data.table::data.table(polya_table_passed)
  data.table::setkey(polya_table_passed,"transcript")

  number_of_samples = length(levels(polya_table_passed$sample_name))
  number_of_transcripts = length(levels(as.factor(polya_table_passed$transcript)))


  ### !<- colorpicker start ->! ###
  ## colorpicker based on Zbyszek Pietras Bar Plot app (https://github.com/zetp/Bar_plot_NC/blob/master/app.R)
  ### here is code to make color scale picker (palette picker) UI input
  ### it is taken from shinyWidgets website: https://dreamrs.github.io/shinyWidgets/articles/palette_picker.html

  # List of palettes
  colors_pal <- lapply(
    X = split(
      x = RColorBrewer::brewer.pal.info,
      f = factor(RColorBrewer::brewer.pal.info$category, labels = c("Diverging", "Qualitative", "Sequential"))
    ),
    FUN = rownames
  )

  # Get all colors given a palette name(s)
  get_brewer_name <- function(name) {
    pals <- RColorBrewer::brewer.pal.info[rownames(RColorBrewer::brewer.pal.info) %in% name, ]
    res <- lapply(
      X = seq_len(nrow(pals)),
      FUN = function(i) {
        RColorBrewer::brewer.pal(n = pals$maxcolors[i], name = rownames(pals)[i])
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

  ### !<- colorpicker end ->! ###

  processing_info_global <- get_nanopolish_processing_info(polya_table,grouping_factor = NA)
  processing_info_global_spread <- processing_info_global %>% tidyr::spread(qc_tag,n)
  processing_info_per_sample <- get_nanopolish_processing_info(polya_table,grouping_factor = "sample_name")
  processing_info_per_sample_spread <- processing_info_per_sample %>% tidyr::spread(qc_tag,n)

  polya_summary_table <- summarize_polya(polya_data = polya_table_passed,summary_factors = "sample_name")
  initial_summary_table <- polya_summary_table %>% dplyr::select(transcript,sample_name,polya_median) %>% tidyr::spread(sample_name,polya_median)
  initial_summary_table_counts <- polya_summary_table %>% dplyr::select(transcript,sample_name,counts) %>% tidyr::spread(sample_name,counts)


  nanotail_shiny_ui <- shinydashboard::dashboardPage(
    skin="blue",
    # dashboard header definition -----------------------------------------------------------
    shinydashboard::dashboardHeader(
      title = paste0("NanoTail - Interactive exploration of polyA lengths estimations ",
                     "done using Nanopore sequencing - version ",
                     packageVersion("nanotail")),
      titleWidth = 900),

    # dashboard sidebar definition -----------------------------------------------------------
    shinydashboard::dashboardSidebar(
      width = 180,
      shinydashboard::sidebarMenu(id="nanotail_menu",
        shinydashboard::menuItem("QC info", icon = icon("info-square"), tabName = "basicInfo"),
        shinydashboard::menuItem("Global polyA distribution", icon = icon("chart-line-down"), tabName = "global_distr"),
        shinydashboard::menuItem("Differential adenylation", icon = icon("dna"), tabName = "diff_polya"),
        shinydashboard::menuItem("Differential expression", icon = icon("dna"), tabName = "diff_exp"),
        shinydashboard::menuItem("Plot settings", icon = icon("settings"), tabName = "plot_settings"),
        shinydashboard::menuItem("About", tabName = "dashboard", icon = icon("dashboard"))
      ),
      selectInput("groupingFactor","Group by",choices=polya_table_passed %>% dplyr::select_if(is.factor) %>% colnames),
      uiOutput("condition1UI"),
      uiOutput("condition2UI"),
      checkboxInput("plot_only_selected_conditions","Plot only selected conditions",value=FALSE)
    ),

    # dashboard body definition ---------------------------------------------------------
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "dashboard",
          fillPage(padding = 0, title = NULL, bootstrap = F, theme = NULL,
             wellPanel(style = "background-color: #ffffff; overflow-y:scroll; max-height: 750px;",
                       includeMarkdown(system.file("extdata", "about.md",package = "nanotail")))
          )),
        shinydashboard::tabItem(tabName = "basicInfo",
          fluidRow(
            infoBoxOutput("numberOfSamples"),
            infoBoxOutput("numberOfTranscripts"),
            infoBoxOutput("numberOfAllReads")),
          fluidRow(
            column(1),
            column(7,fluidRow(h4("Per sample plot:")),
               fluidRow(
                  plotly::plotlyOutput("nanopolish_qc_summary_plot"),
                  checkboxInput("show_frequency_plot_nanopolioshQC",label = "Show QC as frequencies?",value = TRUE))),
            column(1),
            column(3,fluidRow(h4("Information for all samples:")),
               fluidRow(infoBoxOutput("numberOfPassReads",width = NULL)),
               fluidRow(infoBoxOutput("numberOfAdapterReads",width = NULL)),
               fluidRow(infoBoxOutput("numberOfReadLoadFailedReads",width = NULL)),
               fluidRow(infoBoxOutput("numberOfNoregionReads",width = NULL))))),
        shinydashboard::tabItem(tabName = "plot_settings",
          box(
            shinyWidgets::pickerInput(inputId = "col_palette", label = "Colour scale",
                choices = colors_pal, selected = "Set1", width = "100%",
                choicesOpt = list(
                  content = sprintf(
                    "<div style='width:100%%;padding:5px;border-radius:4px;background:%s;color:%s'>%s</div>",
                    unname(background_pals), colortext_pals, names(background_pals)))),
          uiOutput("Colorblind"),
          checkboxInput("reverse", "reverse color scale", value = FALSE),
          sliderInput("scale_limit_low","Plot scale limit low",0,1024,0,10),
          sliderInput("scale_limit_high","Plot scale limit high",50,1024,200,10),
          sliderInput("counts_scatter_low_value","Counts scatter lower counts limit",0,10024,0,50),
          sliderInput("counts_scatter_high_value","Counts scatter high counts limit",100,10024,1000,50))),
        shinydashboard::tabItem(tabName = "global_distr",
          h2("Global comparison"),
          fluidRow(
            box(plotly::plotlyOutput('polya_global') %>% shinycssloaders::withSpinner(type = 4),collapsible=TRUE))),
        shinydashboard::tabItem(tabName = "diff_polya",
          fluidRow(
            box(DT::dataTableOutput('diff_polya'),collapsible = TRUE,height = '90%'),
            tabBox(title = "per transcript plots",id="tabset1", height="500px",
                   tabPanel("boxplot",plotly::plotlyOutput('polya_boxplot') %>% shinycssloaders::withSpinner(type = 4)),
                   tabPanel("distribution plot",plotly::plotlyOutput('polya_distribution') %>% shinycssloaders::withSpinner(type = 4)),
                   tabPanel("volcano plot polya",plotly::plotlyOutput('polya_volcano') %>% shinycssloaders::withSpinner(type = 4)))),
          fluidRow(
            box(actionButton("compute_diff_polya",
                             HTML("Compute Polya statistics"),
                             icon = icon("spinner")),
                checkboxInput("use_dwell_time_for_statistics",label = "Use dwell time instead of calculated polya_length for statistics",value = FALSE),
                selectInput("polya_stat_test","statistical test to use",choices = c("Wilcoxon","KS"),selected = "Wilcoxon"),
                sliderInput("min_reads_for_polya_stats","Minimal count of reads per transcript",0,100,10,1),collapsible = TRUE,width=12))),
        shinydashboard::tabItem(tabName = "diff_exp",
          fluidRow(
            box(DT::dataTableOutput('diff_exp_table')),
                tabBox(title = "plots",id="tabset2",height="500px",
                  #tabPanel("pca_biplot",plotOutput('pca_biplot') %>% shinycssloaders::withSpinner(type = 4)),
                  tabPanel("scatter plot of counts",actionButton("show_scatter_plot",
                                                      HTML("Show scatter plot"),
                                                      icon = icon("spinner")),plotly::plotlyOutput('counts_plot') %>% shinycssloaders::withSpinner(type = 4)),
                  tabPanel("volcano_plot",plotly::plotlyOutput('volcano') %>% shinycssloaders::withSpinner(type = 4)),
                  tabPanel("MA_plot",plotly::plotlyOutput('MAplot') %>% shinycssloaders::withSpinner(type = 4)))),
          fluidRow(actionButton("compute_diff_exp",
                                 HTML("Compute Differential Expression using Binomial Test"),
                                 icon = icon("spinner")))))))


  # call to shiny app and server-siude functions
  shinyApp(ui = nanotail_shiny_ui, server = function(input, output) {


    ## Reactive values definition  ---------------------------------------------------------
    values <- reactiveValues()
    values$polya_table <- polya_table_passed
    values$polya_statistics_summary_table = initial_summary_table
    values$diffexp_summary_table = initial_summary_table_counts
    values$polya_table_summarized <- polya_summary_table
    values$processing_info_global <- processing_info_global
    values$processing_info_global_spread <- processing_info_global_spread
    values$processing_info_per_sample <- processing_info_per_sample
    values$processing_info_per_sample_spread <- processing_info_per_sample_spread
    values$diff_exp_grouping_factor = "sample_name"

    # get polyA data for currently selected transcript
    data_transcript <- shiny::reactive({
      summary_table = values$polya_statistics_summary_table

      selected_row <- input$diff_polya_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript = subset(polya_table_passed, transcript==selected_transcript)
      data_transcript <- data_transcript %>% dplyr::mutate(replicate = gsub(".*(.)$","\\1", sample_name))

      data_transcript
    })


    # UI elements rendering ---------------------------------------------------------

    ## Generate select  inputs for conditions used for comparison, based on selected groupingFactor:  ---------------------------------------------------------
    output$condition1UI <- renderUI({
      group_factor = input$groupingFactor
      selectizeInput("condition1_diff_exp","Condition 1",choices=levels(polya_table_passed[[group_factor]]),selected = levels(polya_table_passed[[group_factor]])[1])

    })
    output$condition2UI <- renderUI({
      group_factor = input$groupingFactor
      selectizeInput("condition2_diff_exp","Condition 2",choices=levels(polya_table_passed[[group_factor]]),selected = levels(polya_table_passed[[group_factor]])[2])

    })


    output$Colorblind <- renderUI({
      if(RColorBrewer::brewer.pal.info %>% subset(rownames(.) == input$col_palette) %>% .$colorblind) {
        txt_ <- "Yes"
      } else {
        txt_ <- "No"
      }
      HTML(paste("<p>Is colour scale colour blind friendly?<b>", txt_, "</b></p>"))
    })



    # Output elements rendering section  ---------------------------------------------------------

    ## Use DT::datatable for tables
    output$diff_polya = DT::renderDataTable(values$polya_statistics_summary_table %>% dplyr::rename_all(funs(stringr::str_replace_all(.,"_"," "))), server = TRUE, selection=list(mode = 'single',selected = 1,target = "row"),options = list(dom = 'ftip'))
    output$diff_exp_table = DT::renderDataTable(values$diffexp_summary_table %>% dplyr::rename_all(funs(stringr::str_replace_all(.,"_"," "))), server = TRUE, selection=list(mode = 'single',selected = 1,target = "row"))


    ## Show global polyA distribution as density plot
    output$polya_global = plotly::renderPlotly({

      if (input$plot_only_selected_conditions) {
        global_distribution_plot <- plot_polya_distribution(polya_data = polya_table_passed,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp)
      }
      else {
        global_distribution_plot <- plot_polya_distribution(polya_data = polya_table_passed,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse)
      }
      plotly::ggplotly(global_distribution_plot)

    })

    ## Differential adenylation analysis section  ---------------------------------------------------------


    # Show boxplot of estimated polya lengths for selected transcript
    output$polya_boxplot = plotly::renderPlotly({
      summary_table = values$polya_statistics_summary_table


      if (length(input$diff_polya_rows_selected)>0) {
        selected_row <- input$diff_polya_rows_selected
        selected_transcript = summary_table[selected_row,]$transcript
        data_transcript = data_transcript()
        if (input$plot_only_selected_conditions) {
          polya_boxplot <- plot_polya_boxplot(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_y_limit_low = input$scale_limit_low,scale_y_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse,plot_title = selected_transcript,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp)
        }
        else {
          polya_boxplot <- plot_polya_boxplot(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_y_limit_low = input$scale_limit_low,scale_y_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse,plot_title = selected_transcript)
        }

        plotly::ggplotly(polya_boxplot)
      }
      else {
        suppressWarnings(plotly::plotly_empty())
      }
    })

    # Show denisty plot of estimated polya lengths for selected transcript
    output$polya_distribution = plotly::renderPlotly({
      summary_table = values$polya_statistics_summary_table


      if (length(input$diff_polya_rows_selected)>0) {
        selected_row <- input$diff_polya_rows_selected
        selected_transcript = summary_table[selected_row,]$transcript
        data_transcript = data_transcript()
        if (input$plot_only_selected_conditions) {
          transcript_distribution_plot <- plot_polya_distribution(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse, plot_title = selected_transcript,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp)
        }
        else {
          transcript_distribution_plot <- plot_polya_distribution(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,reverse_palette = input$reverse, plot_title = selected_transcript)
        }
        plotly::ggplotly(transcript_distribution_plot)
      }
      else {
        suppressWarnings(plotly::plotly_empty())
      }
    })

    # Show volcano plot of differntial adenylation analysis
    output$polya_volcano <- plotly::renderPlotly({

      # show only if differential adenylation analysis was already done
      if("fold_change" %in% colnames(values$polya_table_for_volcano)) {
        volcanoPlot <- plot_volcano(input_data =values$polya_table_for_volcano,color_palette = input$col_palette,reverse_palette = input$reverse)
        plotly::ggplotly(volcanoPlot)
      }
      else {
        showNotification("Launch differential polyadenylation analysis first",closeButton=TRUE,type="warning")
        suppressWarnings(plotly::plotly_empty())
      }

    })

    # Calculate differential adenylation if action button was clicked
    observeEvent(input$compute_diff_polya,
    {
      if(number_of_samples>1) {
        withProgress(message="Computing polya statistics...",
          detail = "This step can take a little while",
          value = 0.05,min=0,max=1,
          {
            condition1 = input$condition1_diff_exp
            condition2 = input$condition2_diff_exp
            incProgress(1/20)
            polya_stats <- calculate_polya_stats(values$polya_table,grouping_factor = input$groupingFactor,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,stat_test=input$polya_stat_test,use_dwell_time=input$use_dwell_time_for_statistics,min_reads = input$min_reads_for_polya_stats)
            incProgress(15/20)
            values$polya_statistics_summary_table <- polya_stats$summary %>% dplyr::select(transcript,!!rlang::sym(paste0(condition1,"_counts")),!!rlang::sym(paste0(condition2,"_counts")),!!rlang::sym(paste0(condition1,"_polya_median")),!!rlang::sym(paste0(condition2,"_polya_median")),p.value,padj)
            values$polya_table_for_volcano <- polya_stats$summary %>% dplyr::select(transcript,fold_change,padj,significance)
            values$condition1 <- input$condition1_diff_exp
            values$condition2 <- input$condition2_diff_exp
            incProgress(3/20)
        })
      }
      else {
        showNotification("Differential analysis is not possible when only single sample is present",closeButton=TRUE,type="error")
      }
    })



    ## Differential expression analysis section  ---------------------------------------------------------

    # Show scatter plot of raw counts
    output$counts_plot <- plotly::renderPlotly({

      counts_scatter_plot <- plot_counts_scatter(polya_data = values$polya_table_summarized_scatter,groupingFactor = values$groupingFactor_scatter,color_palette = input$col_palette,reverse_palette = input$reverse,condition1 = values$condition1_scatter, condition2 = values$condition2_scatter, min_counts = input$counts_scatter_low_value, max_counts = input$counts_scatter_high_value)
      plotly::ggplotly(counts_scatter_plot)


    })

    # Show volcano plot of differential expression analysis
    output$volcano <- plotly::renderPlotly({

      # show only if differential expression analysis was already done
      if("fold_change" %in% colnames(values$diffexp_summary_table)) {
        volcanoPlot <- plot_volcano(input_data =values$diffexp_summary_table,color_palette = input$col_palette,reverse_palette = input$reverse)
        plotly::ggplotly(volcanoPlot)
      }
      else {
        showNotification("Launch differential expression analysis first",closeButton=TRUE,type="warning")
        suppressWarnings(plotly::plotly_empty())
      }

    })


    # Show MA plot of differential expression analysis
    output$MAplot <- plotly::renderPlotly({

      # show only if differential expression analysis was already done
      if("fold_change" %in% colnames(values$diffexp_summary_table)) {
        MAplot <- plot_MA(input_data =values$diffexp_summary_table,color_palette = input$col_palette,reverse_palette = input$reverse)
        plotly::ggplotly(MAplot)
      }
      else {
        showNotification("Launch differential expression analysis first",closeButton=TRUE,type="warning")
        plotly::plotly_empty()
      }

    })

    # Calculate differential expression if action button was clicked
    observeEvent(input$compute_diff_exp,
    {
      if(number_of_samples>1) {
        withProgress(message="Computing Differential epxression using binomial test...",
          detail = "This step can take a little while",
          value = 0.05,min=0,max=1,
          {

            values$polya_table_summarized <- summarize_polya(values$polya_table,summary_factors = c(input$groupingFactor))
            incProgress(4/20)
            values$polya_table_summarized_per_sample <- summarize_polya(values$polya_table,summary_factors = "sample_name")
            incProgress(4/20)

            values$diffexp_summary_table <- calculate_diff_exp_binom(values$polya_table_summarized,grouping_factor = input$groupingFactor,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,summarized_input=TRUE)
            incProgress(10/20)
            values$condition1 <- input$condition1_diff_exp
            values$condition2 <- input$condition2_diff_exp
            values$diff_exp_grouping_factor = input$groupingFactor
        })
      }
      else {
        showNotification("Differential analysis is not possible when only single sample is present",closeButton=TRUE,type="error")
      }
    })

    # show scatter plot of counts if action button was clicked
    observeEvent(input$show_scatter_plot,
    {
      if(number_of_samples>1) {
        withProgress(message="Computing data for scatter counts plot...",
          detail = "This step can take a little while",
          value = 0.05,min=0,max=1,
          {
            values$polya_table_summarized_scatter <- summarize_polya(values$polya_table,summary_factors = c(input$groupingFactor))
            incProgress(18/20)
            values$condition1_scatter <- input$condition1_diff_exp
            values$condition2_scatter <- input$condition2_diff_exp
            values$groupingFactor_scatter = input$groupingFactor
        })
      }
      else {
        showNotification("Scatter plot cannot be shown when only single sample is present",closeButton=TRUE,type="error")
      }
    })



    ## Basic QC Info section  ---------------------------------------------------------


    output$numberOfSamples <- renderInfoBox({
      infoBox(
        "Samples in the analysis:  ", number_of_samples, icon = icon("vials"),
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
        "Reads processed ", sum(values$processing_info_global$n), icon = icon("stream"),
        color = "black"
      )
    })

    output$numberOfPassReads <- renderInfoBox({
      infoBox(
        "Reads passing filters ", values$processing_info_global_spread$PASS, icon = icon("check"),
        color = "green"
      )
    })

    output$numberOfAdapterReads <- renderInfoBox({
      infoBox(
        "Adapter reads ", values$processing_info_global_spread$ADAPTER, icon = icon("times"),
        color = "red"
      )
    })

    output$numberOfReadLoadFailedReads <- renderInfoBox({
      infoBox(
        "Reads failed to load ", values$processing_info_global_spread$READ_FAILED_LOAD, icon = icon("times"),
        color = "red"
      )
    })

    output$numberOfNoregionReads <- renderInfoBox({
      infoBox(
        "NoRegion reads ", values$processing_info_global_spread$NOREGION, icon = icon("times"),
        color = "red"
      )
    })

    #Show Nanopolish QC summary plot
    output$nanopolish_qc_summary_plot <- plotly::renderPlotly({

      nanopolish_qc_summary_plot<-plot_nanopolish_qc(values$processing_info_per_sample,color_palette =  input$col_palette,reverse_palette = input$reverse, frequency=input$show_frequency_plot_nanopolioshQC)
      plotly::ggplotly(nanopolish_qc_summary_plot)

    })

  })

}
