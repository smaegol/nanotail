# TODO ShinyApp: Allow to choose columns shown in the table
# TODO ShinyApp: Allow launching on raw polya data.frame, calculate statistics using selected factor in the app
# TODO ShinyApp: Export report
# TODO allow statistics on more than 2 levels of factor (anova?)
# TODO Add withProgress

#' wrapper for NanoTail Shiny interface
#'
#' @param polya_table polyA predictions table. Can be obtained with \link{read_polya_multiple}
#' @param precomputed_polya_statistics precomputed differential adenylation table (obtained with \link{calculate_polya_stats})
#'
#' @export
#'

nanoTailApp <- function(polya_table,precomputed_polya_statistics=NA,precomputed_annotations=NA) {

  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("NanoTail requires 'shiny'. Please install it using
         install.packages('shiny')")
  }
  if ( !requireNamespace('shiny',quietly = TRUE) ) {
    stop("NanoTail requires 'shinydashboard'. Please install it using
         install.packages('shinydashboard')")
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


  colors_pal <- list()
  colors_pal$ggsci_based_colors <- c("npg", "aaas", "nejm", "lancet", "jama",
                                  "jco", "ucscgb", "d3", "locuszoom",
                                  "igv", "uchicago", "startrek", "tron",
                                  "futurama", "rickandmorty", "simpsons",
                                  "gsea")
  colors_pal$RColorBrewer_qualitative = rownames(RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category=="qual",])


  # Get all colors given a palette name(s)
  get_brewer_name <- function(name) {
    pals <- RColorBrewer::brewer.pal.info[rownames(RColorBrewer::brewer.pal.info) %in% name, ]
    if(nrow(pals)>0) {
      res <- lapply(
        X = seq_len(nrow(pals)),
        FUN = function(i) {
          RColorBrewer::brewer.pal(n = pals$maxcolors[i], name = rownames(pals)[i])
        }
      )
      unlist(res)
    }
    else {
        eval(parse(text = paste0("ggsci::pal_",name,"()")))(7)
    }

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

  colortext_pals <- rep(c("white", "black"), times = sapply(colors_pal, length))

  ### !<- colorpicker end ->! ###

  processing_info_global <- get_nanopolish_processing_info(polya_table,grouping_factor = NA)
  processing_info_global_spread <- processing_info_global %>% tidyr::spread(qc_tag,n)
  processing_info_per_sample <- get_nanopolish_processing_info(polya_table,grouping_factor = "sample_name")
  processing_info_per_sample_spread <- processing_info_per_sample %>% tidyr::spread(qc_tag,n)

  polya_summary_table <- summarize_polya(polya_data = polya_table_passed,summary_factors = "sample_name")
  if (!is.na(precomputed_polya_statistics)) {
    assertthat::assert_that(is.list(precomputed_polya_statistics),msg="Please provide the output of calculate_polya_stats() as precomputed_polya_statistics")
    assertthat::assert_that("summary" %in% names(precomputed_polya_statistics),msg="Please provide the output of calculate_polya_stats() as precomputed_polya_statistics")
    initial_summary_table = precomputed_polya_statistics$summary  %>% dplyr::select(transcript,dplyr::ends_with("_counts"),dplyr::ends_with("_polya_median"),p.value,padj)
    initial_table_for_volcano <- precomputed_polya_statistics$summary %>% dplyr::select(transcript,fold_change,padj,significance)
  }
  else {
    initial_summary_table <- polya_summary_table %>% dplyr::select(transcript,sample_name,polya_median) %>% tidyr::spread(sample_name,polya_median)
    initial_table_for_volcano <- initial_summary_table
  }
  initial_summary_table_counts <- polya_summary_table %>% dplyr::select(transcript,sample_name,counts) %>% tidyr::spread(sample_name,counts)

  grouping_factor_levels <- polya_table_passed %>% dplyr::select_if(is.factor) %>% colnames


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
        shinydashboard::menuItem("QC info", icon = shiny::icon("info-square"), tabName = "basicInfo"),
        shinydashboard::menuItem("Global polyA distribution", icon = shiny::icon("chart-line"), tabName = "global_distr"),
        shinydashboard::menuItem("Differential adenylation", icon = shiny::icon("dna"), tabName = "diff_polya"),
        shinydashboard::menuItem("Differential expression", icon = shiny::icon("dna"), tabName = "diff_exp"),
        shinydashboard::menuItem("Annotation_analysis", icon = shiny::icon("dna"), tabName = "annotation_analysis"),
        shinydashboard::menuItem("Plot settings", icon = shiny::icon("settings"), tabName = "plot_settings"),
        shinydashboard::menuItem("About", tabName = "dashboard", icon = shiny::icon("info"))
      ),
      shiny::selectInput("groupingFactor","Group by",choices=grouping_factor_levels),
      shiny::uiOutput("condition1UI"),
      shiny::uiOutput("condition2UI"),
      shiny::checkboxInput("plot_only_selected_conditions","Plot only selected conditions",value=FALSE)
    ),

    # dashboard body definition ---------------------------------------------------------
    shinydashboard::dashboardBody(
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "dashboard",
          shiny::fillPage(padding = 0, title = NULL, bootstrap = F, theme = NULL,
             shiny::wellPanel(style = "background-color: #ffffff; overflow-y:scroll; max-height: 750px;",
                       shiny::includeMarkdown(system.file("extdata", "about.md",package = "nanotail")))
          )),
        shinydashboard::tabItem(tabName = "basicInfo",
          shiny::fluidRow(
            shinydashboard::infoBoxOutput("numberOfSamples"),
            shinydashboard::infoBoxOutput("numberOfTranscripts"),
            shinydashboard::infoBoxOutput("numberOfAllReads")),
          shiny::fluidRow(
            shiny::column(1),
            shiny::column(7,shiny::fluidRow(shiny::h4("Per sample plot:")),
               shiny::fluidRow(
                  plotly::plotlyOutput("nanopolish_qc_summary_plot"),
                  shiny::checkboxInput("show_frequency_plot_nanopolioshQC",label = "Show QC as frequencies?",value = TRUE))),
            shiny::column(1),
            shiny::column(3,shiny::fluidRow(shiny::h4("Information for all samples:")),
               shiny::fluidRow(shinydashboard::infoBoxOutput("numberOfPassReads",width = NULL)),
               shiny::fluidRow(shinydashboard::infoBoxOutput("numberOfAdapterReads",width = NULL)),
               shiny::fluidRow(shinydashboard::infoBoxOutput("numberOfReadLoadFailedReads",width = NULL)),
               shiny::fluidRow(shinydashboard::infoBoxOutput("numberOfNoregionReads",width = NULL))))),
        shinydashboard::tabItem(tabName = "plot_settings",
          shinydashboard::box(
            shinyWidgets::pickerInput(inputId = "col_palette", label = "Colour scale",
                choices = colors_pal, selected = "Set1", width = "100%",
                choicesOpt = list(
                  content = sprintf(
                    "<div style='width:100%%;padding:5px;border-radius:4px;background:%s;color:%s'>%s</div>",
                    unname(background_pals), colortext_pals, names(background_pals)))),
            #shiny::uiOutput("Colorblind"),

            shiny::sliderInput("scale_limit_low","Plot scale limit low",0,1024,0,10),
            shiny::sliderInput("scale_limit_high","Plot scale limit high",50,1024,200,10),
            shiny::sliderInput("counts_scatter_low_value","Counts scatter lower counts limit",0,10024,0,50),
            shiny::sliderInput("counts_scatter_high_value","Counts scatter high counts limit",100,10024,1000,50),
            shiny::selectInput("center_values_for_distribution_plot",choices=c("none","median","mean","gm_mean"),selected = "median",label = "Show vertical lines with centered values on density plots?")
            )),
        shinydashboard::tabItem(tabName = "global_distr",
          shiny::h2("Global comparison"),
          shiny::fluidRow(
            #shinydashboard::box(plotly::plotlyOutput('polya_global') %>% shinycssloaders::withSpinner(type = 4),collapsible=TRUE))),
            shinydashboard::box(plotly::plotlyOutput('polya_global') %>% shinycssloaders::withSpinner(type = 4),collapsible=TRUE))),
        shinydashboard::tabItem(tabName = "diff_polya",
          shiny::fluidRow(
            shinydashboard::box(DT::dataTableOutput('diff_polya'),collapsible = TRUE,height = '90%'),
            shinydashboard::tabBox(title = "per transcript plots",id="tabset1", height="500px",
                   shiny::tabPanel("boxplot",plotly::plotlyOutput('polya_boxplot') %>% shinycssloaders::withSpinner(type = 4),
                                   shiny::checkboxInput("violin_instead_of_boxplot",value=FALSE,label = "Use violin plot instead of boxplot?")),
                   shiny::tabPanel("distribution plot",plotly::plotlyOutput('polya_distribution') %>% shinycssloaders::withSpinner(type = 4)
                                   ),
                   shiny::tabPanel("isoforms_boxplot",shiny::plotOutput('isoforms_boxplot') %>% shinycssloaders::withSpinner(type = 4)
                   ),
                   shiny::tabPanel("volcano plot polya",plotly::plotlyOutput('polya_volcano') %>% shinycssloaders::withSpinner(type = 4)))),
          shiny::fluidRow(
            shinydashboard::box(shiny::actionButton("compute_diff_polya",
                             shiny::HTML("Compute Polya statistics"),
                             icon = shiny::icon("spinner")),
                shiny::checkboxInput("use_dwell_time_for_statistics",label = "Use dwell time instead of calculated polya_length for statistics",value = FALSE),
                shiny::selectInput("polya_stat_test","statistical test to use",choices = c("Wilcoxon","KS"),selected = "Wilcoxon"),
                shiny::sliderInput("min_reads_for_polya_stats","Minimal count of reads per transcript",0,100,10,1),collapsible = TRUE,width=12))),
        shinydashboard::tabItem(tabName = "diff_exp",
                                shiny::fluidRow(
                                  shinydashboard::box(DT::dataTableOutput('diff_exp_table')),
                                  shinydashboard::tabBox(title = "plots",id="tabset2",height="500px",
                                                         #tabPanel("pca_biplot",plotOutput('pca_biplot') %>% shinycssloaders::withSpinner(type = 4)),
                                                         shiny::tabPanel("scatter plot of counts",shiny::actionButton("show_scatter_plot",
                                                                                                                      shiny::HTML("Show scatter plot"),
                                                                                                                      icon = shiny::icon("spinner")),plotly::plotlyOutput('counts_plot') %>% shinycssloaders::withSpinner(type = 4)),
                                                         shiny::tabPanel("volcano_plot",plotly::plotlyOutput('volcano') %>% shinycssloaders::withSpinner(type = 4)),
                                                         shiny::tabPanel("MA_plot",plotly::plotlyOutput('MAplot') %>% shinycssloaders::withSpinner(type = 4)))),
                                shiny::fluidRow(shiny::actionButton("compute_diff_exp",
                                                                    shiny::HTML("Compute Differential Expression using Binomial Test"),
                                                                    icon = shiny::icon("spinner")))),
        shinydashboard::tabItem(tabName = "annotation_analysis",
          shiny::fluidRow(
            shinydashboard::box(DT::dataTableOutput('annotation_table')),
            shinydashboard::tabBox(title = "annotation_plots",id="tabset3",height="500px",
                  #tabPanel("pca_biplot",plotOutput('pca_biplot') %>% shinycssloaders::withSpinner(type = 4)),
                  shiny::tabPanel("box_plot_of_annotations",plotly::plotlyOutput('annotations_box_plot') %>% shinycssloaders::withSpinner(type = 4),

                                  shiny::uiOutput("select_annotation_factor_levelsUI")
                                  ),
                  shiny::tabPanel("distribution_plot_of_annotations",plotly::plotlyOutput('annotation_distribution_plot') %>% shinycssloaders::withSpinner(type = 4)))),
          shiny::fluidRow(shiny::actionButton("get_annotables_annotation",
                                 shiny::HTML("Annotate input dataset"),
                                 icon = shiny::icon("spinner")),
                          shiny::selectInput("annotables_genome",label="Annotables genome:",choices = c("bdgp6","grch38","grcm38","grch37","rnor6","galgal5","wbcel235","mmul801"),selected="grch38",multiple = FALSE,selectize = TRUE),
                          shiny::uiOutput("select_annotation_factorUI")
                          )))))


  # call to shiny app and server-siude functions
  shiny::shinyApp(ui = nanotail_shiny_ui, server = function(input, output) {


    ## Reactive values definition  ---------------------------------------------------------
    values <- reactiveValues()
    values$polya_table <- polya_table_passed
    values$polya_statistics_summary_table = initial_summary_table
    values$polya_table_for_volcano = initial_table_for_volcano
    values$diffexp_summary_table = initial_summary_table_counts
    values$polya_table_summarized <- polya_summary_table
    values$processing_info_global <- processing_info_global
    values$processing_info_global_spread <- processing_info_global_spread
    values$processing_info_per_sample <- processing_info_per_sample
    values$processing_info_per_sample_spread <- processing_info_per_sample_spread
    values$diff_exp_grouping_factor = "sample_name"
    if(!is.na(precomputed_annotations)) {
      values$polya_table_annotables_annotated <- precomputed_annotations
    }
    else {
      values$polya_table_annotables_annotated <- polya_table_passed %>% dplyr::mutate(annotation_group=factor(gsub("^(.).*","\\1",transcript)))
    }
    # get polyA data for currently selected transcript
    data_transcript <- shiny::reactive({
      summary_table = values$polya_statistics_summary_table

      selected_row <- input$diff_polya_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript = subset(polya_table_passed, transcript==selected_transcript)
      data_transcript <- data_transcript %>% dplyr::mutate(replicate = gsub(".*(.)$","\\1", sample_name))

      data_transcript
    })

    # get polyA data for currently selected transcript
    data_transcript_annot <- shiny::reactive({
      summary_table = values$polya_statistics_summary_table

      selected_row <- input$diff_polya_rows_selected
      selected_transcript = summary_table[selected_row,]$transcript
      data_transcript_annot = subset(values$polya_table_annotables_annotated, transcript==selected_transcript)

      data_transcript_annot
    })

    # UI elements rendering ---------------------------------------------------------

    ## Generate select  inputs for conditions used for comparison, based on selected groupingFactor:  ---------------------------------------------------------
    output$condition1UI <- shiny::renderUI({
      group_factor = input$groupingFactor
      shiny::selectizeInput("condition1_diff_exp","Condition 1",choices=levels(polya_table_passed[[group_factor]]),selected = levels(polya_table_passed[[group_factor]])[1])

    })
    output$condition2UI <- shiny::renderUI({
      group_factor = input$groupingFactor
      shiny::selectizeInput("condition2_diff_exp","Condition 2",choices=levels(polya_table_passed[[group_factor]]),selected = levels(polya_table_passed[[group_factor]])[2])

    })


    output$Colorblind <- shiny::renderUI({
      if(RColorBrewer::brewer.pal.info %>% subset(rownames(.) == input$col_palette) %>% .$colorblind) {
        txt_ <- "Yes"
      } else {
        txt_ <- "No"
      }
      shiny::HTML(paste("<p>Is colour scale colour blind friendly?<b>", txt_, "</b></p>"))
    })



    # Output elements rendering section  ---------------------------------------------------------

    ## Use DT::datatable for tables
    output$diff_polya = DT::renderDataTable(values$polya_statistics_summary_table %>% dplyr::rename_all(dplyr::funs(stringr::str_replace_all(.,"_"," "))), server = TRUE, selection=list(mode = 'single',selected = 1,target = "row"),options = list(dom = 'ftip'))
    output$diff_exp_table = DT::renderDataTable(values$diffexp_summary_table %>% dplyr::rename_all(dplyr::funs(stringr::str_replace_all(.,"_"," "))), server = TRUE, selection=list(mode = 'single',selected = 1,target = "row"))
    output$annotation_table = DT::renderDataTable(data_annotation(), server = TRUE, selection=list(mode = 'multiple',selected = 1,target = "row"))

    annotation_proxy <- DT::dataTableProxy('annotation_table')

    ## Show global polyA distribution as density plot
    output$polya_global = plotly::renderPlotly({

      if (input$plot_only_selected_conditions) {
        global_distribution_plot <- plot_polya_distribution(polya_data = polya_table_passed,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,show_center_values=input$center_values_for_distribution_plot)
      }
      else {
        global_distribution_plot <- plot_polya_distribution(polya_data = polya_table_passed,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,show_center_values=input$center_values_for_distribution_plot)
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
        print(data_transcript)
        if (input$plot_only_selected_conditions) {
          polya_boxplot <- plot_polya_boxplot(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_y_limit_low = input$scale_limit_low,scale_y_limit_high = input$scale_limit_high,color_palette = input$col_palette,plot_title = selected_transcript,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp, violin=input$violin_instead_of_boxplot)
        }
        else {
          polya_boxplot <- plot_polya_boxplot(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_y_limit_low = input$scale_limit_low,scale_y_limit_high = input$scale_limit_high,color_palette = input$col_palette,plot_title = selected_transcript, violin=input$violin_instead_of_boxplot)
        }

        plotly::ggplotly(polya_boxplot)
      }
      else {
        suppressWarnings(plotly::plotly_empty())
      }
    })


    # Show boxplot of estimated polya lengths for selected transcript
    output$isoforms_boxplot = shiny::renderPlot({



      req(values$polya_table_annotables_annotated)

      summary_table = values$polya_statistics_summary_table

      if (length(input$diff_polya_rows_selected)>0) {
        selected_row <- input$diff_polya_rows_selected
        selected_transcript = summary_table[selected_row,]$transcript

        data_transcript = data_transcript_annot()
        print(selected_transcript)
        print(data_transcript)
        polya_boxplot <- plot_polya_boxplot(polya_data = data_transcript,groupingFactor = "ensembl_transcript_id_short",scale_y_limit_low = input$scale_limit_low,scale_y_limit_high = input$scale_limit_high,color_palette = input$col_palette,plot_title = selected_transcript, violin=input$violin_instead_of_boxplot,additional_grouping_factor = input$groupingFactor) + theme(axis.text.x = element_text(angle=90))
        print(polya_boxplot)
        #plotly::ggplotly(polya_boxplot)
      }
      else {
        print("no data")
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
            transcript_distribution_plot <- plot_polya_distribution(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette, plot_title = selected_transcript,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,show_center_values=input$center_values_for_distribution_plot)
          }
          else {
            transcript_distribution_plot <- plot_polya_distribution(polya_data = data_transcript,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette, plot_title = selected_transcript,show_center_values=input$center_values_for_distribution_plot)
          }
          plotly::ggplotly(transcript_distribution_plot)
        }
        else {
          suppressWarnings(plotly::plotly_empty())
        }

    })

    # Show volcano plot of differential adenylation analysis
    output$polya_volcano <- plotly::renderPlotly({

      # show only if differential adenylation analysis was already done
      if("fold_change" %in% colnames(values$polya_table_for_volcano)) {
        volcanoPlot <- plot_volcano(input_data =values$polya_table_for_volcano,color_palette = input$col_palette)
        plotly::ggplotly(volcanoPlot)
      }
      else {
        shiny::showNotification("Launch differential polyadenylation analysis first",closeButton=TRUE,type="warning")
        suppressWarnings(plotly::plotly_empty())
      }

    })

    # Calculate differential adenylation if action button was clicked
    shiny::observeEvent(input$compute_diff_polya,
    {
      if(number_of_samples>1) {
        shiny::withProgress(message="Computing polya statistics...",
          detail = "This step can take a little while",
          value = 0.05,min=0,max=1,
          {
            condition1 = input$condition1_diff_exp
            condition2 = input$condition2_diff_exp
            shiny::incProgress(1/20)
            polya_stats <- calculate_polya_stats(values$polya_table,grouping_factor = input$groupingFactor,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,stat_test=input$polya_stat_test,use_dwell_time=input$use_dwell_time_for_statistics,min_reads = input$min_reads_for_polya_stats,add_summary = TRUE,length_summary_to_show = "median")
            shiny::incProgress(15/20)
            #values$polya_statistics_summary_table <- polya_stats %>% dplyr::select(transcript,!!rlang::sym(paste0(condition1,"_counts")),!!rlang::sym(paste0(condition2,"_counts")),!!rlang::sym(paste0(condition1,"_polya_median")),!!rlang::sym(paste0(condition2,"_polya_median")),p.value,padj)
            values$polya_statistics_summary_table <- polya_stats
            values$polya_table_for_volcano <- polya_stats %>% dplyr::select(transcript,fold_change,padj,significance)
            values$condition1 <- input$condition1_diff_exp
            values$condition2 <- input$condition2_diff_exp
            shiny::incProgress(3/20)
        })
      }
      else {
        shiny::showNotification("Differential analysis is not possible when only single sample is present",closeButton=TRUE,type="error")
      }
    })



    ## Differential expression analysis section  ---------------------------------------------------------

    # Show scatter plot of raw counts
    output$counts_plot <- plotly::renderPlotly({

      counts_scatter_plot <- plot_counts_scatter(polya_data = values$polya_table_summarized_scatter,groupingFactor = values$groupingFactor_scatter,color_palette = input$col_palette,condition1 = values$condition1_scatter, condition2 = values$condition2_scatter, min_counts = input$counts_scatter_low_value, max_counts = input$counts_scatter_high_value)
      plotly::ggplotly(counts_scatter_plot)


    })

    # Show volcano plot of differential expression analysis
    output$volcano <- plotly::renderPlotly({

      # show only if differential expression analysis was already done
      if("fold_change" %in% colnames(values$diffexp_summary_table)) {
        volcanoPlot <- plot_volcano(input_data =values$diffexp_summary_table,color_palette = input$col_palette)
        plotly::ggplotly(volcanoPlot)
      }
      else {
        shiny::showNotification("Launch differential expression analysis first",closeButton=TRUE,type="warning")
        suppressWarnings(plotly::plotly_empty())
      }

    })


    # Show MA plot of differential expression analysis
    output$MAplot <- plotly::renderPlotly({

      # show only if differential expression analysis was already done
      if("fold_change" %in% colnames(values$diffexp_summary_table)) {
        MAplot <- plot_MA(input_data =values$diffexp_summary_table,color_palette = input$col_palette)
        plotly::ggplotly(MAplot)
      }
      else {
        shiny::showNotification("Launch differential expression analysis first",closeButton=TRUE,type="warning")
        plotly::plotly_empty()
      }

    })

    # Calculate differential expression if action button was clicked
    shiny::observeEvent(input$compute_diff_exp,
    {
      if(number_of_samples>1) {
        shiny::withProgress(message="Computing Differential epxression using binomial test...",
          detail = "This step can take a little while",
          value = 0.05,min=0,max=1,
          {

            values$polya_table_summarized <- summarize_polya(values$polya_table,summary_factors = c(input$groupingFactor))
            shiny::incProgress(4/20)
            values$polya_table_summarized_per_sample <- summarize_polya(values$polya_table,summary_factors = "sample_name")
            shiny::incProgress(4/20)

            values$diffexp_summary_table <- calculate_diff_exp_binom(values$polya_table_summarized,grouping_factor = input$groupingFactor,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,summarized_input=TRUE)
            shiny::incProgress(10/20)
            values$condition1 <- input$condition1_diff_exp
            values$condition2 <- input$condition2_diff_exp
            values$diff_exp_grouping_factor = input$groupingFactor
        })
      }
      else {
        shiny::showNotification("Differential analysis is not possible when only single sample is present",closeButton=TRUE,type="error")
      }
    })

    # show scatter plot of counts if action button was clicked
    shiny::observeEvent(input$show_scatter_plot,
    {
      if(number_of_samples>1) {
        shiny::withProgress(message="Computing data for scatter counts plot...",
          detail = "This step can take a little while",
          value = 0.05,min=0,max=1,
          {
            values$polya_table_summarized_scatter <- summarize_polya(values$polya_table,summary_factors = c(input$groupingFactor))
            shiny::incProgress(18/20)
            values$condition1_scatter <- input$condition1_diff_exp
            values$condition2_scatter <- input$condition2_diff_exp
            values$groupingFactor_scatter = input$groupingFactor
        })
      }
      else {
        shiny::showNotification("Scatter plot cannot be shown when only single sample is present",closeButton=TRUE,type="error")
      }
    })



    ## Basic QC Info section  ---------------------------------------------------------


    output$numberOfSamples <- shinydashboard::renderInfoBox({
      shinydashboard::infoBox(
        "Samples in the analysis:  ", number_of_samples, icon = shiny::icon("vials"),
        color = "yellow"
      )
    })

    output$numberOfTranscripts <- shinydashboard::renderInfoBox({
      shinydashboard::infoBox(
        "Transcripts analyzed ", number_of_transcripts, icon = shiny::icon("dna"),
        color = "blue"
      )
    })

    output$numberOfAllReads <- shinydashboard::renderInfoBox({
      shinydashboard::infoBox(
        "Reads processed ", sum(values$processing_info_global$n), icon = shiny::icon("stream"),
        color = "black"
      )
    })

    output$numberOfPassReads <- shinydashboard::renderInfoBox({
      shinydashboard::infoBox(
        "Reads passing filters ", values$processing_info_global_spread$PASS, icon = shiny::icon("check"),
        color = "green"
      )
    })

    output$numberOfAdapterReads <- shinydashboard::renderInfoBox({
      shinydashboard::infoBox(
        "Adapter reads ", values$processing_info_global_spread$ADAPTER, icon = shiny::icon("times"),
        color = "red"
      )
    })

    output$numberOfReadLoadFailedReads <- shinydashboard::renderInfoBox({
      shinydashboard::infoBox(
        "Reads failed to load ", values$processing_info_global_spread$READ_FAILED_LOAD, icon = shiny::icon("times"),
        color = "red"
      )
    })

    output$numberOfNoregionReads <- shinydashboard::renderInfoBox({
      shinydashboard::infoBox(
        "NoRegion reads ", values$processing_info_global_spread$NOREGION, icon = shiny::icon("times"),
        color = "red"
      )
    })

    #Show Nanopolish QC summary plot
    output$nanopolish_qc_summary_plot <- plotly::renderPlotly({

      nanopolish_qc_summary_plot<-plot_nanopolish_qc(values$processing_info_per_sample,color_palette =  input$col_palette, frequency=input$show_frequency_plot_nanopolioshQC)
      plotly::ggplotly(nanopolish_qc_summary_plot)

    })



    ## Annotation analysis section ---------------------------------------------------------


      data_selected_annotation <- shiny::reactive({
        #summary_table = values$polya_table_annotated

        selected_row <- input$annotation_table_rows_selected
        data_annotation = data_annotation()
        data_annotation = as.data.frame(data_annotation)
        column_name = input$annotation_factor
        print(paste0("column name: ",column_name))
        selected_annotation = as.character(data_annotation[selected_row,column_name])
        print("selected annotation:")
        print(selected_annotation)
        data_transcript = subset(values$polya_table_annotables_annotated, eval(parse(text=column_name)) %in% c(selected_annotation))
        data_transcript <- data_transcript %>% dplyr::group_by(.dots = c("read_id",input$groupingFactor)) %>% dplyr::slice(1)
        print("nrow data transc:")
        print(nrow(data_transcript))
        data_transcript
      })


      data_annotation <- shiny::reactive({

        message(input$annotation_factor)

        if (!is.null(input$annotation_factor)) {
        data_annotation_summarized <- values$polya_table_annotables_annotated %>% dplyr::group_by(.dots = c(input$annotation_factor,input$groupingFactor)) %>% dplyr::summarize(polya_gm_mean=gm_mean(polya_length)) %>%
          tidyr::spread(!!rlang::sym(input$groupingFactor),polya_gm_mean)


        data_annotation_summarized

        }
        else {
          shiny::showNotification(ui = "Annotation factor is empty",closeButton = TRUE,type = "message")
        }


      })


      ## Generate selectinputs for annotations:
      output$select_annotation_factorUI <- shiny::renderUI({
        annotation_factor_choices=values$polya_table_annotables_annotated %>% dplyr::select_if(is.factor) %>% dplyr::select(-grouping_factor_levels) %>% colnames
        shiny::selectizeInput("annotation_factor","Column used for annotation comparisons",choices=annotation_factor_choices,selected = annotation_factor_choices[1])

      })
      output$select_annotation_factor_levelsUI <- shiny::renderUI({
        annotation_factor_selected = input$annotation_factor
        annotation_factor_possible_levels = levels(values$polya_table_annotables_annotated[[annotation_factor_selected]])
        shiny::selectInput(inputId = "annotation_factor_levels",label="Which annotations to show on the plot",choices=annotation_factor_possible_levels,selected = annotation_factor_possible_levels,multiple = TRUE)

      })


      shiny::observeEvent(input$annotation_factor,{
        annotation_proxy %>% DT::selectRows(NULL)

      })

    shiny::observeEvent(input$get_annotables_annotation,
                      {

                          shiny::withProgress(message="Annotating input polyA table with annotables...",
                                              detail = "This step can take a little while",
                                              value = 0.05,min=0,max=1,
                                              {

                                                values$polya_table_annotables_annotated <- annotate_with_annotables(values$polya_table,genome = input$annotables_genome)


                                                values$polya_table_annotables_annotated <-data.table::data.table(values$polya_table_annotables_annotated)
                                                shiny::incProgress(18/20)
                                                data.table::setkeyv(values$polya_table_annotables_annotated,c("transcript","biotype","strand","chr"))
                                                shiny::incProgress(1/20)
                                              })

                      })

    # Show denisty plot of estimated polya lengths for selected transcript
    output$annotation_distribution_plot = plotly::renderPlotly({


    #print(input$annotation_table_rows_selected)
      #if (length(input$annotation_table_rows_selected)>0) {
        selected_row <- input$annotation_table_rows_selected
        print(selected_row)
        data_transcript_annotated = data_selected_annotation()
        #print(data_transcript_annotated)
        if (input$plot_only_selected_conditions) {
          transcript_distribution_plot <- plot_polya_distribution(polya_data = data_transcript_annotated,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,condition1 = input$condition1_diff_exp,condition2 = input$condition2_diff_exp,show_center_values=input$center_values_for_distribution_plot)
        }
        else {
          transcript_distribution_plot <- plot_polya_distribution(polya_data = data_transcript_annotated,groupingFactor = input$groupingFactor,scale_x_limit_low = input$scale_limit_low,scale_x_limit_high = input$scale_limit_high,color_palette = input$col_palette,show_center_values=input$center_values_for_distribution_plot)
        }
        plotly::ggplotly(transcript_distribution_plot)
     # }
     # else {
     #   suppressWarnings(plotly::plotly_empty())
     # }
    })


    output$annotations_box_plot = plotly::renderPlotly({

      print(paste0("annot_factor: ",input$annotation_factor))
      print(input$annotation_factor_levels)
        transcript_distribution_plot <- plot_annotations_comparison_boxplot(annotated_polya_data = values$polya_table_annotables_annotated,grouping_factor = input$groupingFactor,annotation_factor = input$annotation_factor, annotation_levels = input$annotation_factor_levels, scale_y_limit_low = input$scale_limit_low,scale_y_limit_high = input$scale_limit_high,color_palette = input$col_palette,violin=FALSE)
        plotly::ggplotly(transcript_distribution_plot) %>% plotly::layout(boxmode = "group")

    })

  })


    }

