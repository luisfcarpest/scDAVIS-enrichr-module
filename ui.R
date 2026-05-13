library(shiny)
library(plotly)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyFiles)
library(shinyjs)
library(shinyhttr)
library(shinycssloaders)
library(esquisse)
library(colourpicker)
library(gtools)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(Matrix)
library(DT)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(ggplot2)
library(dittoSeq)

library(dashboardthemes)
library(shinyalert)

source("R/helperFunctions.R")

# #############################################################################################
# CARGAMOS EL ARCHIVO DEL MODULO

source("R/module_enrichR.R")

# #############################################################################################

shinyUI(
  dashboardPage(
    # skin = "black-light",
    title = "scDAVIS: single-cell Data Analysis VISualization"
    ,header = dashboardHeader(
      title = 
        tagList(
          fluidRow(
            column(1,img(src = "scDAVISIcon.png"
                         ,contentType = "image/png"
                         ,width = "20px", height = "20px"))
            ,column(6,span(class = "logo-lg","scDAVIS"))
          )
        )
      , titleWidth = 250
      , leftUi = tagList(h4(textOutput("t1"),style="color:white; font-size:20px; font-family:Helvetica,Arial,sans-serif; margin-top: 4.5px;"))
    ),
    
    # SIDEBAR MENU
    sidebar = dashboardSidebar(
      useShinyjs(),
      width = 250
      ,sidebarMenu(
        
        id="sidetabs",
        
        menuItem(
          "Home"
          , tabName = "Home"
          , icon = icon("home"))
        , menuItem("Upload Seurat Object"
                   , icon = icon("upload")
                   , tabName = "QCcheck"
                   , fluidRow(
                     column(width = 12
                            , fileInput(
                              inputId = "file"
                              , label = div(h5("From Local File"), style = "color: white;")
                              , accept = '.rds')
                            )
                   )
                   , fluidRow(
                     column(width = 9
                            ,style='padding-right:0px;'
                            ,textInput("getURL",
                                       label = h5( "From URL", style = "color: white;")
                                       , value = "", width = "100%")
                     )
                     ,column(width = 3
                             ,style='margin-left:-48px; padding-top: 47px;'
                             ,align = 'left'
                             ,actionButton(
                               inputId = "GET",
                               label = "GET"
                             )
                     )
                   )
                   ,fluidRow(
                     column(width = 10
                            ,style='margin-left:14px; padding-top:0px;'
                            ,shinyWidgets::progressBar(
                              id = "upb"
                              , value = 50
                              , status = "info"
                              , commas = F
                              , display_pct = TRUE
                              , striped = TRUE)
                     )
                   )
        )
        , menuItem("Load Public Dataset"
                   , tabName = "Public"
                   , icon = icon("database")
        )
        ### DOWNLOAD ANALYSIS ###
        , disabled(
          menuItem("Download Analysis"
                   , icon = icon("download")
                   , id = "data_download"
                   , data_download_ui("module_data_download")
          )
        )
        
        ### MENU HELP
        , menuItem("QC Stats", tabName = "QCcheck", icon = icon("check"))
        , menuItem("Plots", tabName = "Plots", icon = icon("chart-line"))
        , menuItem("Help", tabName = "Help"
                   , icon = icon("question-circle")
        )
      )
    )
    
    ### Plot settings (right bar)
    ,controlbar = dashboardControlbar(
      skin = "light"
      , overlay = FALSE
      , collapsed = F
      , width = 260
      , controlbarMenu(
        controlbarItem(id = 1
                       , icon = icon("desktop")
                       , active = TRUE
                       , p()
                       # switch representation
                       , selectInput("selectRepresentation", label = h5("Dimensional reduction")
                                     , choices = NULL
                                     , selected = NULL)
                       # switch assay
                       , selectInput("select_assay", label = h5("Assay")
                                     , choices = NULL
                                     , selected = NULL)
                       , hr()
                       , fluidRow(
                         column(12,align="left"
                                , radioButtons("radioSelectClust", label = h5("Variable type")
                                               , choices = list("Categorical" = 0, "Numerical" = 1, "Gene" = 2, "Density" = 3)
                                               , selected = 0)
                         )
                       )
                       , hr()
                       , selectInput("selectCategorical", label = h5("Select Categorical Variable")
                                     , choices = NULL
                                     , selected = NULL, multiple = T)
                       , palettePicker(
                         inputId = "categoricalPalette", 
                         label = tags$p("Categorical Colors", style = "font-size: 100%;font-weight: normal;"),
                         choices = getDiscretePals(n=8),
                         selected = "Paired",
                         textColor = c(
                           rep("white", 5), rep("black", 4) 
                         )
                       )
                       , fluidRow(
                         column(5, checkboxInput(inputId = "revCatPal",label = "reverse",value = F)),
                         column(5,disabled(div(id="color_button",style = "position:absolute;right:3em;"
                                               ,dropdownButton(size = "sm"
                                                               , icon = icon("fill-drip"), up = F, right = T, circle = F, width = 125
                                                               , selectInput("color_var", label = h5("Select Category")
                                                                             , choices = NULL
                                                                             , selected = NULL)
                                                               , colourpicker::colourInput("color_cat", "Select color", value = "white", palette = "limited"))
                                               
                         ))
                         )
                         
                       )
                       , hr()
                       , selectInput("selectNumerical", label = h5("Select Numerical Variable")
                                     , choices = NULL
                                     , selected = NULL, multiple = T)
                       , palettePicker(
                         inputId = "numericalPalette", 
                         label = tags$p("Numerical Palette", style = "font-size: 100%;font-weight: normal;"),
                         choices = getContinuousPals(n=10), 
                         selected = "BYR",
                         textColor = c(
                           rep("white", 5), rep("black", 4) 
                         )
                       )
                       , checkboxInput(inputId = "revNumPal",label = "reverse",value = F)
                       , hr()
                       , selectizeInput("selectGene", label = h5("Select Gene")
                                        , choices = NULL
                                        , selected = NULL
                                        , multiple = T)
                       , palettePicker(
                         inputId = "genePalette", 
                         label = tags$p("Gene Palette", style = "font-size: 100%;font-weight: normal;"),
                         choices = getContinuousPals(n=10),
                         selected = "BYR",
                         textColor = c(
                           rep("white", 5), rep("black", 4) 
                         )
                       )
                       , checkboxInput(inputId = "revGenPal",label = "reverse",value = F)
                       , hr()
                       , sliderInput("selectDensity", label = h5("Select Density Width")
                                     , min = 0.5, max = 50, value = 2, step = 0.5
                       )
                       , palettePicker(
                         inputId = "densityPalette", 
                         label = tags$p("Density Palette", style = "font-size: 100%;font-weight: normal;"),
                         choices = getContinuousPals(n=10), 
                         selected = "BYR",
                         textColor = c(
                           rep("white", 5), rep("black", 4) 
                         )
                       )
                       , checkboxInput(inputId = "revDensPal",label = "reverse",value = F)
                       
        )
        , controlbarItem(id = 2
                         , icon = icon("adjust")
                         , p()
                         , data_contrasts_ui("module_data_contrasts")
        )
        , controlbarItem(id = 3
                         , icon = icon("filter")
                         , p()
                         , data_filter_ui("module_data_filter")
                         
        )
        , controlbarItem(id = 4
                         , icon = icon("edit")
                         , p()
                         , h5("Selected Cells")
                         , textOutput("brushCluster")
                         , br()
                         , actionButton("gateCells", label = "Gate To Cells")
                         , actionButton("undoGateCells", label = "Undo")
                         , hr()
                         , textInput("nameCellCluster", label = h5("Name of Cluster"), value = "")
                         , textOutput("cwarnings") #Cluster Warnings
                         , actionButton("saveCellCluster", label = "Save Cluster")
                         , br()
                         , useShinyjs()
                         , DT::dataTableOutput("cellCluster")
                         , textOutput('cellClust')
                         , br()
                         , hr()
                         , textInput("nameNewClustering", label = h5("Name of New Clustering Set"), value = "")
                         , textOutput("cswarnings") #Clustering Set Warnings
                         , div(id = "clusterControl"
                               , mainPanel(
                                 radioButtons("radioNewClusterControl", label = NULL
                                              , choices = list("Yes" = 0, "No" = 1)
                                              , selected = 0)
                               )
                         ) %>% shinyjs::hidden() 
                         , div(id = "invisibleVariable"
                               , mainPanel(
                                 radioButtons("invvar", label = NULL
                                              , choices = list("?" = 0, "Continue" = 1) 
                                              , selected = 0)
                               )
                         ) %>% shinyjs::hidden() 
                         , actionButton("saveNewClustering", label = "Save Manual Clustering")
                         , checkboxInput("checkboxKeepIdentCluster", label = "Keep old ident for not labelled cells", value = F)
                         , hr()
        )
      )
    )
    
    ### OVERVIEW INFO PANEL
    , body = dashboardBody(
      # includeScript("keycloak.js"),
      # includeScript("kcjs.js"),
      # tags$script("const kc = initKeycloak();"),
      tabItems(
        ### ABOUT TAB
        tabItem("Home",
                div(
                  id = "aboutDiv",
                  fluidPage(
                    fluidRow(
                      box(
                        width = 12
                        , title = h2(HTML(
                          '<span style="color:#3c72b3; font-weight: bold">s</span>ingle-<span
                              style="color:#3c72b3; font-weight: bold">c</span>ell <span
                              style="color:#3c72b3; font-weight: bold">D</span>ata <span
                              style="color:#3c72b3; font-weight: bold">A</span>nalysis and <span
                              style="color:#3c72b3; font-weight: bold">VIS</span>ualization <span
                              style="color:#3c72b3; font-weight: bold"> scDAVIS </span>'
                        )
                        )
                        , br()
                        , br()
                        , br()
                        , fluidRow(
                          column(width = 3,
                                 img(src="SpatialEmbryoHeart.png",width="100%"))
                          ,column(width = 7,
                                  p("IMPaCT-Data es el programa de IMPaCT que persigue apoyar el desarrollo de un sistema común, interoperable e integrado, de recogida y análisis de datos clínicos y moleculares aportando para ello el conocimiento y los recursos disponibles en el Sistema Español de Ciencia y Tecnología. Este desarrollo permitirá dar respuesta a preguntas de investigación a partir de los diferentes sistemas de información clínica y molecular disponibles. Fundamentalmente, persigue que los investigadores puedan disponer de una perspectiva poblacional basada en datos individuales.", style = "font-size:18px;") 
                          )
                        )
                      )
                      ,box(width = 12,title = NULL,collapsible = F,headerBorder = F
                           , fluidRow(
                             column(width = 4," ")
                             , column(width = 1,
                                      img(src="scRNA_icon.png",width="80%",align = "center"))
                             , column(width = 1,
                                      img(src="Imaris_icon.png",width="80%",align = "center"))
                             , column(width = 1,
                                      img(src="Cyt_icon.png",width="80%",align = "center"))
                             , column(width = 1,
                                      img(src="CSV_icon.png",width="80%",align = "center"))
                             , column(width = 4," ")
                             , br()
                           )
                           ,br()
                           ,fluidRow(
                             column(width = 2)
                             ,column(width = 8
                                     , p("scDAVIS is a web-based tool for the analysis and visualization of different modalities of single-cell omics data:",align = "center",style = "font-size:18px;")
                                     , p("IMPACT scDAVIS is coupled with access control for IMPACT users only. Contact IMPACT CNIC Team to get access.",style = "font-size:18px;",align="center")
                             )
                             ,column(width = 2)
                           )
                           , br()
                      )
                      , br()
                      ,box(width = 12,title = NULL,collapsible = F,headerBorder = F
                           , p("In scDAVIS you can (1) Upload previously processed data. (2) Load published analysis (currently, 4 public datasets), searchable through a keyword-based search engine. (3) Download analysis, exporting data for future reanalysis. (4) QC-Stats with basic information about cells and features/genes profiled (5) Plots for visualization of results: Dim-Plots, Violin-Plots, Bar-Plots, Heatmaps, Scatter-Plots and Dot-Plots. (6) Finally, scDAVIS incorporates interactive tools for results representation, statistical contrasts, data filtering, and manual annotation/correction.",align = "center",style = "font-size:18px;")
                      )
                      ,br()
                      ,box(width = 12,title = NULL,collapsible = F,headerBorder = F
                           ,fluidRow(
                             column(width = 2,
                                    img(src="data-repository-icon-5.png",width = "80%",align="right"))
                             ,column(width = 8,
                                     br(),
                                     br()
                                     , p("Importantly, scDAVIS will be also a repository of publicly available and IMPACT resource single-cell experiments, allowing further exploration and reuse of the data, in line with the compliance of FAIR principles. Also, it can serve as a repository to share single-cell data with the community and to ease integration of different datasets.",style = "font-size:18px;",align="left")
                                     , p("To permanently add a dataset to IMPACT scDAVIS repository and share it with other project members contact IMPACT CNIC Team for instructions.",style = "font-size:18px;",align="left")
                             )
                           )
                      )
                      , br()
                      , br()
                      ,box(width = 12,title = NULL,collapsible = F,headerBorder = F
                           , p("This tool was created by the CNIC Bioinformatics Unit",style = "font-size:18px;",align="left")
                      )
                    )
                  )
                )
        )
        ,tabItem("QCcheck",
                 qc_check_ui("module_qc_check")
        )
        ,tabItem("Public",
                 public_data_ui("module_public_data")
        )
        ,tabItem("Plots", 
                 fluidPage(
                   box(
                     title = "Plots"
                     ,width = 12
                     ,collapsible = TRUE
                     ,closable = F
                     ,fluidRow(
                       tabBox(title = "", id="plot_tab"
                              , width = 12
                              
                              ### DIMPLOTS TAB
                              , tabPanel(title = "DimPlot",
                                        plot_reductions_ui("module_plot_reductions"))
                              
                              ### VLN PLOTS TAB
                              ,tabPanel("Violin Plots",
                                        plot_vln_ui("module_plot_vln")
                              )
                              
                              ### BARPLOTS TAB (WIP)
                              ,tabPanel("BarPlots",
                                        plot_bar_ui("module_plot_bar")
                              )
                              
                              ### HEATMAP TAB
                              ,tabPanel("HeatMaps",
                                        plot_heatmap_ui("module_plot_heatmap")
                              )
                              
                              
                              ### SCATTERPLOTS TAB (WIP)
                              ,tabPanel("ScatterPlots",
                                        plot_scatter_ui("module_plot_scatter"))
                              
                              
                              ### DOTPLOT TAB (WIP)
                              , tabPanel("DotPlots"
                                         ,plot_dot_ui("module_plot_dot"))
                              
                              ####################################################################################
                              
                              # TABPANEL DE ENRICHR
                              
                              , tabPanel("EnrichR"
                                         ,enrichR_ui("module_enrichr"))
                              
                              ####################################################################################
                              
                       )
                     )
                     
                   )
                   ,box(
                     title = "Diff Expression Table"
                     ,width = 12
                     ,collapsible = F
                     ,collapsed = F
                     ,closable = F
                     ,box(
                       collapsible = F,
                       closable = F,
                       width = 12
                       , fluidRow(
                         div(style = "position:absolute;right:6.9em;"
                             , downloadButton(outputId = "downloadDataMks"
                                              , label = ""
                                              , class = "btn-sm"
                             )
                         )
                         ,div(style = "position:absolute;right:3em;"
                              , dropdownButton(
                                icon = icon(name = "upload"), size = "sm", up = F, right = T, circle = F, width = 4
                                , fileInput(inputId = "uploadTableFile"
                                            , label = "Upload a Table"
                                            , width= '100%'
                                            , accept = c('txt','csv','tsv')
                                )
                              )
                         )
                         ,div(style = "position:absolute;right:10em;"
                              , actionButton("clearSelect"
                                             , up = F, right = F
                                             , label = "Clear Selection"
                                             , class = "btn-sm")
                         )
                       )
                       , br()
                       , br()
                       , DT::dataTableOutput('df.mks')
                     )
                   )
                 )
        )
        ### HELP TAB
        ,tabItem("Help",
                 fluidPage(
                   fluidRow(
                     box(
                       h3("Session Info:")
                       , width = 8
                       , br()
                       , htmlOutput("sessionInfo")
                     )
                   )
                 )
        )
      )
    )
  )
)
