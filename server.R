# Library loading
# library(jsonlite)
suppressWarnings(suppressMessages(library(shiny)))
suppressWarnings(suppressMessages(library(shinyjs)))
suppressWarnings(suppressMessages(library(shinyFiles)))
suppressWarnings(suppressMessages(library(shinyWidgets)))
suppressWarnings(suppressMessages(library(shinycssloaders)))
suppressWarnings(suppressMessages(library(shinydashboard)))
suppressWarnings(suppressMessages(library(shinydashboardPlus)))
suppressWarnings(suppressMessages(library(DT)))
suppressWarnings(suppressMessages(library(htmlTable)))
suppressWarnings(suppressMessages(library(colourpicker)))
suppressWarnings(suppressMessages(library(gtools)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(SeuratDisk)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(plotly)))
suppressWarnings(suppressMessages(library(ComplexHeatmap)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(circlize)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(dittoSeq)))
suppressWarnings(suppressMessages(library(tidyverse)))

suppressWarnings(suppressMessages(library(Seurat)))

suppressWarnings(suppressMessages(library(ggpointdensity)))

##########################################################################################

suppressWarnings(suppressMessages(library(enrichR)))
suppressWarnings(suppressMessages(library(purrr)))

##########################################################################################

options(shiny.maxRequestSize = 5120*1024^2)

shinyServer(function(input, output, session) {
  
  ### Reactive variables
  # Reactive values are updated automatically when relevant user input changes
  rv <- reactiveValues(sc_obj=NULL
                       ,sc_obj_backup=NULL
                       ,cellCluster=c()
                       ,manualClusters=NULL
                       ,df.cluster.all=NULL
                       ,geneset = NULL
                       ,identset = NULL
                       ,drset = NULL
                       ,palette_cat=NULL
                       ,subset_cells=NULL
                       ,subset_cells_backup=NULL
                       ##########################################################
                       ,enrich_table=NULL
                       ##########################################################
                       
  )
  
  ### Basic page info
  
  output$t1 <- renderText({"scDAVIS: single-cell Data Analysis VISualization"})
  
  output$sessionInfo <- renderTable(colnames=F,{
    print(capture.output(sessionInfo()))   
  })
  
  shinyjs::hide(id = "upb")
  ### OBJECT CONTROL
  {
  ### LOADING AN UPLOADED (USER-PROVIDED) DATASET
  {
    updateSelectInput(session = session, inputId = "selectDataset", choices = as.list(list.files("PublicDatasets")),label = "Select a dataset")
    
    observeEvent(input$file, {
      if (is.null(input$file)) {
        return(NULL)
      }
      else {
        sc_obj <- input$file$datapath
        rv$filename <- basename(input$file$name)
        load_object(sc_obj = sc_obj, rv = rv, output = output, session = session, input = input)
        updateTabsetPanel(session
                          ,inputId = "sidetabs"
                          ,selected = "QCcheck")
        updateBox("module_qc_check-box_hgenes", action = "toggle")
        updateBox("module_qc_check-box_hreads", action = "toggle")
        updateBox("module_qc_check-box_info", action = "toggle")
      }
    })
    
    observeEvent(input$GET, {
      tryCatch({
        mytempdir <- tempdir(check=TRUE)
        mytempdir <- gsub("\\","/",mytempdir,fixed=T)
        mytempfile <- tempfile(tmpdir = mytempdir,fileext = ".rds")
        myfile <- file.path(mytempdir,basename(input$getURL))
        message(myfile)
        message(paste("getting:",myfile))
        shinyjs::show(id = "upb")
        GET(url = input$getURL,
            write_disk(mytempfile,overwrite = T),
            progress(session, id = "upb")
          )
        message(paste("after get"))
        file.copy(mytempfile, myfile)
        rv$filename <- myfile
        unlink(mytempfile)
        shinyjs::hide("upb")
        load_object(sc_obj = myfile, rv = rv, output = output, session = session, input = input)
        updateTabsetPanel(session
                          ,inputId = "sidetabs"
                          ,selected = "QCcheck")
        updateBox("module_qc_check-box_hgenes", action = "toggle")
        updateBox("module_qc_check-box_hreads", action = "toggle")
        updateBox("module_qc_check-box_info", action = "toggle")
        return(myfile)
      },
      error=function(cond) {
        message(paste("URL does not seem to exist:"))
        message("Here's the original error message:")
        message(cond)
        unlink(myfile)
        unlink(mytempfile)
        # Choose a return value in case of error
        return(NA)
      }
      )
    }
    )
  }
  
  ### LOADING A PRELOADED DATASET
  # Create seuratSC and fill basic info
  {
  public_data_server(id="module_public_data",rv)
    
  observeEvent(rv$publicfilename, {
    message(rv$publicfilename)
    sc_obj <- rv$publicfilename
    rv$filename <- basename(sc_obj)
    load_object(sc_obj = sc_obj, rv = rv, output = output, session = session, input = input)
    updateTabsetPanel(session
                      ,inputId = "sidetabs"
                      ,selected = "QCcheck")
    updateBox("module_qc_check-box_hgenes", action = "toggle")
    updateBox("module_qc_check-box_hreads", action = "toggle")
    updateBox("module_qc_check-box_info", action = "toggle")
  })
  }
    
  ### QC CHECKS
  {
    qc_check_server(id = "module_qc_check", rv )    
  }
  ###
    
  ### SWITCHING ASSAYS
    {
      observeEvent(input$select_assay, {
        if(!is.null(rv$sc_obj)){
          DefaultAssay(object = rv$sc_obj)<-input$select_assay
          rv$geneset <- rownames(rv$sc_obj)
          # updateSelectizeInput(session = session, inputId = "selectGene",choices = rv$geneset,server = T)
        }
      })
    }
    {
      observeEvent(rv$geneset, {
        if(!is.null(rv$sc_obj)){
          updateSelectizeInput(session = session, inputId = "selectGene",choices = rv$geneset,server = T)
        }
      })
    }
    
  
  ### DATA FILTER
  data_filter_server(id="module_data_filter",rv=rv)
  
  ### DOWNLOAD OBJECT
  data_download_server(id="module_data_download",rv=rv)
  }

  ### COLOR CONTROL
  {
  #### NEW Add Categorical Palette colors to plot ####
  observeEvent(input$selectCategorical,{
                 nlevels <- 8
                 if (!is.null(rv$sc_obj)) {
                   Idents(object = rv$sc_obj) <- getIdent(scObject = rv$sc_obj, identsColumns = input$selectCategorical)
                   ident_levels <- levels(Idents(object = rv$sc_obj))
                   nlevels <- length(ident_levels)
                 }
                 updatePalettePicker(session = session,
                                     inputId = "categoricalPalette",
                                     choices = getDiscretePals(n=nlevels,rev=input$revCatPal))
                
                                     
                 
                 updateSelectInput(session = session, inputId = "color_var",
                                   choices = c("",ident_levels)) # HE INTRODUCIDO MIXEDSORT AQUI
                 

                 rv$palette_cat<-update_palette(palette_updated=rv$palette_cat,input = input, regen=TRUE, rv=rv, ident_levels=ident_levels)
                 
                 
                 
  })
  
  observeEvent(input$categoricalPalette,
               {
                 colourpicker::updateColourInput(session,"color_cat",allowedCols= c("#FFFFFF",getDiscretePals(name=input$categoricalPalette,n=31)))
                 if (!is.null(rv$sc_obj)) {
                   ident_levels <- levels(Idents(object = rv$sc_obj))
                   nlevels <- length(ident_levels)
                   rv$palette_cat<-update_palette(palette_updated=rv$palette_cat,input = input, regen=TRUE, rv=rv, ident_levels=ident_levels)
                 }
               })
  
  observeEvent(input$revCatPal,
               {
                 nlevels <- 8
                 if (!is.null(rv$sc_obj)) {
                   ident_levels <- levels(Idents(object = rv$sc_obj))
                   nlevels <- length(ident_levels)
                 }
                 updatePalettePicker(session = session,
                                     inputId = "categoricalPalette",
                                     choices = getDiscretePals(n=nlevels,rev=input$revCatPal,default = input$categoricalPalette),
                 )
               }
  )

  #### NEW Add Numerical Palette colors to plot ####
  observeEvent(input$revNumPal,
               {
                 updatePalettePicker(session = session,
                                     inputId = "numericalPalette",
                                     choices = getContinuousPals(n=10,rev=input$revNumPal, default = input$numericalPalette),
                 )
               }
  )

  #### NEW Add Gene Palette colors to plot ####
  observeEvent(input$revGenPal,
               {
                 updatePalettePicker(session = session,
                                     inputId = "genePalette",
                                     choices = getContinuousPals(n=10,rev=input$revGenPal, default = input$genePalette),
                 )
               }
  )
  
  #### NEW Add Dense Palette colors to plot ####
  
  observeEvent(input$revDensPal,
               {
                 updatePalettePicker(session = session,
                                     inputId = "densityPalette",
                                     choices = getContinuousPals(n=10,rev=input$revDensPal, default = input$densityPalette),
                 )
               }
  )
  
  observe({
    if(!is.null(rv$sc_obj) & ((input$plot_tab %in% c("DimPlot","BarPlots") & input$radioSelectClust==0)| input$plot_tab=="Violin Plots" & input$radioSelectClust!=0)){
      enable("color_button")
    }else{
      disable("color_button")
    }
  })
  
  observeEvent(input$color_cat, {
    req(input$color_var,rv$sc_obj)
    ident_levels <- levels(Idents(object = rv$sc_obj))
    nlevels <- length(ident_levels)
    rv$palette_cat<-update_palette(palette_updated=rv$palette_cat,input = input, regen=FALSE, rv=rv, ident_levels=ident_levels)
  })
  }
  
  ### CONTRASTS
  {
  data_contrasts_server(id="module_data_contrasts",
                        rv=rv,
                        selectCategorical=reactive(input$selectCategorical),
                        selected_genes=reactive(input$selectGene)
                        )
    
    observeEvent(rv$numericaln, {
      newElement <- setdiff(rv$numericaln,rv$numerical)
      rv$numerical <- rv$numericaln
      updateSelectInput(session = session, inputId = "selectNumerical", choices = rv$numerical, selected = newElement)
      updateSelectInput(session = session, inputId = "radioSelectClust",selected = 1)
    })
    
    ### MARKERS TABLE
    {
      # Upload a diff table
      observeEvent(input$uploadTableFile, {
          if(!is.null(rv$sc_obj)){
            cdtf <- input$uploadTableFile$datapath
            cdt <- read.delim(cdtf,sep = "\t",header = T,quote = "\"",stringsAsFactors = F)
            if ("avg_log2FC" %in% colnames(cdt)) cdt <- cdt %>% dplyr::rename(avg_logFC=avg_log2FC)
            rv$df.cluster.all <- cdt
            rv$df.cluster <- cdt[cdt$p_val_adj < 0.05,]
            rv$df.contrast <- input$uploadTableFile$name
          }
        })
      
      # Fill table for mks cluster
      observeEvent(rv$df.cluster.all, {
        if (exists("dtproxy")) selectRows(proxy = dtproxy, selected = NULL)
        output$df.mks <- renderDataTable({
          if(!is.null(rv$sc_obj)){
            
            shiny::validate(
              need(rv$df.cluster.all, "Please run a contrast in Find Markers tab to create a table"))
            
            DT::datatable(
              rv$df.cluster.all,
              rownames = F,
              filter="top",
              escape = F,
              selection = 'multiple',
              caption = paste("Table: Contrast",rv$df.contrast),
              options = list(pageLength = 10
                             ,lengthMenu = c(10, 20, 50)
                             ,lengthChange=FALSE
                             ,stateSave = TRUE
                             ,scrollX = TRUE
              )) %>%
              formatRound(columns = c(3:ncol(rv$df.cluster.all)), digits = 2) %>%
              formatStyle(columns = c(1:ncol(rv$df.cluster.all)), 'text-align' = 'center')
          }
        })
        dtproxy <- dataTableProxy("df.mks")
      })
      
      # Clear table selections when modifying the state
      observeEvent(input$clearSelect, {
       
        #################################################################################################
        # selectRows(proxy = dtproxy, selected = NULL)
        
        if (exists("dtproxy")) {
          DT::selectRows(proxy = dtproxy, selected = NULL)
        }
        
        rv$df.cluster.all <- NULL
        rv$df.cluster <- NULL
        rv$df.contrast <- NULL
        
        output$df.mks <- DT::renderDataTable({
          data.frame()
        })
        #################################################################################################
        
      })
      
      # Update plot with genes selected from table 
      observeEvent(input$df.mks_rows_selected, {
        updateSelectizeInput(session = session
                             , inputId = "selectGene"
                             , choices = rv$geneset
                             , selected = rv$df.cluster.all[input$df.mks_rows_selected,grep("^gene$",colnames(rv$df.cluster.all))]
                             , server = T
        )
      })
      
      # Button for download mks in cluster tab
      output$downloadDataMks <- downloadHandler(
        filename = function() {
          paste(paste(sub(".rds","",rv$filename),gsub("[.]", "_", input$selectCategorical), 
                      rv$df.contrast,"mks", sep="_"), ".csv", sep = "")
        },
        content = function(file) {
          write.table(rv$df.cluster.all, file, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
      )

    }
  }
  
  ### MANUAL CLUSTERING
  {
    
    output$brushCluster <- renderText({
      length(rv$cellCluster)
    })
    
    observeEvent(event_data("plotly_selected"), {
      d <- event_data("plotly_selected")
      rv$cellCluster <- d$key
    }
    )
    
    shinyInput <- function(FUN, len, id, ...) {
      inputs <- character(len)
      for (i in seq_len(len)) {
        inputs[i] <- as.character(FUN(paste0(id, i), ...))
      }
      inputs
    }
    
    # Gate Manually Selected Cells
    observeEvent(input$gateCells, {
      rv$subset_cells_backup[[length(rv$subset_cells_backup)+1]] <- rv$subset_cells
      rv$subset_cells <- rv$cellCluster
    })
    
    observeEvent(input$undoGateCells, {
      if (length(rv$subset_cells_backup)) {
        rv$subset_cells <- rv$subset_cells_backup[[length(rv$subset_cells_backup)]]
        rv$subset_cells_backup[[length(rv$subset_cells_backup)]] <- NULL
      }
    })
    
    # Save Manually Selected Cell Cluster
    observeEvent(input$saveCellCluster, {
      if (is.null(rv$manualClusters)) {rv$manualClusters <- list()}
      if (input$nameCellCluster == "") {
        output$cwarnings <- renderText({c("Cluster name is not set")})
      } else if (input$nameCellCluster %in% names(rv$manualClusters)) {
        output$cwarnings <- renderText({c("Cluster name is repeated")})
      } else {
        output$cwarnings <- renderText("")
        rv$manualClusters[[input$nameCellCluster]] <- unlist(rv$cellCluster)
        output$cellCluster <- DT::renderDataTable(
          data.frame(
            Cluster=names(rv$manualClusters),
            Cells=unlist(lapply(rv$manualClusters,length))
            , Remove=shinyInput(actionButton, length(rv$manualClusters), 'button_', label = "", icon = icon("times-circle"), onclick = 'Shiny.onInputChange(\"remove_button\", [this.id, Math.random()])'))
          , rownames = FALSE
          , server = FALSE
          , escape = FALSE
          , editable = TRUE
          , selection = 'multiple'
          , options = list(info = FALSE #Hides the info al the bottom ("Showing 1 to N of N entries")
                           #dom = 't' #This option would hide all the extra elements (shows only the table). However we want to show
                           #the page number
                           ,sDom  = '<"top">lrt<"bottom">ip' #Hides search option
                           ,pageLength = 10 
                           ,lengthChange=FALSE #Introducing these parameters, the number of entries per page cannot be modifiable and
                           #this option will not appear in the app (better visualization)
                           ,columnDefs = list(list(className = 'dt-center', targets="_all")))
        )
        rv$cellCluster <- c()
      }
    })
    
    # Remove selected row in Cluster Table (Select Cells)
    observeEvent(input$remove_button, {
      selectedRow <- as.numeric(strsplit(input$remove_button, "_")[[1]][2])
      clusterName <- names(rv$manualClusters)[selectedRow]
      rv$manualClusters[clusterName] <- NULL
      if (length(rv$manualClusters) == 0) {
        rv$manualClusters <- NULL
        output$cellCluster <- DT::renderDataTable(
          expr = data.frame()
          , rownames = FALSE
          , server = FALSE
          , options = list(dom = 't' #This option would hide all the extra elements (shows only the table). However we want to show
                           #the page number
                           ,pageLength = 0 
                           ,lengthChange=FALSE
          ))
      } else {
        output$cellCluster <- NULL
        output$cellCluster <- DT::renderDataTable(
          expr = data.frame(
            Cluster=names(rv$manualClusters)
            , Cells=unlist(lapply(rv$manualClusters,length))
            , Remove=shinyInput(actionButton, length(rv$manualClusters), 'button_', label = "", icon = icon("times-circle"), onclick = 'Shiny.onInputChange(\"remove_button\", [this.id, Math.random()])')
          )
          , rownames = FALSE
          , server = FALSE
          , escape = FALSE
          , editable = TRUE
          , selection = 'single'
          , options = list(info = FALSE #Hides the info al the bottom ("Showing 1 to N of N entries")
                           #dom = 't' #This option would hide all the extra elements (shows only the table). However we want to show
                           #the page number
                           ,sDom  = '<"top">lrt<"bottom">ip' #Hides search option
                           ,pageLength = 10 
                           ,lengthChange=FALSE #Introducing these parameters, the number of entries per page cannot be modifiable and
                           #this option will not appear in the app (better visualization)
                           ,columnDefs = list(list(className = 'dt-center', targets="_all")))
        )  
      }
    })
    
    # Save new Clustering Set
    observeEvent(input$saveNewClustering, {
      shiny::validate(
        need(rv$sc_obj, "Please load a seurat object"))
      if (input$nameNewClustering == ""){
        output$cswarnings <- renderText({c("Please introduce a name for the new Clustering Set")})
      } else {
        if(input$nameNewClustering %in% rv$categorical){
          output$cswarnings <- renderText({c("There is already a Clustering Set with this name. Do you want to overwrite the existing cluster?")})
          updateRadioButtons(session, "invvar", label = NULL
                             , choices = list("?" = 0, "Continue" = 1) #Invisible control variable. Will be 1 in case of overwritting Clustering Set Name
                             , selected = 1)
          shinyjs::toggle("clusterControl")
        }
        if (input$invvar == 1 || !(input$nameNewClustering %in% rv$categorical)){
          if(input$radioNewClusterControl == 1){
            output$cswarnings <- NULL
            updateRadioButtons(session, "invvar", label = NULL
                               , choices = list("?" = 0, "Continue" = 1)
                               , selected = 0)
            updateRadioButtons(session, "radioNewClusterControl", label = NULL
                               , choices = list("Yes" = 0, "No" = 1)
                               , selected = 0)
            shinyjs::hide("clusterControl")
          }
          if(input$radioNewClusterControl == 0){
            output$cswarnings <- NULL
            withProgress(message = "Saving Manual Clustering",{
              mix <- unique(sort(unlist(rv$manualClusters)[duplicated(unlist(rv$manualClusters))])) # get repeated cells
              if(input$checkboxKeepIdentCluster){
                newIdent <- data.frame(newIdent = as.character(rv$sc_obj@meta.data[, input$selectCategorical]),stringsAsFactors = F)
                rownames(newIdent) <- rownames(rv$sc_obj@meta.data)
              } else {
                newIdent <- data.frame(newIdent = rep("Mix",dim(rv$sc_obj@meta.data)[1]),stringsAsFactors = F)
                rownames(newIdent) <- rownames(rv$sc_obj@meta.data)
              }
              for (c in names(rv$manualClusters)) {
                cellNames <- rv$manualClusters[[c]][!rv$manualClusters[[c]] %in% mix]
                newIdent[cellNames,"newIdent"] <- c
              }
              rv$sc_obj@meta.data[,input$nameNewClustering] <- newIdent$newIdent
              rv$categorical <- c(rv$categorical,input$nameNewClustering)
              
              updateCat <- input$nameNewClustering
              
              updateSelectInput(session = session, inputId = "selectCategorical", choices = rv$categorical, selected = updateCat)
              
              rv$manualClusters <- NULL
              updateTextInput(session = session, inputId = "nameCellCluster",value = "")
              updateTextInput(session = session, inputId = "nameNewClustering",value = "")
              rv$cellCluster <- c()
              output$cellCluster <- DT::renderDataTable( 
                expr = data.frame()
                , rownames = FALSE
                , server = FALSE
                , escape = FALSE
                , editable = FALSE
                , selection = 'single'
                , options = list(dom = 't' #This option would hide all the extra elements (shows only the table). However we want to show
                                 #the page number
                                 ,pageLength = 0 
                                 ,lengthChange=FALSE
                )
              )
              updateRadioButtons(session, "invvar", label = NULL
                                 , choices = list("?" = 0, "Continue" = 1)
                                 , selected = 0)
              shinyjs::hide("clusterControl")
              output$cswarnings <- renderText({c("Cluster correctly created")})
            })
          }
        }
      }
    })
    
    # Save new Ident
    observeEvent(input$saveIdent, {
      rv$gene <- input$selectGene
      mix <- intersect(rv$cell1, rv$cell2)
      cell1 <- rv$cell1[!rv$cell1 %in% mix]
      cell2 <- rv$cell2[!rv$cell2 %in% mix]
      mix <- c(mix, colnames(x = rv$sc_obj)[!colnames(x = rv$sc_obj) %in% c(cell1, cell2)])
      if(input$checkboxKeepIdent){
        Idents(object = rv$sc_obj)<-rv$sc_obj@meta.data[, input$nameClusteringIdent]
        mix <- WhichCells(rv$sc_obj, ident = "Mix") 
        mix <- mix[!mix %in% c(cell1, cell2)]
      }
      rv$sc_obj@meta.data[cell1, input$nameClusteringIdent] = input$nameClusteringCell1
      rv$sc_obj@meta.data[cell2, input$nameClusteringIdent] = input$nameClusteringCell2
      rv$sc_obj@meta.data[mix, input$nameClusteringIdent] = "Mix"
      rv$identset <- colnames(rv$sc_obj@meta.data)
      updateSelectInput(session = session, inputId = "selectCategorical", choices = rv$identset, selected = rv$ident)
    })
  }
  
  #### PLOT TABS
  {
    ### DIMPLOT TAB 
    {
      plot_reductions_server(id="module_plot_reductions",
                             rv=rv,
                             selectRepresentation=reactive(input$selectRepresentation),
                             radioSelectClust=reactive(input$radioSelectClust),
                             categorical_list =reactive(rv$categorical),
                             selectedCategory=reactive(input$selectCategorical),
                             selectedNumeric=reactive(input$selectNumerical),
                             selectedGene=reactive(input$selectGene),
                             selectDensity=reactive(input$selectDensity),
                             categoricalPalette=reactive(input$categoricalPalette),
                             revCatPal=reactive(input$revCatPal),
                             numericalPalette=reactive(input$numericalPalette),
                             revNumPal=reactive(input$revNumPal),
                             genePalette=reactive(input$genePalette),
                             revGenPal=reactive(input$revGenPal),
                             densityPalette=reactive(input$densityPalette),
                             revDensPal=reactive(input$revDensPal),
                             palette_cat=reactive(rv$palette_cat)
      )
    }
    
    ### VLN PLOT TAB
    {
      plot_vln_server(id="module_plot_vln"
                      , sc_obj=reactive(rv$sc_obj)
                      , selectedCategory=reactive(input$selectCategorical)
                      , radio_select=reactive(input$radioSelectClust)
                      , selectedGene =reactive(input$selectGene)
                      , selectedNumeric=reactive(input$selectNumerical)
                      , categoricalPalette=reactive(rv$palette_cat)
      )
    }

    ### HEATMAP TAB
    {
      plot_heatmap_server(id="module_plot_heatmap"
                          , rv = rv
                          , sc_obj=reactive(rv$sc_obj)
                          , category=reactive(input$selectCategorical)
                          , palette_cat=reactive(input$categoricalPalette)
                          , palettec_rev=reactive(input$revCatPal)
                          , palette_gene=reactive(input$genePalette)
                          , paletteg_rev=reactive(input$revGenPal)
                          , selected_genes=reactive(input$selectGene)
                          , cat_vars = reactive(rv$categorical)
                          , num_vars = reactive(rv$numerical)
      )
    }
    
    ### BARPLOT TAB
    {
      plot_bar_server("module_plot_bar",sc_obj=reactive(rv$sc_obj),category=reactive(input$selectCategorical), categorical_list =reactive(rv$categorical), radio_select=reactive(input$radioSelectClust), palette_bplot=reactive(rv$palette_cat))
    }
    
    ### SCATTERPLOT TAB
    {
      plot_scatter_server(id="module_plot_scatter",
                          rv=rv,
                          numerical=reactive(input$selectNumerical),
                          gene_select=reactive(input$selectGene),
                          radio_select=reactive(input$radioSelectClust),
                          palette_cat=reactive(input$categoricalPalette),
                          palettec_rev=reactive(input$revCatPal),
                          palette_numerical=reactive(input$numericalPalette),
                          paletten_rev=reactive(input$revNumPal),
                          palette_gene=reactive(input$genePalette),
                          paletteg_rev=reactive(input$revGenPal),
                          numerical_list=reactive(rv$numerical),
                          categorical_list=reactive(rv$categorical),
                          gene_list=reactive(rv$geneset),
                          # DJC
                          palette_dens=reactive(input$densityPalette),
                          paletted_rev = reactive(input$revDensPal),
                          slider_dens=reactive(input$selectDensity),
                          #######
                          ## DJC
                          cellCluster=reactive(rv$cellCluster),
                          assay=reactive(input$select_assay)
                          ## DJC
                          #######
      )
    }
    
    ### DOTPLOT TAB
    {
      plot_dot_server("module_plot_dot",sc_obj=reactive(rv$sc_obj),category=reactive(input$selectCategorical), 
                      gene_select =reactive(input$selectGene), radio_select=reactive(input$radioSelectClust), 
                      palette_gene=reactive(input$genePalette),palette_rev=reactive(input$revGenPal),
                      df_table=reactive(rv$df.cluster.all))
    }
    
    
    ##################################################################################################################
    
    ### ENRICHR TAB
    
    {
      enrichr_res <- enrichR_server(
        id = "module_enrichr",
        gene_select = reactive(input$selectGene),
        dif_table = reactive(rv$df.cluster.all),
        radio_select = reactive(input$radioSelectClust)
      )
      
      observeEvent(enrichr_res$enrich_table(), {
        rv$enrich_table <- enrichr_res$enrich_table()
      })
    }
    
    ##################################################################################################################
    
  }
  
})
  
