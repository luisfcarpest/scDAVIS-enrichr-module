
# ======================================================================================================
# FILE MODULE_ENRICHR.R 
# 
# Módulo de enriquecimiento funcional con EnrichR integrado en scDAVIS
# ======================================================================================================

enrichR_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
  
  # Ajuste visual para que los nombres largos de las bases de datos no choquen con la flecha
    
  tags$style(HTML(sprintf("
    #%s .selectize-input {
    padding-right: 32px !important;
  }

  #%s .selectize-input .item {
    display: block;
    max-width: calc(100%% - 18px);
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }
", ns("database_container"), ns("database_container")))),
    
    fluidRow(
      div(
        style = "position:absolute;right:3em;",
        
        dropdownButton(
          size = "sm",
          icon = icon(name = "sliders-h"),
          up = FALSE,
          right = TRUE,
          circle = FALSE,
          width = 700,
          
          # ==============================================
          # Download Options
          # ==============================================
          column(
            width = 4,
            h4(strong("Download Options")),
            numericInput(ns("oWidth"), "Output Width", min = 4, max = 20, value = 10),
            numericInput(ns("aRatio"), "Aspect Ratio", min = 0.4, max = 2, value = 0.8, step = 0.1),
            numericInput(ns("oRes"), "Output Resolution", min = 72, max = 600, value = 150),
            
            h5(strong("Download Plot")),
            
            fluidRow(
              column(
                6,
                downloadButton(
                  ns("download_enrichplot"),
                  "Download plot",
                  class = "btn-sm"
                )
              ),
              column(
                6,
                selectInput(
                  ns("format"),
                  label = NULL,
                  choices = c("jpeg", "pdf"),
                  selected = "jpeg"
                )
              )
            ),
            
            br(),
            
            h5(strong("Download Table")),
            
            downloadButton(
              ns("download_enrichtable"),
              "Download table",
              class = "btn-sm"
            )
          ),
          
          # ==============================================
          # Plot Options
          # ==============================================
          column(
            width = 4,
            sliderInput(ns("n_terms"), "Number of enriched terms", min = 5, max = 30, value = 10, step = 1),
            sliderInput(ns("xLabelsSize"), "X Font Size", min = 6, max = 20, value = 10, step = 1),
            sliderInput(ns("yLabelsSize"), "Y Font Size", min = 6, max = 20, value = 10, step = 1),
            sliderInput(ns("yLabelsWrap"), "Wrap term labels at", min = 20, max = 80, value = 40, step = 5)
          ),
          
          
          # ==============================================
          # Enrichment Options
          # ==============================================
          column(
            width = 4, 
            radioButtons(
              ns("gene_source"),
              label = h4("Gene source"),
              choices = list(
                "My geneset" = "selected",
                "Upload gene list (.txt)" = "upload",
                "Top genes from diff table" = "diff"
              ),
              selected = "selected"
            ),
            
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'upload'", ns("gene_source")),
              fileInput(
                ns("upload_gene_list"),
                "Upload gene list (.txt)",
                accept = c(".txt")
              )
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == 'diff'", ns("gene_source")),
              
              sliderInput(
                ns("n_diff_genes"),
                "Number of diff genes",
                min = 5,
                max = 100,
                value = 20,
                step = 5
              )
            ),
            
            div(
              id = ns("database_container"),
              selectizeInput(
                ns("database"),
                "Database",
                choices = NULL,
                selected = "GO_Biological_Process_2023",
                multiple = TRUE,
                width = "100%",
                options = list(
                  placeholder = "Search Enrichr databases"
                )
              )
            ),
            
            br(),
            
            actionButton(
              ns("run_enrich"),
              "Run enrichment",
              class = "btn-sm"
            ),
            
            br(),
            
            
            checkboxInput(ns("show_table"), "Show enrichment table", value = FALSE)
          )
        )
      )
    ),
    
    br(),
    br(),
    
    fluidRow(
      column(
        width = 12, 
        plotOutput(ns("plot_enrich"), width = "95%", height = "auto") %>% withSpinner())
      ),
    
    conditionalPanel(
      condition = sprintf("input['%s'] == true", ns("show_table")),
      fluidRow(
        column(
          width = 12,
          DT::dataTableOutput(ns("table_enrich"))
        )
      )
    )
    
  )
}


enrichR_server <- function(id, gene_select, dif_table, radio_select) {
  
  moduleServer(id, function(input, output, session) {
    
    
    message("Loading module_enrichR")
    
    
    # ===============================================================
    # 1. Bases de datos disponibles en Enrichr
    # ===============================================================
    available_dbs <- reactive({
      shiny::validate(
        need(requireNamespace("enrichR", quietly = TRUE), "Package 'enrichR' is not installed")
      )
      
      dbs <- enrichR::listEnrichrDbs()
      
      shiny::validate(
        need(!is.null(dbs), "Could not retrieve Enrichr databases"),
        need(nrow(dbs) > 0, "No Enrichr databases available"),
        need("libraryName" %in% colnames(dbs), "Enrichr database list does not contain 'libraryName'")
      )
      
      dbs
    })
    
    
    # ===============================================================
    # 2. Actualizar selector de bases de datos
    # ===============================================================
    observe({
      dbs <- available_dbs()
      db_names <- dbs$libraryName
      
      # Bases recomendadas por defecto
      preferred <- c(
        "GO_Biological_Process_2023",
        "GO_Molecular_Function_2023",
        "GO_Cellular_Component_2023",
        "KEGG_2021_Human",
        "Reactome_2022",
        "WikiPathways_2024_Human"
      )
      
      preferred <- preferred[preferred %in% db_names]
      ordered_choices <- c(preferred, sort(setdiff(db_names, preferred)))
      
      shiny::validate(
        need(length(ordered_choices) > 0, "No Enrichr databases available to populate selector")
      )
      
      updateSelectizeInput(
        session = session,
        inputId = "database",
        choices = ordered_choices,
        selected = if ("GO_Biological_Process_2023" %in% ordered_choices) {
          "GO_Biological_Process_2023"
        } else {
          ordered_choices[1]
        },
        server = TRUE
      )
    })
    
    
    # ===============================================================
    # 3. Función auxiliar para ejecutar enrichR
    # ===============================================================
    run_enrichr <- function(genes, dbs) {
      genes <- unique(genes)
      
      if (length(genes) < 2) {
        return(NULL)
      }
      
      res_list <- tryCatch(
        enrichR::enrichr(genes, databases = dbs),
        error = function(e) NULL
      )
      res_list <- res_list %>% discard(~ is.null(.x)) %>% discard(~ dim(.x)[1] == 0)
      
      if (is.null(res_list) | length(res_list) == 0) {
        return(NULL)
      }
      
      res <- dplyr::bind_rows(res_list, .id = "database")
      
      if (nrow(res) == 0) {
        return(NULL)
      }
      
      res
    }
    
    
    # ===============================================================
    # 4. Selección de genes para enriquecimiento normal
    # ===============================================================
    selected_gene_sets <- reactive({
      shiny::validate(
        need(radio_select() == 2, "Please select the Variable Type Gene for EnrichR representation")
      )
      
      if (input$gene_source == "selected") {
        shiny::validate(
          need(gene_select(), "Please select one or more genes in scDAVIS")
        )
        
        genes <- unique(gene_select())
        
        shiny::validate(
          need(length(genes) > 0, "No genes available from selected geneset")
        )
        
        return(list("My geneset" = genes))
      }
      
      if (input$gene_source == "upload") {
        shiny::validate(
          need(input$upload_gene_list$datapath, "Please upload a .txt gene list")
        )
        
        genes <- read.delim(
          input$upload_gene_list$datapath,
          header = FALSE,
          stringsAsFactors = FALSE
        )[, 1]
        
        genes <- unique(trimws(genes))
        genes <- genes[genes != ""]
        
        shiny::validate(
          need(length(genes) > 0, "No genes found in uploaded file")
        )
        
        return(list("Uploaded gene list" = genes))
      }
      
      if (input$gene_source == "diff") {
        shiny::validate(
          need(dif_table(), "Please generate or upload a differential table first")
        )
        
        df <- dif_table()
        
        shiny::validate(
          need("gene" %in% colnames(df), "The differential table must contain a 'gene' column"),
          need("avg_logFC" %in% colnames(df), "The differential table must contain an 'avg_logFC' column")
        )
        
        if ("cluster" %in% colnames(df)) {
          split_df <- df %>%
            dplyr::group_by(cluster) %>%
            dplyr::group_split()
          
          gene_sets <- purrr::map(split_df, function(cluster_df) {
            cluster_df %>%
              dplyr::arrange(desc(avg_logFC)) %>%
              dplyr::slice_head(n = input$n_diff_genes) %>%
              dplyr::pull(gene) %>%
              unique()
          })
          
          names(gene_sets) <- df %>%
            dplyr::group_by(cluster) %>%
            dplyr::group_keys() %>%
            dplyr::pull(cluster)
          
          return(gene_sets)
        }
        
        genes <- df %>%
          dplyr::arrange(desc(avg_logFC)) %>%
          dplyr::slice_head(n = input$n_diff_genes) %>%
          dplyr::pull(gene) %>%
          unique()
        
        return(list("Top diff genes" = genes))
      }
    })
    
    
    # ===============================================================
    # 5. Tabla de enriquecimiento normal
    # ===============================================================
    enrich_table <- eventReactive(input$run_enrich, {
      shiny::validate(
        need(requireNamespace("enrichR", quietly = TRUE), "Package 'enrichR' is not installed")
      )
      
      last_gene_source(input$gene_source)
      
      dbs <- input$database
      gene_sets <- selected_gene_sets()
      
      shiny::validate(
        need(!is.null(dbs) && length(dbs) > 0, "Please select at least one Enrichr database"),
        need(length(gene_sets) > 0, "No gene sets available for enrichment")
      )
      
      res_list <- purrr::imap(gene_sets, function(genes, gene_set_name) {
        res <- run_enrichr(genes, dbs)
        
        if (is.null(res)) {
          return(NULL)
        }
        
        res <- res %>%
          dplyr::mutate(
            gene_set = gene_set_name,
            minus_log10_padj = -log10(Adjusted.P.value + 1e-300),
            overlap_count = as.numeric(sub("/.*", "", Overlap))
          )
        return(res) 
      
      }) %>%
        discard(~ is.null(.x)) %>%
        discard(~ dim(.x)[1] == 0)
      
      
      res <- dplyr::bind_rows(res_list)
      
      shiny::validate(
        need(nrow(res) > 0, "No enrichment results found")
      )
      
      res %>%
        dplyr::select(gene_set, database, dplyr::everything()) %>%
        dplyr::arrange(gene_set, database, Adjusted.P.value)
    })
    
    
    
    # ===============================================================
    # 6. Datos preparados para el plot normal
    # ===============================================================
    initial_message <- reactive({
      if (input$gene_source == "selected") {
        return("Select one or more genes in the right panel, choose one or more databases, and click Run enrichment.")
      }
      
      if (input$gene_source == "upload") {
        return("Upload a .txt file with one gene per line, choose one or more databases, and click Run enrichment.")
      }
      
      if (input$gene_source == "diff") {
        return("Run mks first in the Contrasts tab, choose one or more databases, and click Run enrichment.")
      }
    })
    
    enrich_plot_data <- reactive({
      
      df <- selected_enrich_block()
      
      df <- df %>%
        dplyr::slice_head(n = input$n_terms) %>%
        dplyr::mutate(
          GeneCount = as.numeric(sub("/.*", "", Overlap)),
          ShortTerm = stringr::str_wrap(Term, width = input$yLabelsWrap)
        )
      
      df$ShortTerm <- factor(df$ShortTerm, levels = rev(unique(df$ShortTerm)))
      
      df
    })
    
    
    # ===============================================================
    # 7. Construcción del plot normal
    # ===============================================================
    enrich_plot <- reactive({
      df <- enrich_plot_data()
      
      ggplot(df, aes(x = GeneCount, y = ShortTerm, fill = P.value)) +
        geom_col(width = 0.82, color = "white", linewidth = 0.2) +
        scale_fill_gradientn(
          colours = c("#F46D43", "#D948A1", "#4C6EDB"),
          name = expression(italic("P value"))
        ) +
        labs(
          title = "Enrichment Analysis by EnrichR",
          subtitle = paste(unique(df$gene_set), "-", unique(df$database)),
          x = "Gene count",
          y = "Enriched terms"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.background = element_rect(fill = "#FFFFFF", color = NA),
          panel.background = element_rect(fill = "#FFFFFF", color = NA),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "#DCE8F2", linewidth = 0.4),
          axis.title.x = element_text(color = "#223447", face = "bold", size = input$xLabelsSize),
          axis.title.y = element_text(color = "#223447", face = "bold", size = input$yLabelsSize),
          axis.text.x = element_text(color = "#35506A", size = input$xLabelsSize),
          axis.text.y = element_text(color = "#35506A", size = input$yLabelsSize),
          plot.title = element_text(face = "bold", size = 16, color = "#1E3A5F"),
          plot.subtitle = element_text(color = "#5A748C", size = 11)
        )
    })
    
    
    # ===============================================================
    # 8. Render plot
    # ===============================================================
    last_gene_source <- reactiveVal(NULL)
    
    observeEvent(input$gene_source, {
      last_gene_source(NULL)
    })
    
    output$plot_enrich <- renderPlot({
      
      if (is.null(last_gene_source()) || last_gene_source() != input$gene_source) {
        par(mar = c(0, 0, 0, 0))
        plot.new()
        text(
          0.5, 0.5,
          initial_message(),
          cex = 1.1,
          col = "gray30"
        )
        return(invisible(NULL))
      }
      
      print(enrich_plot())
      
    }, height = function() {
      h <- 120 + input$n_terms * 35
      if (h < 400) h <- 400
      h
    })
    
    
    # ===============================================================
    # 9. Render table
    # ===============================================================
    output$table_enrich <- DT::renderDataTable({
      df <- enrich_table()
      
      if (input$gene_source %in% c("selected", "upload")) {
        df_show <- df %>%
          dplyr::select(-gene_set)
      } else {
        df_show <- df %>%
          dplyr::rename (Cluster = gene_set)      }
      
      DT::datatable(
        df_show,
        rownames = FALSE,
        filter = "top",
        selection = "single",
        options = list(
          pageLength = 10,
          lengthChange = FALSE,
          scrollX = TRUE
        )
      ) %>%
        DT::formatRound(
          columns = intersect(c("P.value", "Adjusted.P.value", "Combined.Score"), colnames(df_show)),
          digits = 4
        )
    })
    
    
    selected_enrich_index <- reactiveVal(1)
    
    observeEvent(enrich_table(), {
      selected_enrich_index(1)
    })
    
    observeEvent(input$table_enrich_rows_selected, {
      if (!is.null(input$table_enrich_rows_selected) &&
          length(input$table_enrich_rows_selected) > 0) {
        selected_enrich_index(input$table_enrich_rows_selected[1])
      }
    })
  
    
    selected_enrich_block <- reactive({
      df <- enrich_table()
      
      shiny::validate(
        need(!is.null(df), "Please run enrichment first"),
        need(nrow(df) > 0, "No enrichment results available")
      )
      
      selected_index <- selected_enrich_index()
      
      if (is.null(selected_index) || selected_index > nrow(df)) {
        selected_index <- 1
      }
      
      selected_row <- df[selected_index, , drop = FALSE]
      
      df %>%
        dplyr::filter(
          gene_set == selected_row$gene_set[1],
          database == selected_row$database[1]
        )
    })
    
    
    # ===============================================================
    # 10. Download plot
    # ===============================================================
    output$download_enrichplot <- downloadHandler(
      filename = function() {
        ext <- ifelse(input$format == "pdf", ".pdf", ".jpeg")
        paste0("EnrichR_", input$database, "_", input$gene_source, ext)
      },
      
      content = function(file) {
        if (input$format == "pdf") {
          pdf(file, width = input$oWidth, height = input$oWidth * input$aRatio)
          print(enrich_plot())
          dev.off()
        } else {
          jpeg(
            file,
            res = input$oRes,
            units = "in",
            width = input$oWidth,
            height = input$oWidth * input$aRatio
          )
          print(enrich_plot())
          dev.off()
        }
      },
      
      contentType = "application/pdf"
    )
    
    
    # ===============================================================
    # 11. Download enrichment table
    # ===============================================================
    output$download_enrichtable <- downloadHandler(
      
      filename = function() {
        paste0("EnrichR_table_", input$gene_source, ".csv")
      },
      
      content = function(file) {
        df <- enrich_table()
        
        write.csv(
          df,
          file,
          row.names = FALSE
        )
      }
    )
    
    
    # ===============================================================
    # 12. Devolver reactivos al server principal
    # ===============================================================
    return(list(
      enrich_table = enrich_table,
      enrich_plot = enrich_plot,
      selected_gene_sets = selected_gene_sets
    ))
    
  })
}

