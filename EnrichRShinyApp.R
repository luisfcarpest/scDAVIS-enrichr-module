# =========================
# LIBRERÍAS
# =========================

# shinythemes: para aplicar temas visuales de shiny
# DT: muestra tablas interactivas
# ggplot2: genera gráficos
# dplyr: manipulación de datos
# enrichR: consulta enriquecimiento funcional en EnrichR

library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(dplyr)
library(enrichR)


# =========================
# INTERFAZ DE USUARIO (UI)
# =========================

ui <- fluidPage(
  
  # Tema visual general
  theme = shinytheme("flatly"),
  
  # CSS simplificado
  tags$head(
    tags$style(HTML("
      body {
        background-color: #F5F7FA;
        font-family: Arial, sans-serif;
      }

      h1, h2, h3, h4 {
        color: #2C3E50;
      }

      .well, .tab-content {
        background-color: white;
        border-radius: 10px;
        padding: 15px;
      }

      .nav-tabs > li > a {
        background-color: #EAEFF5;
        color: #2C3E50;
        border-radius: 6px 6px 0 0;
      }

      .nav-tabs > li.active > a {
        background-color: #2C7FB8 !important;
        color: white !important;
      }

      .btn-primary {
        background-color: #2C7FB8;
        border: none;
        border-radius: 6px;
      }

      table.dataTable thead {
        background-color: #EAEFF5;
      }
    "))
  ),
  
  # Título principal de la app
  titlePanel("Single-Cell Functional Enrichment Explorer"),
  
  sidebarLayout(
    
# =========================
# PANEL LATERAL
# =========================

    sidebarPanel(
      width = 4,
      
      h4("Gene Input"),
      
      # Selección de origen de genes
      radioButtons(
        "gene_source",
        "Select gene source",
        choices = c(
          "Use EnrichR example genes" = "example",
          "Upload or paste custom set genes" = "custom"
        ),
        selected = "example"
      ),
      
      # Este panel solo aparece si el usuario elige genes personalizados
      conditionalPanel(
        condition = "input.gene_source == 'custom'",
        
        # Subida de archivo de texto
        fileInput(
          "gene_file",
          "Upload a .txt file with genes",
          accept = ".txt"
        ),
        
        # Entrada manual de genes
        textAreaInput(
          "manual_genes",
          "Or paste gene symbols (one per line or comma-separated)",
          rows = 8,
          placeholder = "TP53\nEGFR\nBRCA1\nBRCA2\nMYC"
        )
      ),
      
      br(),
      h4("EnrichR Databases"),
      
      # Selector dinámico de bases de datos
      uiOutput("db_selector"),
      
      # Número de términos a mostrar en la visualización
      numericInput(
        "top_n",
        "Number of top terms to display",
        value = 20,
        min = 5,
        max = 50
      ),
      
      # Botón para ejecutar el análisis
      actionButton("run_enrich", "Run Enrichment")
    ),
    
    
# =========================
# PANEL PRINCIPAL
# =========================

    mainPanel(
      width = 8,
      
      tabsetPanel(
        
        # Pestaña para ver genes cargados
        tabPanel(
          "Loaded Genes",
          br(),
          verbatimTextOutput("gene_source_info"),
          br(),
          DTOutput("gene_table")
        ),
        
        # Pestaña para ver resultados del enriquecimiento
        tabPanel(
          "Results",
          br(),
          uiOutput("result_tabs")
        ),
        
        # Pestaña para visualizar resultados en gráfico
        tabPanel(
          "Visualization",
          br(),
          plotOutput("enrich_plot", height = "650px")
        )
      )
    )
  )
)


# =========================
# SERVIDOR
# =========================

server <- function(input, output, session) {
  
  # -------------------------
  # 1. Cargar bases de datos disponibles de EnrichR
  # -------------------------
  enrichr_dbs <- reactive({
    tryCatch(
      listEnrichrDbs(),
      error = function(e) NULL
    )
  })
  
  
  # -------------------------
  # 2. Mostrar selector de bases de datos
  # -------------------------
  output$db_selector <- renderUI({
    dbs <- enrichr_dbs()
    
    if (is.null(dbs) || nrow(dbs) == 0) {
      helpText("Could not load EnrichR databases.")
    } else {
      selectizeInput(
        "selected_dbs",
        "Select one or more databases",
        choices = dbs$libraryName,
        selected = c(
          "GO_Molecular_Function_2023",
          "GO_Cellular_Component_2023",
          "GO_Biological_Process_2023"
        ),
        multiple = TRUE,
        options = list(
          placeholder = "Choose EnrichR libraries",
          maxOptions = 500
        )
      )
    }
  })
  
  
  # -------------------------
  # 3. Cargar y limpiar genes
  # -------------------------
  genes_reactive <- reactive({
    
    # Si el usuario usa genes de ejemplo
    if (input$gene_source == "example") {
      data("input", package = "enrichR", envir = environment())
      genes <- get("input", envir = environment())
    } else {
      # Si el usuario usa genes personalizados
      genes <- character(0)
      
      # Leer genes desde archivo
      if (!is.null(input$gene_file)) {
        genes <- c(genes, readLines(input$gene_file$datapath, warn = FALSE))
      }
      
      # Leer genes desde texto manual
      if (nzchar(input$manual_genes)) {
        manual_genes <- unlist(strsplit(input$manual_genes, "[,\n\r\t; ]+"))
        genes <- c(genes, manual_genes)
      }
    }
    
    # Limpiar: quitar espacios, vacíos y duplicados
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    genes <- unique(genes)
    
    genes
  })
  
  
  # -------------------------
  # 4. Mostrar resumen de genes cargados
  # -------------------------
  output$gene_source_info <- renderText({
    source_text <- ifelse(
      input$gene_source == "example",
      "Using EnrichR example genes:",
      "Using custom input:"
    )
    
    paste(source_text, length(genes_reactive()), "genes loaded.")
  })
  
  
  # -------------------------
  # 5. Mostrar tabla de genes
  # -------------------------
  output$gene_table <- renderDT({
    datatable(
      data.frame(Gene = genes_reactive()),
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })
  
  
  # -------------------------
  # 6. Ejecutar enriquecimiento al pulsar el botón
  # -------------------------
  enrich_results <- eventReactive(input$run_enrich, {
    req(length(genes_reactive()) > 0)
    req(input$selected_dbs)
    
    tryCatch(
      enrichr(genes_reactive(), input$selected_dbs),
      error = function(e) NULL
    )
  })
  
  
  # -------------------------
  # 7. Crear pestañas dinámicas de resultados
  # -------------------------
  
  output$result_tabs <- renderUI({
    res <- enrich_results()
    req(res)
    
    tabs <- lapply(names(res), function(db_name) {
      tabPanel(
        title = db_name,
        DTOutput(paste0("table_", make.names(db_name)))
      )
    })
    
    do.call(tabsetPanel, tabs)
  })
  
  
  # -------------------------
  # 8. Renderizar tablas de resultados para cada base de datos
  # -------------------------
  
  observe({
    res <- enrich_results()
    req(res)
    
    lapply(names(res), function(db_name) {
      output[[paste0("table_", make.names(db_name))]] <- renderDT({
        df <- res[[db_name]]
        
        if (is.null(df) || nrow(df) == 0) {
          return(datatable(data.frame(Message = "No results available")))
        }
        
        datatable(
          df,
          options = list(pageLength = 8, scrollX = TRUE)
        )
      })
    })
  })
  
  
  # -------------------------
  # 9. Preparar datos para la gráfica
  # -------------------------
  plot_data <- reactive({
    res <- enrich_results()
    req(res)
    
    # Se usa la primera base de datos seleccionada para la visualización
    first_db <- names(res)[1]
    df <- res[[first_db]]
    
    validate(
      need(!is.null(df) && nrow(df) > 0, "No results to display.")
    )
    
    df %>%
      mutate(
        # Extrae el número de genes desde la columna Overlap
        GeneCount = as.numeric(sub("/.*", "", Overlap)),
        
        # Acorta nombres de términos muy largos
        ShortTerm = ifelse(
          nchar(Term) > 55,
          paste0(substr(Term, 1, 55), "..."),
          Term
        )
      ) %>%
      arrange(P.value) %>%
      slice_head(n = input$top_n) %>%
      mutate(
        ShortTerm = factor(ShortTerm, levels = rev(ShortTerm))
      )
  })
  
  
  # -------------------------
  # 10. Mostrar gráfico de enriquecimiento
  # -------------------------
  
  output$enrich_plot <- renderPlot({
    res <- enrich_results()
    req(res)
    
    first_db <- names(res)[1]
    df_plot <- plot_data()
    
    ggplot(df_plot, aes(x = GeneCount, y = ShortTerm, fill = P.value)) +
      geom_col(width = 0.82, color = "white", linewidth = 0.2) +
      scale_fill_gradientn(
        colours = c("#F46D43", "#D948A1", "#4C6EDB"),
        name = expression(italic("P value"))
      ) +
      labs(
        title = "Enrichment Analysis by Enrichr",
        subtitle = paste("Top terms from", first_db),
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
        axis.title = element_text(color = "#223447", face = "bold"),
        axis.text = element_text(color = "#35506A"),
        plot.title = element_text(face = "bold", size = 16, color = "#1E3A5F"),
        plot.subtitle = element_text(color = "#5A748C", size = 11)
      )
  })
}


# =========================
# EJECUTAR APP
# =========================

shinyApp(ui = ui, server = server)

