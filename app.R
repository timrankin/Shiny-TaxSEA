#
# TBA
#

library(shiny)
library(tidyverse)
library(ggrepel)
library(bslib)
library(bsicons)
library(openxlsx2)
library(DT)
library(TaxSEA)

# Helper to display P values and FDR values in scientific notation, without changing the dataframe values to text formatted in scientific notation.
# At the time of writing, there is an open issue for this functionality in DT: https://github.com/rstudio/DT/issues/938
jsExponential <- c(
  "function(row, data, displayNum, index){",
  "  for (var i = 3; i <= 4; i++) {",
  "    var x = data[i];",
  "    if (!isNaN(parseFloat(x))) {",
  "      $('td:eq(' + i + ')', row).html(parseFloat(x).toExponential(2));",
  "   }",
  "  }",
  "}"
)

# Helper to detect header rows in supplied (or example) data
has_header <- function(file_path) {
  if (tools::file_ext(file_path) == "csv") {
    first_row <- read.csv(file_path, header = FALSE, nrows = 1, stringsAsFactors = FALSE)
  } else if (tools::file_ext(file_path) %in% c("xlsx", "xlsm")) {
    first_row <- read_xlsx(file_path, col_names = FALSE, rows = 1)
  }
  
  # Check if all columns in the first row are numeric, except for the first column
  is_numeric <- sapply(first_row[-1], is.numeric)
  
  if(all(is_numeric)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

valid_input_format <- function(suppliedData) {
  if (ncol(suppliedData) != 4) {
    showModal(modalDialog(
      title = "Input Error",
      "Incorrect number of input columns. Expecting exactly 4; Taxa, log 2-fold changes, P value, and Padj or FDR."
    ))
    return(FALSE)
  } else if (!is.character(suppliedData[[1]])) {
    showModal(modalDialog(
      title = "Input Error",
      "First column (Taxa) must be text"
    ))
    return(FALSE)
  } else if (!is.numeric(suppliedData[[2]])) {
    showModal(modalDialog(
      title = "Input Error",
      "Second column (log 2-fold change) must be numeric"
    ))
    return(FALSE)
  } else if (!is.numeric(suppliedData[[3]]) || min(suppliedData[[3]]) < 0 || max(suppliedData[[3]]) > 1) {
    showModal(modalDialog(
      title = "Input Error",
      "Third column (P value) must be a numeric value between 0 and 1"
    ))
    return(FALSE)
  } else if (!is.numeric(suppliedData[[4]]) || min(suppliedData[[4]]) < 0 || max(suppliedData[[4]]) > 1) {
    showModal(modalDialog(
      title = "Input Error",
      "Fourth column (Padj / FDR) must be a numeric value between 0 and 1."
    ))
    return(FALSE)
  } else {
    return(TRUE)
  }
}

ui <- page_sidebar(
  title = "Shiny TaxSEA",
  sidebar = sidebar(
    width = "310px",
    fileInput(
      "file",
      label = tooltip(
        trigger = list(
          "Differential abundance data for analysis",
          bs_icon("question-circle")
        ), "Browse for a Microsoft Excel️ (R) or CSV file with columns for: Taxa, log 2-fold change, P value, and Padj / FDR"
      ) 
    ),
    
    selectInput(
      "database_selection",
      "Select TaxSEA results to display",
      choices = list("Metabolites" = "Metabolite_producers", "Health Associations" = "Health_associations", "BugSigDB" = "BugSigDB")
    ),
    
    actionButton(
      "loadSampleData",
      icon = icon("bolt"),
      "Analyse sample data",
    ),
    downloadButton(
      "downloadSampleData",
      icon = icon("cloud-download"),
      "Download sample data"
    ),
    div(
      class = "text-center",
      tags$a(
        href = "https://github.com/timrankin/Shiny-TaxSEA/issues",
        target = "_blank",
        class = "d-flex align-items-center justify-content-center text-decoration-none",
        span("Report a bug", class = "me-2"),
        bs_icon("bug-fill")
      )
    )
  ),
  
  layout_columns(card(
    card_header(
      class = "d-flex align-items-center",
      "Bar Plot",
      tags$span(style = "margin-left: 5px"),
      tooltip(
        bs_icon("info-circle"),
        "Select up to 8 Taxon Sets below to display in this plot. If no selection is made, the top 8 Taxon Sets by -log10 FDR value will be displayed"
      ),
      uiOutput("downloadBarPlotUi"),
    ),
    plotOutput("barPlot")
  ),
  card(
    card_header(
      class = "d-flex align-items-center",
      "Volcano Plot",
      tags$span(style = "margin-left: 5px"),
      tooltip(
        bs_icon("info-circle"),
        "The last (most recent) selection you make below will be used to title the plot and highlight taxon set members of interest. If no selection is made, the top result by -log10 FDR will be displayed"
      ),
      uiOutput("downloadVolcanoPlotUi"),
    ),
    plotOutput("volcanoPlot")
  )), card(
    card_header(
      class = "d-flex align-items-center",
      "TaxSEA Results",
      uiOutput("downloadDataTableUi"),
      # actionButton(
      #   "downloadTable",
      #   icon = icon("cloud-download"),
      #   "Download",
      #   class = "btn-sm btn-primary ms-auto"
      # )
    ),
    DTOutput("table")
  )
)

server <- function(input, output, session) {
  # Notifications variable, used if > 8 rows are selected in table
  notificationIds <- NULL
  
  # Reactive value to store data from either user supplied file or example data
  data <- reactiveVal(NULL)
  
  barPlotReady <- reactiveVal(FALSE)
  volcanoPlotReady <- reactiveVal(FALSE)
  dataTableReady <- reactiveVal(FALSE)
  
  # TODO: implementing a check to display a warning if there are no taxon sets w/FDR <0.02. This probably belongs better in the TaxSEA function, otherwise DT
  
  # Read uploaded file
  observeEvent(input$file, {
    req(input$file)
    
    if (tools::file_ext(input$file$datapath) == "csv") {
      suppliedData <- read.csv(input$file$datapath, header = has_header(input$file$datapath))
      
    } else if (tools::file_ext(input$file$datapath) %in% c("xlsx", "xlsm")) {
      suppliedData <- read_xlsx(input$file$datapath, col_names = has_header(input$file$datapath))
    }
    
    if(valid_input_format(suppliedData)) {
      return(data(suppliedData))
    } else {
      return()
    }
  })

  # Handle example data button click & load example data
  observeEvent(input$loadSampleData, {
    exampleData <- read.csv("test_input.csv", header = has_header("test_input.csv"))
    data(exampleData)
  })
    
  # Run TaxSEA on data
  taxseaResults <- reactive({
    # Make sure data has been supplied in a valid format
    req(data())
    
    # Get taxon ranks from user supplied data
    taxonRanks <- setNames(data()[[2]], data()[[1]])
    
    results <- TaxSEA(taxonRanks)
    
    # Drop column 5 from all results.
    results$Metabolite_producers <- results$Metabolite_producers[, -5]
    results$Health_associations <- results$Health_associations[, -5]
    results$BugSigDB <- results$BugSigDB[, -5]
    
    # Make results presentable
    # TODO: Only apply this function to disease & metabolites - bsdb requires its own to remove hyphens etc.
    results <- lapply(results, function(df) {
      # Remove rownames
      rownames(df) <- NULL
      
      # Change to readable column names
      colnames(df) <- c("Taxon Set", "Median rank of set members", "P Value", "FDR", "Taxon Set Members")
      
      # Make taxon set names presentable
      df$`Taxon Set` <- df$`Taxon Set` %>%
        str_replace("_(.)", ": \\1") %>%
        str_replace("(: )(.)", toupper) %>%
        str_replace_all("_", " ")
      
      # Make taxon set members presentable
      df$`Taxon Set Members` <- df$`Taxon Set Members` %>%
        str_replace_all("_", " ")
      
      return(df)
    })
    return(results)
  })
  
  # 'Debounce' the clicks on table rows, so we don't redraw the plots too frequently
  debounced_selection <- debounce(
    reactive(input$table_rows_selected),
    millis = 500  # 500 ms delay
  )
  
  # Render the bar plot
  barPlot <- reactive({
    # Make sure data has been supplied and in a valid format
    req(taxseaResults())
    
    # Store selected rows, even if this is 0
    selected_rows <- debounced_selection()
    
    if (is.null(selected_rows)) {
      # Plot top 6 results based on -log10 FDR values
      dataForPlot <- taxseaResults()[[input$database_selection]] %>%
        mutate(negativeLog10FDR = -log10(taxseaResults()[[input$database_selection]][[4]])) %>%
        arrange(desc(negativeLog10FDR)) %>%
        slice_head(n = 8)
    } else if (length(selected_rows) <= 8) {
      # Plot the selected rows
      dataForPlot <- taxseaResults()[[input$database_selection]][selected_rows, ] %>%
        mutate(negativeLog10FDR = -log10(taxseaResults()[[input$database_selection]][selected_rows, ][[4]])) %>%
        arrange(desc(negativeLog10FDR))
    } else {
      # Display the first 8 user selected taxon sets, together with a warning
      dataForPlot <- taxseaResults()[[input$database_selection]][selected_rows[1:8], ] %>%
        mutate(negativeLog10FDR = -log10(taxseaResults()[[input$database_selection]][selected_rows[1:8], ][[4]])) %>%
        arrange(desc(negativeLog10FDR))
      
      notificationId <<- showNotification(
        "⚠️ Bar plot limited to a max of 8 Taxon Sets. Your first 8 selections are displayed, deselect some in order to add more.",
        duration = 10,
        closeButton = TRUE,
        type = "warning"
      )
    }
    
    plot <- ggplot(dataForPlot, aes(x = negativeLog10FDR, y = reorder(str_wrap(`Taxon Set`, width = 34), negativeLog10FDR))) +
      geom_col(fill = "#00aedb") +
      labs(
        x = expression(-log[10] ~ FDR),
        y = "Taxon Sets"
      ) +
      geom_vline(xintercept = -log10(0.1), linetype = 5) +
      theme_classic() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
      )
    barPlotReady(TRUE)
    return(plot)
  })
  
  # Output the bar plot
    output$barPlot <- renderPlot({
      barPlot()
    })
  
  # Render the volcano plot
    volcanoPlot <- reactive({
      req(data())
      req(taxseaResults())
      
      dataForPlot <- data()
      
      # Store last selected row, it will be NULL if none are selected
      last_row_selected <- debounced_selection()[length(debounced_selection())]
      
      # If no selection is made, take the taxon set members from the top TaxSEA result. Otherwise, use the last selection.
      if (is.null(last_row_selected)) {
        taxa_of_interest <- unlist(strsplit(taxseaResults()[[input$database_selection]][[5]][1], ", "))
        plot_title <- taxseaResults()[[input$database_selection]][[1]][1]
      } else {
        taxa_of_interest <- unlist(strsplit(taxseaResults()[[input$database_selection]][[5]][last_row_selected], ", "))
        plot_title <- taxseaResults()[[input$database_selection]][[1]][last_row_selected]
      }
      
      # TODO: Support hyphens as well as underscores (input may also be supplied with spaces already)
      taxa_of_interest <- str_replace(taxa_of_interest, " ", "_")
      
      dataForPlot$is_of_interest <- dataForPlot[[1]] %in% taxa_of_interest
      
      label_data <- dataForPlot[dataForPlot$is_of_interest != FALSE &
                                  dataForPlot[[3]] < 0.05, ]
      label_data$Taxa <- gsub("_", " ", label_data$Taxa)
      label_data$Abbreviated_taxa <- sub("^([A-Za-z])[a-z]+\\s", "\\1. ", label_data$Taxa)
      
      plot <- ggplot(
        dataForPlot,
        aes(
          x = dataForPlot[, 2],
          y = -log10(dataForPlot[, 3]),
          color = is_of_interest,
          alpha = is_of_interest
        )
      ) +
        geom_point(aes(size = is_of_interest)) +
        theme_classic() +
        labs(title = plot_title,
             x = "Input Ranks",
             y = expression(-log[10] ~ FDR)) +
        scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
        scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue")) +
        scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
        geom_text_repel(
          data = label_data,
          aes(
            x = label_data[, 2],
            y = -log10(label_data[, 3]),
            label = Abbreviated_taxa,
            fontface = "italic"
          ),
          size = 5
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.position = "none"
        ) +
        geom_vline(xintercept = 0, linetype = 5)
      
      volcanoPlotReady(TRUE)
      return(plot)
  })
  
    # Output the volcano plot
    output$volcanoPlot <- renderPlot({
      volcanoPlot()
    })
    
  # Render the data table
  output$table = DT::renderDataTable({
    req(taxseaResults())
    updateActionButton(session, "downloadTable", label = "Download", disabled = FALSE)
    table <- datatable(
      taxseaResults()[[input$database_selection]],
      filter = "top",
      options = list (
        searching = FALSE,
        paging = FALSE,
        columnDefs = list(
          list(
            targets = 2, width = "125px"
          ),
          list(
            targets = c(3, 4), width = "70px"
          )
        ),
        rowCallback = JS(jsExponential),
        selection = 'multiple'
      )
    ) %>%
      formatRound(
        columns = "Median rank of set members",
        digits = 3 # TODO: Consider whether it's best to reference the textual column title or by number
        ) %>% 
      formatStyle(
        columns = 5,
        fontStyle = "italic"
      )
    
    dataTableReady(TRUE)
    return(table)
  })
  
  # Handle sample data download
  output$downloadSampleData <- downloadHandler(
    filename = function() {"Shiny-TaxSEA_sample_data.csv"},
    content = function(file) {
      file.copy("Shiny-TaxSEA_sample_data.csv", file)
    }
  )
  
  # Render bar plot download button only when ready
  output$downloadBarPlotUi <- renderUI({
    req(barPlotReady())
    downloadButton(
      "downloadBarPlot",
      icon = icon("cloud-download"),
      "Download",
      class = "btn-sm btn-primary ms-auto"
    )
  })
  
  # Handle bar plot download
  output$downloadBarPlot <- downloadHandler(
    filename = function() {
      paste0("Shiny-TaxSEA_", gsub("_(.)", "_\\U\\1", input$database_selection, perl = TRUE), "_Bar_Plot_", format(Sys.time(), "%Y%m%d_%H%M%S"),".png")
    },
    content = function(file){
      png(file, width = 800, height = 400)
      print(barPlot())
      dev.off()
    }
  )
  
  # Render volcano plot download button only when ready
  output$downloadVolcanoPlotUi <- renderUI({
    req(volcanoPlotReady())
    downloadButton(
      "downloadVolcanoPlot",
      icon = icon("cloud-download"),
      "Download",
      class = "btn-sm btn-primary ms-auto"
    )
  })
  
  # Handle bar plot download
  output$downloadVolcanoPlot <- downloadHandler(
    filename = function() {
      paste0("Shiny-TaxSEA_", gsub("_(.)", "_\\U\\1", input$database_selection, perl = TRUE), "_Volcano_Plot_", format(Sys.time(), "%Y%m%d_%H%M%S"),".png")
    },
    content = function(file){
      png(file, width = 800, height = 400)
      print(volcanoPlot())
      dev.off()
    }
  )
  
  # Render table download button only when ready
  output$downloadDataTableUi <- renderUI({
    req(dataTableReady())
    downloadButton(
      "downloadTable",
      icon = icon("cloud-download"),
      "Download",
      class = "btn-sm btn-primary ms-auto"
    )
  })
  
  # Handle table download
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste0("Shiny-TaxSEA_", gsub("_(.)", "_\\U\\1", input$database_selection, perl = TRUE), "_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(taxseaResults()[[input$database_selection]], file, row.names = FALSE)
    },
    # contentType = "text/csv"
  )
}

# Run the app
shinyApp(ui = ui, server = server)