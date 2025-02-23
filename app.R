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
# There is an open issue for this functionality in DT: https://github.com/rstudio/DT/issues/938
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
      "Incorrect number of input columns. Expecting exactly 4; Taxa, log 2-fold change, P value, and Padj or FDR."
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
  title = div (
    class = "d-flex justify-content-between align-items-center w-100",
    span("Shiny TaxSEA"),
    tags$a(
      href = "https://github.com/timrankin/Shiny-TaxSEA/issues",
      target = "_blank",
      class = "d-flex align-items-center text-decoration-none",
      span("Report a bug", class = "me-2"),
      bs_icon("bug-fill"),
    )
  ),
  
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
      choices = list("Metabolites" = "Metabolite_producers", "Health Associations" = "Health_associations", "BugSigDB" = "BugSigdB")
    ),
    
    actionButton(
      "loadExample",
      icon = icon("bolt"),
      "Analyse sample data",
    ),
    actionButton(
      "downloadExample",
      icon = icon("cloud-download"),
      "Download sample data",
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
      actionButton(
        "downloadBarPlot",
        icon = icon("cloud-download"),
        "Download",
        class = "btn-sm btn-primary ms-auto"
      )
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
      actionButton(
        "downloadVolcanoPlot",
        icon = icon("cloud-download"),
        "Download",
        class = "btn-sm btn-primary ms-auto"
      )
    ),
    plotOutput("volcanoPlot")
  )), card(
    card_header(
      class = "d-flex align-items-center",
      "TaxSEA Results",
      # TODO: Hide / grey out the button until we have data to download.
      actionButton(
        "downloadTable",
        icon = icon("cloud-download"),
        "Download",
        class = "btn-sm btn-primary ms-auto"
      )
    ),
    DTOutput("table")
  )
)

server <- function(input, output, session) {
  # Disable all downlolad buttons initially
  updateActionButton(session, "downloadTable", label = "Download", disabled = TRUE)
  updateActionButton(session, "downloadBarPlot", label = "Download", disabled = TRUE)
  updateActionButton(session, "downloadVolcanoPlot", label = "Download", disabled = TRUE)
  
  # Notifications variable, used if > 8 rows are selected in table
  notificationIds <- NULL
  
  # Reactive value to store data from either user supplied file or example data
  data <- reactiveVal(NULL)
  
  is_sufficiently_enriched <- reactiveVal(NULL)
  
  # Handle example data button click & load example data
  observeEvent(input$loadExample, {
    exampleData <- read.csv("test_input.csv", header = has_header("test_input.csv"))
    data(exampleData)
  })
  
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
  
  # Run TaxSEA on data
  taxseaResults <- reactive({
    # Make sure data has been supplied in a valid format
    req(data())
    
    # Get taxon ranks from user supplied data
    taxonRanks <- setNames(data()[[2]], data()[[1]])
    
    results <- TaxSEA(taxonRanks)
    
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
  
  
  # TODO: Remove - it isn't needed.
  # 'Debounce' the clicks on table rows, so we don't redraw the plots too frequently
  debounced_selection <- debounce(
    reactive(input$table_rows_selected),
    millis = 0  # 500 ms delay
  )
  
  
  # Render the bar plot
  output$barPlot <- renderPlot({
    # Make sure data has been supplied and in a valid format
    req(taxseaResults())
    
    # Store selected rows, even if this is 0
    selected_rows <- debounced_selection()
    
    if(is.null(selected_rows)) {
      # Plot top 6 results based on -log10 FDR values
      dataForPlot <- taxseaResults()[[input$database_selection]] %>%
        mutate(negativeLog10FDR = -log10(taxseaResults()[[input$database_selection]][[4]])) %>%
        arrange(desc(negativeLog10FDR)) %>%
        slice_head(n = 8)
    } else if (length(selected_rows) <=8) {
      # Plot the selected rows
      dataForPlot <- taxseaResults()[[input$database_selection]][selected_rows, ] %>%
        mutate(negativeLog10FDR = -log10(taxseaResults()[[input$database_selection]][selected_rows, ][[4]])) %>%
        arrange(desc(negativeLog10FDR))
    } else {
      # Display the first 8 the user selected if possible, together with a warning
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
    
    ggplot(dataForPlot, aes(x = negativeLog10FDR, y = reorder(str_wrap(`Taxon Set`, width = 34), negativeLog10FDR))) +
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
  })
  
  # Render the volcano plot
  output$volcanoPlot <- renderPlot({
    req(data())
    req(taxseaResults())

    dataForPlot <- data()
    
    # Store last selected row, it will be NULL if none are selected
    last_row_selected <- debounced_selection()[length(debounced_selection())]
    
    # Debugging
    print(last_row_selected)
    
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
    
    label_data <- dataForPlot[dataForPlot$is_of_interest != FALSE & dataForPlot[[3]] < 0.05,]
    label_data$Taxa <- gsub("_", " ", label_data$Taxa)
    label_data$Abbreviated_taxa <- sub("^([A-Za-z])[a-z]+\\s", "\\1. ", label_data$Taxa)
    
    # Debugging
    print(label_data)
    print(head(dataForPlot))

    ggplot(dataForPlot, aes(x = dataForPlot[,2], y = -log10(dataForPlot[,3]), color = is_of_interest, alpha=is_of_interest)) +
      geom_point(aes(
        size = is_of_interest
      )) + 
      theme_classic() +
      labs(
        title = plot_title,
        x = "Input Ranks",
        y = expression(-log[10] ~ FDR)
      ) +
      scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
      scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue")) +
      scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
      geom_text_repel(data = label_data,
                      aes(
                        x = label_data[,2],
                        y = -log10(label_data[,3]),
                        label = Abbreviated_taxa,
                        fontface = "italic"
                        ),
                      size = 5) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "none") +
        geom_vline(xintercept = 0, linetype = 5)
  })
  
  # Render the data table
  output$table = DT::renderDataTable({
    req(taxseaResults())
    datatable(
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
  })
  

  # Handle downloads
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste(input$database_selection, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(taxseaResults()[[input$database_selection]], file, row.names = FALSE)
    }
  )
}

# JavaScript to disable and enable buttons
jsCode <- "
Shiny.addCustomMessageHandler('disable_button', function(id) {
  document.getElementById(id).setAttribute('disabled', 'true');
});

Shiny.addCustomMessageHandler('enable_button', function(id) {
  document.getElementById(id).removeAttribute('disabled');
});
"

# Include JS in UI
ui <- tagList(
  tags$script(jsCode),
  ui
)

# Run the application 
shinyApp(ui = ui, server = server)