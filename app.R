#
# TBA
#

library(shiny)
library(shinyjs)
library(tidyverse)
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
        ), "Browse for a Microsoft Excel ®️ or csv file with columns for: Taxa, log 2-fold change, P value, and Padj / FDR"
      ) 
    ),
    
    # TODO: Delete if auto detection works well
    # checkboxInput(
    #   "hasHeaders",
    #   "The supplied data has headers",
    #   value = FALSE
    # ),
    
    selectInput(
      "database_selection",
      "Select TaxSEA results to display",
      choices = list("Metabolites" = "Metabolite_producers", "Health Associations" = "Health_associations", "BugSigDB" = "BugSigdB")
    ),
    
    actionButton(
      "loadExample",
      "Load example data"
    )
  ),
  
  layout_columns(
    card(
      card_header(
        "Bar Plot",
        tooltip(
          bs_icon("info-circle"),
          "Select up to 8 Taxon Sets below to display in this plot. If no selection is made, the top 8 Taxon Sets by -log10 FDR value will be displayed"
        )),
      plotOutput("barPlot")
    ),
    card(
      card_header(
        "Volcano Plot",
        tooltip(
          bs_icon("info-circle"),
          "The last selection you make below will be used to title the plot and highlight taxon set members of interest. If no selection is made, the top result by -log10 FDR will be displayed"
        )),
      plotOutput("volcanoPlot")
    )
  ),
  
  card(
    card_header(
      class = "d-flex align-items-center",
      "TaxSEA Results",
      downloadButton(
        "downloadButton",
        "Download",
        class = "btn-sm btn-primary ms-auto"
      )
    ),
    # "Table goes here"
    DTOutput("table")
  )
)

server <- function(input, output, session) {
  # Disable download button until there are results available to download
  shinyjs::disable("downloadButton")
  
  # Notifications variable, used if > 8 rows are selected in table
  notificationIds <- NULL
  
  # Reactive value to store data from either user supplied file or example data
  data <- reactiveVal(NULL)
  
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
      
      if (ncol(suppliedData) != 4 ||
          !is.character(suppliedData[[1]]) ||
          !is.numeric(suppliedData[[2]]) ||
          !is.numeric(suppliedData[[3]]) ||
          !is.numeric(suppliedData[[4]])) {
        # TODO: Rather than crashing, turn this into a warning message in the UI.
        stop("The supplied file must have 4 columns: Taxa, log 2-fold change, P value, and Padj or FD")
      }
      
      data(suppliedData)
    } else if (tools::file_ext(input$file$datapath) %in% c("xlsx", "xlsm")) {
      suppliedData <- read_xlsx(input$file$datapath, col_names = has_header(input$file$datapath))
    }
    
    data(suppliedData)
  })
  
  # Run TaxSEA on data
  taxseaResults <- reactive({
    # Make sure data has been supplied in a valid format
    req(data())
    
    # TODO: catch errors and display something meaningful in the UI
    
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
    
    shinyjs::enable("downloadButton")
    
    return(results)
  })
  
  # TODO: see if this is actually necessary (it probably isn't)
  # 'Debounce' clicks on table rows, so we don't redraw the plots too frequently
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
    # selected_rows <- input$table_rows_selected
    
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
    
    ggplot(dataForPlot, aes(x = negativeLog10FDR, y = reorder(str_wrap(`Taxon Set`, width = 30), negativeLog10FDR))) +
      geom_col(fill = "steelblue") +
      labs(
        title = "-log10 FDR values", # @Feargal, is this needed?
        x = "-log10 FDR value",
        y = "Taxon Sets"
      ) +
      geom_vline(xintercept = -log10(0.1), linetype = 5) +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 12),
      )
  })
  
  # Render the volcano plot
  output$volcanoPlot <- renderPlot({
    req(data())

    dataForPlot <- data()

    # Store selected rows, it will be NULL if none are selected
    last_row_selected <- debounced_selection()[length(debounced_selection())]

    if (is.null(last_row_selected)) {
      # Taxon set members
    }

    dataForPlot$is_of_interest <- dataForPlot[[1]] %in% c("Bacteroides_vulgatus")

    ggplot(dataForPlot, aes(x = dataForPlot[,2], y = -log10(dataForPlot[,3]), color = is_of_interest)) +
             geom_point()+theme_minimal()
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
            targets = c(3, 4), width = "75px"
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
  output$downloadButton <- downloadHandler(
    filename = function() {
      paste(input$database_selection, "_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(taxseaResults()[[input$database_selection]], file, row.names = FALSE)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)