#
# TBA
#

library(shiny)
library(tidyverse)
library(bslib)
library(bsicons)
library(DT)
library(TaxSEA)

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

ui <- page_sidebar(
  title = "Shiny TaxSEA",
  
  sidebar = sidebar(
    fileInput(
      "file",
      label = tooltip(
        trigger = list(
          "Data for analysis",
          bs_icon("question-circle")
        ), "This is a Microsoft Excel ®️ or csv file with columns for: Taxa, log 2-fold change, P value, and Padj / FDR"
      ) 
    ),
    
    checkboxInput(
      "hasHeaders",
      "The supplied data has headers",
      value = FALSE
    ),
    
    selectInput(
      "database_selection",
      "Select results to display in the table",
      choices = list("Metabolites" = "Metabolite_producers", "Health Associations" = "Health_associations", "BugSigDB" = "BugSigdB")
    ),
    
    hr(),
    
    "TODO: Taxa of interest, a searchable drop-down kinda widget would be great here",
    
    actionButton(
      "exampleData",
      "Load example data"
    )
  ),
  
  layout_columns(
    card(
      card_header("Plot 1"),
      plotOutput("barPlot")
    ),
    card(
      card_header("Plot 2"),
      plotOutput("volcanoPlot")
    )
  ),
  
  card(
    card_header(
      class = "d-flex align-items-center",
      "Table",
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
  # Read uploaded file
  suppliedData <- reactive({
    req(input$file)
    data <- read.csv(input$file$datapath, header = input$hasHeaders)
    
    # Validate the structure of the data
    if (ncol(data) != 4 ||
        !is.character(data[[1]]) ||
        !is.numeric(data[[2]]) ||
        !is.numeric(data[[3]]) ||
        !is.numeric(data[[4]])) {
      # TODO: Rather than crashing, turn this into a warning message in the UI.
      stop("The supplied file must have 4 columns: Taxa, log 2-fold change, P value, and Padj or FD")
    }
    
    return(data)
  })
  
  # Run TaxSEA on data
  taxseaResults <- reactive({
    # Make sure data has been supplied in a valid format
    req(suppliedData())
    
    # TODO: catch errors and display something meaningful in the UI
    
    # Get taxon ranks from user supplied data
    taxonRanks <- setNames(suppliedData()[[2]], suppliedData()[[1]])
    
    results <- TaxSEA(taxonRanks)
    
    # Make results presentable
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
  
  # Render the histogram
  output$barPlot <- renderPlot({
    # Make sure data has been supplied in a valid format
    req(taxseaResults())
    
    dataForPlot <- taxseaResults()$Metabolite_producers %>%
      mutate(negativeLog10PValue = -log10(taxseaResults()$Metabolite_producers[[3]])) %>%
      mutate(`Taxon Set` = factor(`Taxon Set`, levels = `Taxon Set`[order(negativeLog10PValue, decreasing = FALSE)])) %>%
      arrange(desc(negativeLog10PValue)) %>%
      slice_head(n = 6)
    
    ggplot(dataForPlot, aes(x = negativeLog10PValue, y = `Taxon Set`)) +
      geom_col(fill = "steelblue") +
      labs(
        title = "Top 6 -log10 P values",
        x = "-log10 P value",
        y = "gutMGene taxon sets"
      ) +
      theme_minimal()
    
  })
  
  # TODO: Render the volcano plot
  
  # output$volcanoPlot < renderPlot({
  #   
  # })
  
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
        rowCallback = JS(jsExponential)
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
  
  # observeEvent(
  #   
  # )
  
  # Load example data
  
  
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
