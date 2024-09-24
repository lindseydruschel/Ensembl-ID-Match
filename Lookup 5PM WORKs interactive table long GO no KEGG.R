
# Load necessary packages
if (!require("shiny")) install.packages("shiny")
library(shiny)
library(biomaRt)
library(httr)
library(jsonlite)
library(writexl)
library(DT)  # For interactive tables

# Available Ensembl versions for Rnor6.0 and mRatBN7.2 in the past 10 years
ensembl_versions <- c("Rnor6.0 (version 80)" = 80, 
                      "Rnor6.0 (version 90)" = 90, 
                      "Rnor6.0 (version 100)" = 100,
                      "Rnor6.0 (version 101)" = 101,
                      "mRatBN7.2 (version 102)" = 102,
                      "mRatBN7.2 (version 103)" = 103,
                      "mRatBN7.2 (version 104)" = 104)

# Define UI
ui <- fluidPage(
  
  # Add custom CSS to make the background pink
  tags$head(
    tags$style(HTML("
      body {
        background-color: #ffe6f0;
      }
    "))
  ),
  
  # App title
  titlePanel("Ensembl ID and Gene Symbol Search"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # Sidebar panel for user input
    sidebarPanel(
      # User selects the input type (Reordered as per request)
      radioButtons("input_type", "Select input type:",
                   choices = list("Gene Symbol" = "gene_symbol",
                                  "Current Ensembl ID" = "current_ensembl",
                                  "Outdated Ensembl ID" = "outdated_ensembl")),
      
      # Text input for comma-separated gene IDs
      textInput("user_input", "Enter comma-separated gene IDs:", ""),
      
      # Text input for the Excel file name
      textInput("file_name", "Enter file name for download:", "gene_query_results"),
      
      # Display instructions
      helpText("Enter multiple gene IDs or symbols separated by commas (e.g., 'TP53, BRCA1').")
    ),
    
    # Main panel to display output
    mainPanel(
      # Output progress bar
      uiOutput("progress_bar"),
      
      # Output interactive table with results
      DT::dataTableOutput("result_table"),
      
      # Download button
      downloadButton("download_data", "Download Results"),
      
      # Display additional gene information when a gene is clicked
      h3("Gene Information"),
      verbatimTextOutput("gene_info")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive value to store results
  results <- reactiveVal(data.frame())
  
  # Show modal only when "Outdated Ensembl ID" is selected
  observeEvent(input$input_type, {
    if (input$input_type == "outdated_ensembl") {
      showModal(modalDialog(
        title = "Ensembl Version Selection",
        "This app is currently running Rnor6.0 version 80. Do you want to continue with that reference database?",
        selectInput("ensembl_version", "Select Ensembl version:", choices = ensembl_versions),
        footer = tagList(
          actionButton("continue_ensembl", "Continue with Selected Version"),
          modalButton("Cancel")
        ),
        easyClose = TRUE
      ))
    }
  })
  
  # Store selected Ensembl version
  ensembl_version <- reactiveVal(80)  # Default version is 80
  
  # When the user clicks continue, update the Ensembl version
  observeEvent(input$continue_ensembl, {
    ensembl_version(input$ensembl_version)
    removeModal()
  })
  
  # Function to process the comma-separated input and perform the selected query
  output$result_table <- DT::renderDataTable({
    req(input$user_input)  # Ensure input is not empty
    
    # Split the input string into a list of gene IDs or symbols (remove spaces)
    gene_ids <- strsplit(gsub("\\s+", "", input$user_input), ",")[[1]]
    
    # Initialize progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Query in progress...", value = 0)
    on.exit(progress$close())  # Close progress bar when done
    
    # Number of queries to process
    total_queries <- length(gene_ids)
    
    # Conditional logic based on input type
    if (input$input_type == "outdated_ensembl") {
      # Use the selected Ensembl version to retrieve gene symbols
      selected_version <- ensembl_version()
      ensembl_v <- useEnsembl(biomart = "ensembl", version = selected_version, dataset = "rnorvegicus_gene_ensembl")
      
      results_old <- lapply(seq_along(gene_ids), function(i) {
        progress$inc(1 / total_queries, detail = paste("Processing", gene_ids[i]))
        result <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'go_id', 'name_1006', 'kegg_pathway_id', 'kegg_pathway_description'), 
                        filters = 'ensembl_gene_id', 
                        values = gene_ids[i], 
                        mart = ensembl_v)
        
        # Limit GO terms to top 5 and collapse into a single string
        go_terms <- paste(head(unique(result$name_1006), 5), collapse = ", ")
        result$GO_Terms <- ifelse(nzchar(go_terms), go_terms, NA)
        
        # Collapse KEGG pathways into a single string
        result$KEGG_Pathways <- paste(unique(result$kegg_pathway_description), collapse = ", ")
        
        # Reorder the columns
        result <- result[, c('ensembl_gene_id', 'external_gene_name', 'GO_Terms', 'KEGG_Pathways')]
        colnames(result) <- c('Ensembl_ID', 'Gene_Symbol', 'GO_Terms', 'KEGG_Pathways')  # Adjust column names
        return(result)
      })
      results_old <- do.call(rbind, results_old)
      results(results_old)  # Store results for download
      return(results_old)
      
    } else if (input$input_type == "current_ensembl") {
      # Use the latest Ensembl version to map to new Ensembl IDs
      ensembl_latest <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
      
      # Query for gene symbols based on current Ensembl IDs
      results_new <- lapply(seq_along(gene_ids), function(i) {
        progress$inc(1 / total_queries, detail = paste("Processing", gene_ids[i]))
        result <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'go_id', 'name_1006', 'kegg_pathway_id', 'kegg_pathway_description'), 
                        filters = 'ensembl_gene_id', 
                        values = gene_ids[i], 
                        mart = ensembl_latest)
        
        # Limit GO terms to top 5 and collapse into a single string
        go_terms <- paste(head(unique(result$name_1006), 5), collapse = ", ")
        result$GO_Terms <- ifelse(nzchar(go_terms), go_terms, NA)
        
        # Collapse KEGG pathways into a single string
        result$KEGG_Pathways <- paste(unique(result$kegg_pathway_description), collapse = ", ")
        
        # Reorder the columns
        result <- result[, c('ensembl_gene_id', 'external_gene_name', 'GO_Terms', 'KEGG_Pathways')]
        colnames(result) <- c('Ensembl_ID', 'Gene_Symbol', 'GO_Terms', 'KEGG_Pathways')  # Adjust column names
        return(result)
      })
      results_new <- do.call(rbind, results_new)
      results(results_new)  # Store results for download
      return(results_new)
      
    } else if (input$input_type == "gene_symbol") {
      gene_info_results <- do.call(rbind, lapply(seq_along(gene_ids), function(i) {
        progress$inc(1 / total_queries, detail = paste("Processing", gene_ids[i]))
        gene_symbol <- gene_ids[i]
        url <- paste0("https://mygene.info/v3/query?q=", gene_symbol, 
                      "&fields=ensembl.gene,alias,symbol,go,KEGG&species=rat")
        response <- GET(url)
        content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
        
        ensembl_id <- if (!is.null(content$hits$ensembl$gene)) content$hits$ensembl$gene[1] else NA
        aliases <- if (!is.null(content$hits$alias)) paste(content$hits$alias, collapse = ", ") else NA
        
        # Extract top 5 GO terms (if available)
        go_terms <- if (!is.null(content$hits$go$BP)) paste(head(sapply(content$hits$go$BP, `[[`, "term"), 5), collapse = ", ") else NA
        
        # Extract KEGG pathways
        kegg_pathways <- if (!is.null(content$hits$KEGG)) paste(content$hits$KEGG, collapse = ", ") else NA
        
        result <- data.frame(Gene_Symbol = gene_symbol, Ensembl_ID = ensembl_id, Aliases = aliases, GO_Terms = go_terms, KEGG_Pathways = kegg_pathways)
        
        return(result)
      }))
      results(gene_info_results)  # Store results for download
      return(gene_info_results)
    }
  })
  
  # Handle gene selection from the table and display additional info
  observeEvent(input$result_table_rows_selected, {
    selected_row <- input$result_table_rows_selected
    if (length(selected_row) > 0) {
      selected_gene <- results()[selected_row, "external_gene_name"]
      
      # Query gene info from mygene.info or Ensembl
      if (!is.null(selected_gene)) {
        url <- paste0("https://mygene.info/v3/query?q=", selected_gene, "&fields=ensembl.gene,alias,symbol,summary")
        response <- GET(url)
        content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
        
        # Check if content contains the necessary fields
        if (!is.null(content$hits) && length(content$hits) > 0 && !is.null(content$hits[[1]])) {
          
          # Handle aliases: Remove NULLs and only show unique aliases
          alias_data <- content$hits[[1]]$alias
          unique_aliases <- if (!is.null(alias_data)) unique(alias_data[!is.null(alias_data)]) else "No aliases available"
          
          # Handle summary: Check if summary is present
          summary_info <- if (!is.null(content$hits[[1]]$summary)) content$hits[[1]]$summary else "No summary available"
          
          gene_info <- paste0(
            "Gene: ", selected_gene, "\n",
            "Aliases: ", if (length(unique_aliases) > 0) paste(unique_aliases, collapse = ", ") else "No aliases available", "\n",
            "Summary: ", summary_info
          )
          
        } else {
          gene_info <- "No information found for this gene."
        }
        
        # Display the gene info in the UI
        output$gene_info <- renderText(gene_info)
        
      } else {
        output$gene_info <- renderText("No gene selected.")
      }
    }
  })
  
  # Download handler for results
  output$download_data <- downloadHandler(
    filename = function() {
      paste(input$file_name, ".xlsx", sep = "")  # Use the custom name provided by the user
    },
    content = function(file) {
      write_xlsx(results(), file)  # Save results as an Excel file
    }
  )
  
  # Display progress bar
  output$progress_bar <- renderUI({
    tags$div(style = "margin-top: 10px;")
  })
}

# Run the application
shinyApp(ui = ui, server = server)
