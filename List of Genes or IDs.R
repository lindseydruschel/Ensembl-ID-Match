# Load necessary packages
if (!require("shiny")) install.packages("shiny")
library(shiny)
library(biomaRt)
library(httr)
library(jsonlite)

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
      # User selects the input type
      radioButtons("input_type", "Select input type:",
                   choices = list("Outdated Ensembl ID" = "outdated_ensembl", 
                                  "Current Ensembl ID" = "current_ensembl",
                                  "Gene Symbol" = "gene_symbol")),
      
      # Text input for comma-separated gene IDs
      textInput("user_input", "Enter comma-separated gene IDs:", ""),
      
      # Display instructions
      helpText("Enter multiple gene IDs or symbols separated by commas (e.g., 'TP53, BRCA1').")
    ),
    
    # Main panel to display output
    mainPanel(
      # Output text
      tableOutput("result")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Function to process the comma-separated input and perform the selected query
  output$result <- renderTable({
    req(input$user_input)  # Ensure input is not empty
    
    # Split the input string into a list of gene IDs or symbols (remove spaces)
    gene_ids <- strsplit(gsub("\\s+", "", input$user_input), ",")[[1]]
    
    # Conditional logic based on input type
    if (input$input_type == "outdated_ensembl") {
      # Use Ensembl version 80 to retrieve gene symbols
      ensembl_v80 <- useEnsembl(biomart = "ensembl", version = 80, dataset = "rnorvegicus_gene_ensembl")
      
      # Query for gene symbols based on outdated Ensembl IDs
      results_old <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                           filters = 'ensembl_gene_id', 
                           values = gene_ids, 
                           mart = ensembl_v80)
      
      return(results_old)
      
    } else if (input$input_type == "current_ensembl") {
      # Use the latest Ensembl version to map to new Ensembl IDs
      ensembl_latest <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
      
      # Query for gene symbols based on current Ensembl IDs
      results_new <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                           filters = 'ensembl_gene_id', 
                           values = gene_ids, 
                           mart = ensembl_latest)
      
      return(results_new)
      
    } else if (input$input_type == "gene_symbol") {
      # Use mygene.info to query for Ensembl IDs and aliases based on gene symbols
      gene_info_results <- do.call(rbind, lapply(gene_ids, function(gene_symbol) {
        url <- paste0("https://mygene.info/v3/query?q=", gene_symbol, 
                      "&fields=ensembl.gene,alias,symbol&species=rat")
        response <- GET(url)
        content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
        
        # Extract the Ensembl ID and alias if available
        ensembl_id <- if (!is.null(content$hits$ensembl$gene)) content$hits$ensembl$gene[1] else NA
        aliases <- if (!is.null(content$hits$alias)) paste(content$hits$alias, collapse = ", ") else NA
        
        return(data.frame(Gene_Symbol = gene_symbol, Ensembl_ID = ensembl_id, Aliases = aliases))
      }))
      
      return(gene_info_results)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)

