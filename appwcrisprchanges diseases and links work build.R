# Load necessary packages
load_packages <- function(pkgs) {
  sapply(pkgs, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  })
}

packages <- c("shiny", "writexl", "DT", "readxl", "Autoseed", "httr", "jsonlite", "dplyr")
load_packages(packages)

# Load data from the Autoseed package
data("drugbank")
data("edgar")
data("mala")

# Function to query disease-causing mutations (e.g., ClinVar)
get_disease_mutations <- function(gene) {
  # Example with shiny::tags$a() to open the link in a new tab in the default browser
  mutations <- paste(shiny::tags$a("Mutation 1", href = "https://www.ncbi.nlm.nih.gov/clinvar/variation/12345", target = "_blank"),
                     shiny::tags$a("Mutation 2", href = "https://www.ncbi.nlm.nih.gov/clinvar/variation/67890", target = "_blank"),
                     sep = ", ")
  return(mutations)
}

# Function to query gene expression in disease context (using GTEx or MalaCards)
get_gene_expression_disease <- function(gene) {
  # Placeholder for gene expression data
  expression <- "Upregulated in Disease A, Downregulated in Disease B"  # Example output
  return(expression)
}

# Function to query gene function and pathways (using KEGG, Reactome, or MalaCards)
get_gene_function_pathways <- function(gene) {
  pathways <- paste(shiny::tags$a("Pathway X", href = "https://www.kegg.jp/pathway/hsa04610", target = "_blank"),
                    shiny::tags$a("Pathway Y", href = "https://reactome.org/PathwayBrowser/#/R-HSA-2029480", target = "_blank"),
                    sep = ", ")
  return(pathways)
}

# Function to get gene-related information including mutations, expression, and pathways
get_gene_related_info <- function(gene, species) {
  # Filter the datasets by the provided gene
  drugbank_results <- subset(drugbank, gene == drugbank[, 1])
  edgar_results <- subset(edgar, gene == edgar[, 1])
  malacards_results <- subset(mala, gene == mala[, 2])
  
  # Combine results into a single data frame
  combined_results <- data.frame(
    Gene = gene,
    Disease = c(drugbank_results[, 2],
                edgar_results[, 2],
                malacards_results[, 1])
  )
  
  # Remove duplicate rows
  combined_results <- combined_results %>% distinct()
  
  # Get additional information for the gene
  combined_results$Mutations <- get_disease_mutations(gene)
  combined_results$Expression <- get_gene_expression_disease(gene)
  combined_results$Pathways <- get_gene_function_pathways(gene)
  
  return(combined_results)
}

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
  titlePanel("Gene-Disease Association Dashboard with Mutations, Expression, and Pathways"),
  
  # Tabs for species selection
  tabsetPanel(
    id = "species",
    tabPanel("Human", value = "human",
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_input_human", "Enter gene symbol:", ""),
                 textInput("file_name_human", "Enter file name for download:", "gene_disease_results")
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_human"),
                 downloadButton("download_data_human", "Download Full Results")
               )
             )
    ),
    tabPanel("Mouse", value = "mouse", 
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_input_mouse", "Enter gene symbol:", ""),
                 textInput("file_name_mouse", "Enter file name for download:", "gene_disease_results")
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_mouse"),
                 downloadButton("download_data_mouse", "Download Full Results")
               )
             )
    ),
    tabPanel("Rat", value = "rat", 
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_input_rat", "Enter gene symbol:", ""),
                 textInput("file_name_rat", "Enter file name for download:", "gene_disease_results")
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_rat"),
                 downloadButton("download_data_rat", "Download Full Results")
               )
             )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive values to store results for each species
  results_human <- reactiveVal(data.frame())
  results_mouse <- reactiveVal(data.frame())
  results_rat <- reactiveVal(data.frame())
  
  # Function to process input for each species and perform the query
  process_input <- function(gene_input, species) {
    if (gene_input == "") {
      return(data.frame(Gene = NA, Disease = NA, Mutations = NA, Expression = NA, Pathways = NA))
    }
    
    # Initialize progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Query in progress...", value = 0)
    on.exit(progress$close())  # Close progress bar when done
    
    # Query gene-related diseases, mutations, expression, and pathways
    gene_info_results <- get_gene_related_info(gene_input, species)
    
    progress$inc(1, detail = paste("Processing", gene_input))
    
    return(gene_info_results)
  }
  
  # Reactive outputs for each species
  output$result_table_human <- DT::renderDataTable({
    result_data <- process_input(input$gene_input_human, "human")
    results_human(result_data)
    DT::datatable(result_data, escape = FALSE)  # escape = FALSE to allow links
  })
  
  output$result_table_mouse <- DT::renderDataTable({
    result_data <- process_input(input$gene_input_mouse, "mouse")
    results_mouse(result_data)
    DT::datatable(result_data, escape = FALSE)
  })
  
  output$result_table_rat <- DT::renderDataTable({
    result_data <- process_input(input$gene_input_rat, "rat")
    results_rat(result_data)
    DT::datatable(result_data, escape = FALSE)
  })
  
  # Download handlers for each species
  output$download_data_human <- downloadHandler(
    filename = function() {
      paste(input$file_name_human, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(results_human(), file)
    }
  )
  
  output$download_data_mouse <- downloadHandler(
    filename = function() {
      paste(input$file_name_mouse, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(results_mouse(), file)
    }
  )
  
  output$download_data_rat <- downloadHandler(
    filename = function() {
      paste(input$file_name_rat, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(results_rat(), file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)


  
  