# Load necessary packages
if (!require("shiny")) install.packages("shiny")
if (!require("biomaRt")) install.packages("biomaRt")
if (!require("httr")) install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("writexl")) install.packages("writexl")
if (!require("DT")) install.packages("DT")
if (!require("KEGGREST")) install.packages("KEGGREST")
if (!require("readxl")) install.packages("readxl")

library(shiny)
library(biomaRt)
library(httr)
library(jsonlite)
library(writexl)
library(DT)
library(KEGGREST)
library(readxl)

# Function to query KEGG pathways using KEGGREST and limit to top 5 pathways
get_kegg_pathways <- function(gene_symbol, species) {
  kegg_species_code <- switch(species, 
                              "human" = "hsa", 
                              "mouse" = "mmu", 
                              "rat" = "rno")
  
  kegg_info <- tryCatch({
    kegg_res <- keggGet(paste0(kegg_species_code, ":", gene_symbol))
    if (length(kegg_res) > 0 && !is.null(kegg_res[[1]]$PATHWAY)) {
      pathways <- paste(head(kegg_res[[1]]$PATHWAY, 5), collapse = ", ")  # Limit to top 5
      return(pathways)
    } else {
      return(NA)
    }
  }, error = function(e) {
    return(NA)  # Return NA if an error occurs
  })
  
  return(kegg_info)
}

# Function to limit GO terms based on keyword and prioritize them
limit_go_terms <- function(go_terms, keyword) {
  go_terms <- unlist(go_terms)  # Flatten the list of terms into a vector
  go_terms <- unique(go_terms)  # Remove duplicates
  
  # Prioritize terms with the keyword
  keyword_matches <- grep(keyword, go_terms, ignore.case = TRUE, value = TRUE)
  other_terms <- setdiff(go_terms, keyword_matches)
  
  prioritized_terms <- c(keyword_matches, other_terms)
  
  # Limit to top 5 terms
  return(paste(head(prioritized_terms, 5), collapse = ", "))
}

# Function to query MyGene.info API (adjusted for correct species names)
get_gene_info <- function(gene_symbol, keyword, species) {
  species_name <- switch(species, 
                         "human" = "human", 
                         "mouse" = "mouse", 
                         "rat" = "rat")
  
  url <- paste0("https://mygene.info/v3/query?q=", gene_symbol, "&fields=ensembl.gene,alias,symbol,go.BP,KEGG&species=", species_name)
  response <- GET(url)
  content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
  
  # Extract Ensembl ID
  ensembl_id <- if (!is.null(content$hits$ensembl$gene)) content$hits$ensembl$gene[1] else NA
  
  # Flatten Aliases (remove NULLs and keep all aliases as a string)
  aliases <- if (!is.null(content$hits$alias) && length(content$hits$alias) > 0) {
    paste(unique(unlist(content$hits$alias)), collapse = ", ")
  } else NA
  
  # Limit top 5 GO terms (Biological Process) with keyword prioritization
  go_terms <- if (!is.null(content$hits$go$BP)) {
    terms <- sapply(content$hits$go$BP, `[[`, "term")
    limit_go_terms(terms, keyword)  # Prioritize based on keyword
  } else NA
  
  # Query KEGG pathways
  kegg_pathways <- get_kegg_pathways(gene_symbol, species)
  
  return(data.frame(Gene_Symbol = gene_symbol, Ensembl_ID = ensembl_id, Aliases = aliases, GO_Terms = go_terms, KEGG_Pathways = kegg_pathways))
}

# Function to get gene symbols from current Ensembl IDs using biomaRt
get_ensembl_to_gene_symbol <- function(ensembl_ids, species) {
  dataset <- switch(species, 
                    "human" = "hsapiens_gene_ensembl", 
                    "mouse" = "mmusculus_gene_ensembl", 
                    "rat" = "rnorvegicus_gene_ensembl")
  
  ensembl <- useMart("ensembl", dataset = dataset)
  
  results <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                   filters = 'ensembl_gene_id', 
                   values = ensembl_ids, 
                   mart = ensembl)
  
  # Ensure we keep all input Ensembl IDs
  full_results <- merge(data.frame(Ensembl_ID = ensembl_ids), results, by.x = "Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)
  
  return(full_results)
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
  titlePanel("Gene Symbol and Ensembl ID Search with Species Selection"),
  
  # Tabs for species selection
  tabsetPanel(
    id = "species",
    tabPanel("Human", value = "human",
             sidebarLayout(
               sidebarPanel(
                 radioButtons("input_type_human", "Select Gene ID type:",
                              choices = list("Gene Symbol" = "gene_symbol", 
                                             "Current Ensembl ID" = "ensembl")),
                 radioButtons("input_source_human", "Select input source:",
                              choices = list("Text Input" = "text", 
                                             "Excel Input" = "excel")),
                 conditionalPanel(
                   condition = "input.input_source_human == 'text' && input.input_type_human == 'gene_symbol'",
                   textInput("user_input_human", "Enter comma-separated gene symbols:", "")
                 ),
                 conditionalPanel(
                   condition = "input.input_source_human == 'text' && input.input_type_human == 'ensembl'",
                   textInput("ensembl_input_human", "Enter comma-separated Ensembl IDs:", "")
                 ),
                 conditionalPanel(
                   condition = "input.input_source_human == 'excel'",
                   fileInput("file_input_human", "Upload Excel file", accept = c(".xlsx"))
                 ),
                 textInput("file_name_human", "Enter file name for download:", "gene_query_results"),
                 textInput("keyword_human", "Keyword Search (GO Only):", ""),
                 helpText("Enter gene symbols or Ensembl IDs separated by commas, or upload an Excel file.")
               ),
               mainPanel(
                 uiOutput("progress_bar_human"),
                 DT::dataTableOutput("result_table_human"),
                 downloadButton("download_data_human", "Download Full Results"),
                 downloadButton("download_filtered_human", "Download Filtered Results"),
                 h3("Gene Information"),
                 verbatimTextOutput("gene_info_human")
               )
             )
    ),
    tabPanel("Mouse", value = "mouse", 
             sidebarLayout(
               sidebarPanel(
                 radioButtons("input_type_mouse", "Select Gene ID type:",
                              choices = list("Gene Symbol" = "gene_symbol", 
                                             "Current Ensembl ID" = "ensembl")),
                 radioButtons("input_source_mouse", "Select input source:",
                              choices = list("Text Input" = "text", 
                                             "Excel Input" = "excel")),
                 conditionalPanel(
                   condition = "input.input_source_mouse == 'text' && input.input_type_mouse == 'gene_symbol'",
                   textInput("user_input_mouse", "Enter comma-separated gene symbols:", "")
                 ),
                 conditionalPanel(
                   condition = "input.input_source_mouse == 'text' && input.input_type_mouse == 'ensembl'",
                   textInput("ensembl_input_mouse", "Enter comma-separated Ensembl IDs:", "")
                 ),
                 conditionalPanel(
                   condition = "input.input_source_mouse == 'excel'",
                   fileInput("file_input_mouse", "Upload Excel file", accept = c(".xlsx"))
                 ),
                 textInput("file_name_mouse", "Enter file name for download:", "gene_query_results"),
                 textInput("keyword_mouse", "Keyword Search (GO Only):", ""),
                 helpText("Enter gene symbols or Ensembl IDs separated by commas, or upload an Excel file.")
               ),
               mainPanel(
                 uiOutput("progress_bar_mouse"),
                 DT::dataTableOutput("result_table_mouse"),
                 downloadButton("download_data_mouse", "Download Full Results"),
                 downloadButton("download_filtered_mouse", "Download Filtered Results"),
                 h3("Gene Information"),
                 verbatimTextOutput("gene_info_mouse")
               )
             )
    ),
    tabPanel("Rat", value = "rat", 
             sidebarLayout(
               sidebarPanel(
                 radioButtons("input_type_rat", "Select Gene ID type:",
                              choices = list("Gene Symbol" = "gene_symbol", 
                                             "Current Ensembl ID" = "ensembl")),
                 radioButtons("input_source_rat", "Select input source:",
                              choices = list("Text Input" = "text", 
                                             "Excel Input" = "excel")),
                 
                 # Text Input section
                 conditionalPanel(
                   condition = "input.input_source_rat == 'text'",
                   conditionalPanel(
                     condition = "input.input_type_rat == 'gene_symbol'",
                     textInput("user_input_rat", "Enter comma-separated gene symbols:", "")
                   ),
                   conditionalPanel(
                     condition = "input.input_type_rat == 'ensembl'",
                     textInput("ensembl_input_rat", "Enter comma-separated Ensembl IDs:", "")
                   ),
                   textInput("file_name_rat", "Enter file name for download:", "gene_query_results"),
                   textInput("keyword_rat", "Keyword Search (GO Only):", "")
                 ),
                 
                 # Excel Input section
                 conditionalPanel(
                   condition = "input.input_source_rat == 'excel'",
                   fileInput("file_input_rat", "Upload Excel file", accept = c(".xlsx")),
                   textInput("file_name_rat", "Enter file name for download:", "gene_query_results"),
                   textInput("keyword_rat", "Keyword Search (GO Only):", "")
                 ),
                 
                 helpText("Enter gene symbols or Ensembl IDs separated by commas, or upload an Excel file.")
               ),
               mainPanel(
                 uiOutput("progress_bar_rat"),
                 DT::dataTableOutput("result_table_rat"),
                 downloadButton("download_data_rat", "Download Full Results"),
                 downloadButton("download_filtered_rat", "Download Filtered Results"),
                 h3("Gene Information"),
                 verbatimTextOutput("gene_info_rat")
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
  process_input <- function(species, input_type, input_source, user_input, ensembl_input, file_input, keyword) {
    if (input_source == "text") {
      if (input_type == "gene_symbol") {
        gene_ids <- strsplit(gsub("\\s+", "", user_input), ",")[[1]]
      } else if (input_type == "ensembl") {
        gene_ids <- strsplit(gsub("\\s+", "", ensembl_input), ",")[[1]]
      }
    } else if (input_source == "excel") {
      req(file_input)
      file_data <- read_excel(file_input$datapath)
      gene_ids <- as.vector(file_data[[1]])
    }
    
    # Initialize progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Query in progress...", value = 0)
    on.exit(progress$close())  # Close progress bar when done
    
    if (input_type == "gene_symbol") {
      gene_info_results <- do.call(rbind, lapply(gene_ids, function(gene_symbol) {
        progress$inc(1 / length(gene_ids), detail = paste("Processing", gene_symbol))
        result <- get_gene_info(gene_symbol, keyword, species)
        return(result)
      }))
    } else if (input_type == "ensembl") {
      ensembl_to_gene <- get_ensembl_to_gene_symbol(gene_ids, species)
      gene_info_results <- do.call(rbind, lapply(ensembl_to_gene$external_gene_name, function(gene_symbol) {
        progress$inc(1 / length(ensembl_to_gene$external_gene_name), detail = paste("Processing", gene_symbol))
        result <- get_gene_info(gene_symbol, keyword, species)
        return(result)
      }))
    }
    
    return(gene_info_results)
  }
  
  # Reactive values to handle different species results
  output$result_table_human <- DT::renderDataTable({
    result_data <- process_input("human", input$input_type_human, input$input_source_human, 
                                 input$user_input_human, input$ensembl_input_human, 
                                 input$file_input_human, input$keyword_human)
    results_human(result_data)
    DT::datatable(result_data)
  })
  
  output$result_table_mouse <- DT::renderDataTable({
    result_data <- process_input("mouse", input$input_type_mouse, input$input_source_mouse, 
                                 input$user_input_mouse, input$ensembl_input_mouse, 
                                 input$file_input_mouse, input$keyword_mouse)
    results_mouse(result_data)
    DT::datatable(result_data)
  })
  
  output$result_table_rat <- DT::renderDataTable({
    result_data <- process_input("rat", input$input_type_rat, input$input_source_rat, 
                                 input$user_input_rat, input$ensembl_input_rat, 
                                 input$file_input_rat, input$keyword_rat)
    results_rat(result_data)
    DT::datatable(result_data)
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
  
  # Download handlers for filtered results
  output$download_filtered_human <- downloadHandler(
    filename = function() {
      paste(input$file_name_human, "_filtered.xlsx", sep = "")
    },
    content = function(file) {
      filtered_data <- results_human()[input$result_table_human_rows_all, ]
      write_xlsx(filtered_data, file)
    }
  )
  
  output$download_filtered_mouse <- downloadHandler(
    filename = function() {
      paste(input$file_name_mouse, "_filtered.xlsx", sep = "")
    },
    content = function(file) {
      filtered_data <- results_mouse()[input$result_table_mouse_rows_all, ]
      write_xlsx(filtered_data, file)
    }
  )
  
  output$download_filtered_rat <- downloadHandler(
    filename = function() {
      paste(input$file_name_rat, "_filtered.xlsx", sep = "")
    },
    content = function(file) {
      filtered_data <- results_rat()[input$result_table_rat_rows_all, ]
      write_xlsx(filtered_data, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)


                                             