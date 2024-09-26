# Load necessary packages
if (!require("shiny")) install.packages("shiny")
if (!require("biomaRt")) install.packages("biomaRt")
if (!require("httr")) install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("writexl")) install.packages("writexl")
if (!require("DT")) install.packages("DT")
if (!require("readxl")) install.packages("readxl")

library(shiny)
library(biomaRt)
library(httr)
library(jsonlite)
library(writexl)
library(DT)
library(readxl)

# Function to query Ensembl REST API for a new Ensembl ID
get_new_ensembl_id <- function(gene_symbol, species) {
  species_map <- list(
    "human" = "homo_sapiens",
    "mouse" = "mus_musculus",
    "rat" = "rattus_norvegicus"
  )
  species_name <- species_map[[species]]
  base_url <- paste0("https://rest.ensembl.org/xrefs/symbol/", species_name, "/", gene_symbol, "?content-type=application/json")
  
  result <- tryCatch({
    response <- httr::GET(base_url)
    if (httr::status_code(response) == 200) {
      parsed_content <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))
      if (length(parsed_content) > 0 && !is.null(parsed_content$id)) {
        return(parsed_content$id[1])
      } else {
        return(NA)
      }
    } else {
      return(NA)
    }
  }, error = function(e) {
    return(NA)
  })
  
  return(result)
}

# Function to process outdated Ensembl IDs with version selection
process_outdated_ensembl <- function(ensembl_ids, species, version, progress) {
  
  # Ensure the version is selected properly
  if (is.null(version) || version == "") {
    stop("No Ensembl version selected.")
  }
  
  # Dataset selection based on species
  dataset <- switch(species, 
                    "human" = "hsapiens_gene_ensembl", 
                    "mouse" = "mmusculus_gene_ensembl", 
                    "rat" = "rnorvegicus_gene_ensembl")
  
  # Ensure the version is passed as an integer
  version <- as.integer(version)
  
  # Query the appropriate Ensembl version
  ensembl_v <- useEnsembl(biomart = "ensembl", version = version, dataset = dataset)
  
  # Progress message
  progress$set(message = paste0("Querying Ensembl v", version, " for ", species, "..."), value = 0)
  
  # Query for gene symbols using old Ensembl IDs
  results_old <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                       filters = 'ensembl_gene_id', 
                       values = ensembl_ids, 
                       mart = ensembl_v)
  progress$inc(0.3, detail = paste0("Ensembl IDs queried from version ", version))
  
  return(results_old)
}

# Function to get gene info using MyGene.info
get_gene_info <- function(gene_symbol, species) {
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
  
  # Return data frame with Ensembl ID and aliases
  return(data.frame(Gene_Symbol = gene_symbol, Ensembl_ID = ensembl_id, Aliases = aliases))
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
  titlePanel("Outdated and Current Ensembl ID Search with Species Selection"),
  
  # Tabs for species selection
  tabsetPanel(
    id = "species",
    
    # Human tab
    tabPanel("Human", value = "human",
             sidebarLayout(
               sidebarPanel(
                 # Dropdown for selecting the Ensembl version for Human
                 selectInput("ensembl_version_human", "Select Ensembl Version:",
                             choices = list("Ensembl v76" = 76, 
                                            "Ensembl v80" = 80, 
                                            "Ensembl v85" = 85, 
                                            "Ensembl v90" = 90, 
                                            "Ensembl v95" = 95, 
                                            "Ensembl v100" = 100, 
                                            "Ensembl v105" = 105, 
                                            "Ensembl v110" = 110),
                             selected = 100),
                 
                 # Text or file input for current Ensembl IDs
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
                 
                 # File input for outdated Ensembl IDs
                 fileInput("outdated_file_human", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                 textInput("outdated_file_name_human", "Enter file name for download:", "outdated_gene_query_results")
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_human"),
                 textOutput("percent_mapped_human"),
                 downloadButton("download_data_human", "Download Full Results"),
                 downloadButton("download_filtered_human", "Download Filtered Results"),
                 downloadButton("download_outdated_human", "Download Outdated Ensembl Results")
               )
             )
    ),
    
    # Mouse tab
    tabPanel("Mouse", value = "mouse",
             sidebarLayout(
               sidebarPanel(
                 # Dropdown for selecting the Ensembl version for Mouse
                 selectInput("ensembl_version_mouse", "Select Ensembl Version:",
                             choices = list("Ensembl v76" = 76, 
                                            "Ensembl v80" = 80, 
                                            "Ensembl v85" = 85, 
                                            "Ensembl v90" = 90, 
                                            "Ensembl v95" = 95, 
                                            "Ensembl v100" = 100, 
                                            "Ensembl v105" = 105, 
                                            "Ensembl v110" = 110),
                             selected = 100),
                 
                 # Text or file input for current Ensembl IDs
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
                 
                 # File input for outdated Ensembl IDs
                 fileInput("outdated_file_mouse", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                 textInput("outdated_file_name_mouse", "Enter file name for download:", "outdated_gene_query_results")
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_mouse"),
                 textOutput("percent_mapped_mouse"),
                 downloadButton("download_data_mouse", "Download Full Results"),
                 downloadButton("download_filtered_mouse", "Download Filtered Results"),
                 downloadButton("download_outdated_mouse", "Download Outdated Ensembl Results")
               )
             )
    ),
    
    # Rat tab
    tabPanel("Rat", value = "rat",
             sidebarLayout(
               sidebarPanel(
                 # Dropdown for selecting the Ensembl version for Rat
                 selectInput("ensembl_version_rat", "Select Ensembl Version:",
                             choices = list("Ensembl v76 (Rnor6.0)" = 76, 
                                            "Ensembl v80 (Rnor6.0)" = 80, 
                                            "Ensembl v85 (Rnor6.0)" = 85, 
                                            "Ensembl v90 (Rnor6.0)" = 90, 
                                            "Ensembl v95 (Rnor6.0)" = 95, 
                                            "Ensembl v100 (Mrat7.2)" = 100, 
                                            "Ensembl v101 (Mrat7.2)" = 101, 
                                            "Ensembl v104 (Mrat7.2)" = 104, 
                                            "Ensembl v110 (Mrat7.2)" = 110),
                             selected = 100),
                 
                 # Text or file input for current Ensembl IDs
                 radioButtons("input_type_rat", "Select Gene ID type:",
                              choices = list("Gene Symbol" = "gene_symbol", 
                                             "Current Ensembl ID" = "ensembl")),
                 radioButtons("input_source_rat", "Select input source:",
                              choices = list("Text Input" = "text", 
                                             "Excel Input" = "excel")),
                 conditionalPanel(
                   condition = "input.input_source_rat == 'text' && input.input_type_rat == 'gene_symbol'",
                   textInput("user_input_rat", "Enter comma-separated gene symbols:", "")
                 ),
                 conditionalPanel(
                   condition = "input.input_source_rat == 'text' && input.input_type_rat == 'ensembl'",
                   textInput("ensembl_input_rat", "Enter comma-separated Ensembl IDs:", "")
                 ),
                 conditionalPanel(
                   condition = "input.input_source_rat == 'excel'",
                   fileInput("file_input_rat", "Upload Excel file", accept = c(".xlsx"))
                 ),
                 
                 # File input for outdated Ensembl IDs
                 fileInput("outdated_file_rat", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                 textInput("outdated_file_name_rat", "Enter file name for download:", "outdated_gene_query_results")
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_rat"),
                 textOutput("percent_mapped_rat"),
                 downloadButton("download_data_rat", "Download Full Results"),
                 downloadButton("download_filtered_rat", "Download Filtered Results"),
                 downloadButton("download_outdated_rat", "Download Outdated Ensembl Results")
               )
             )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive values to store results for each species
  outdated_results_human <- reactiveVal(data.frame())
  outdated_results_mouse <- reactiveVal(data.frame())
  outdated_results_rat <- reactiveVal(data.frame())
  results_human <- reactiveVal(data.frame())
  results_mouse <- reactiveVal(data.frame())
  results_rat <- reactiveVal(data.frame())
  
  # Process outdated Ensembl IDs for Human
  observeEvent(input$outdated_file_human, {
    req(input$outdated_file_human)
    ensembl_ids <- read_excel(input$outdated_file_human$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])  # Assuming the Ensembl IDs are in the first column
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    version_selected <- input$ensembl_version_human
    result_data <- process_outdated_ensembl(ensembl_ids, "human", version_selected, progress)
    outdated_results_human(result_data)
    
    output$result_table_human <- DT::renderDataTable({
      DT::datatable(outdated_results_human())
    })
  })
  
  # Process current Ensembl IDs for Human
  observeEvent(input$user_input_human, {
    req(input$user_input_human)
    gene_ids <- strsplit(gsub("\\s+", "", input$user_input_human), ",")[[1]]
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing current Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    gene_info_results <- do.call(rbind, lapply(gene_ids, function(gene_symbol) {
      progress$inc(1 / length(gene_ids), detail = paste("Processing", gene_symbol))
      result <- get_gene_info(gene_symbol, "human")
      return(result)
    }))
    
    results_human(gene_info_results)
    
    output$result_table_human <- DT::renderDataTable({
      DT::datatable(results_human())
    })
  })
  
  # Process outdated Ensembl IDs for Mouse
  observeEvent(input$outdated_file_mouse, {
    req(input$outdated_file_mouse)
    ensembl_ids <- read_excel(input$outdated_file_mouse$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])  # Assuming the Ensembl IDs are in the first column
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    version_selected <- input$ensembl_version_mouse
    result_data <- process_outdated_ensembl(ensembl_ids, "mouse", version_selected, progress)
    outdated_results_mouse(result_data)
    
    output$result_table_mouse <- DT::renderDataTable({
      DT::datatable(outdated_results_mouse())
    })
  })
  
  # Process current Ensembl IDs for Mouse
  observeEvent(input$user_input_mouse, {
    req(input$user_input_mouse)
    gene_ids <- strsplit(gsub("\\s+", "", input$user_input_mouse), ",")[[1]]
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing current Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    gene_info_results <- do.call(rbind, lapply(gene_ids, function(gene_symbol) {
      progress$inc(1 / length(gene_ids), detail = paste("Processing", gene_symbol))
      result <- get_gene_info(gene_symbol, "mouse")
      return(result)
    }))
    
    results_mouse(gene_info_results)
    
    output$result_table_mouse <- DT::renderDataTable({
      DT::datatable(results_mouse())
    })
  })
  
  # Process outdated Ensembl IDs for Rat
  observeEvent(input$outdated_file_rat, {
    req(input$outdated_file_rat)
    ensembl_ids <- read_excel(input$outdated_file_rat$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])  # Assuming the Ensembl IDs are in the first column
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    version_selected <- input$ensembl_version_rat
    result_data <- process_outdated_ensembl(ensembl_ids, "rat", version_selected, progress)
    outdated_results_rat(result_data)
    
    output$result_table_rat <- DT::renderDataTable({
      DT::datatable(outdated_results_rat())
    })
  })
  
  # Process current Ensembl IDs for Rat
  observeEvent(input$user_input_rat, {
    req(input$user_input_rat)
    gene_ids <- strsplit(gsub("\\s+", "", input$user_input_rat), ",")[[1]]
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing current Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    gene_info_results <- do.call(rbind, lapply(gene_ids, function(gene_symbol) {
      progress$inc(1 / length(gene_ids), detail = paste("Processing", gene_symbol))
      result <- get_gene_info(gene_symbol, "rat")
      return(result)
    }))
    
    results_rat(gene_info_results)
    
    output$result_table_rat <- DT::renderDataTable({
      DT::datatable(results_rat())
    })
  })
  
  # Download handlers for each species
  output$download_outdated_human <- downloadHandler(
    filename = function() {
      paste(input$outdated_file_name_human, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(outdated_results_human(), file)
    }
  )
  
  output$download_outdated_mouse <- downloadHandler(
    filename = function() {
      paste(input$outdated_file_name_mouse, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(outdated_results_mouse(), file)
    }
  )
  
  output$download_outdated_rat <- downloadHandler(
    filename = function() {
      paste(input$outdated_file_name_rat, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(outdated_results_rat(), file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
