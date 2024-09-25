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

# Function to query MyGene.info API
get_gene_info <- function(gene_symbol, species) {
  species_name <- switch(species, 
                         "human" = "human", 
                         "mouse" = "mouse", 
                         "rat" = "rat")
  
  url <- paste0("https://mygene.info/v3/query?q=", gene_symbol, "&fields=ensembl.gene&species=", species_name)
  response <- httr::GET(url)
  content <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))
  
  # Extract Ensembl ID
  ensembl_id <- if (!is.null(content$hits$ensembl$gene)) content$hits$ensembl$gene[1] else NA
  
  return(data.frame(Ensembl_ID = ensembl_id))
}

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

# Advanced process_outdated_ensembl function
process_outdated_ensembl <- function(ensembl_ids, species, progress) {
  # Print ensembl_ids to troubleshoot
  print("Ensembl IDs:")
  print(ensembl_ids)
  
  # Step 1: Query old Ensembl IDs from version 80
  ensembl_v80 <- useEnsembl(biomart = "ensembl", version = 80, dataset = switch(species, 
                                                                                "human" = "hsapiens_gene_ensembl", 
                                                                                "mouse" = "mmusculus_gene_ensembl", 
                                                                                "rat" = "rnorvegicus_gene_ensembl"))
  progress$set(message = "Querying Ensembl version 80...", value = 0)
  
  # Query for gene symbols using old Ensembl IDs
  results_old <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                       filters = 'ensembl_gene_id', 
                       values = ensembl_ids, 
                       mart = ensembl_v80)
  progress$inc(0.3, detail = "Old Ensembl IDs queried")
  print("Results from Ensembl v80:")
  print(results_old)
  
  # Step 2: Query the latest Ensembl IDs from biomaRt using gene symbols
  ensembl_latest <- useMart("ensembl", dataset = switch(species, 
                                                        "human" = "hsapiens_gene_ensembl", 
                                                        "mouse" = "mmusculus_gene_ensembl", 
                                                        "rat" = "rnorvegicus_gene_ensembl"))
  gene_symbols <- results_old$external_gene_name  # Extract gene symbols
  print("Gene symbols from v80:")
  print(gene_symbols)
  
  # Query for new Ensembl IDs using gene symbols from biomaRt
  results_new_biomart <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                               filters = 'external_gene_name', 
                               values = gene_symbols, 
                               mart = ensembl_latest)
  progress$inc(0.3, detail = "New Ensembl IDs from biomaRt queried")
  print("Results from biomaRt:")
  print(results_new_biomart)
  
  # Step 3: Query MyGene.info for new Ensembl IDs using gene symbols
  results_new_mygene <- do.call(rbind, lapply(gene_symbols, function(symbol) {
    result <- get_gene_info(symbol, species)
    return(data.frame(Gene_Symbol = symbol, New_Ensembl_ID_MyGene = result$Ensembl_ID))
  }))
  progress$inc(0.3, detail = "New Ensembl IDs from MyGene.info queried")
  print("Results from MyGene.info:")
  print(results_new_mygene)
  
  # Step 4: Query Ensembl REST API for new Ensembl IDs using gene symbols
  results_new_ensembl_rest <- do.call(rbind, lapply(gene_symbols, function(symbol) {
    new_ensembl_id <- get_new_ensembl_id(symbol, species)
    return(data.frame(Gene_Symbol = symbol, New_Ensembl_ID_EnsemblRest = new_ensembl_id))
  }))
  progress$inc(0.3, detail = "New Ensembl IDs from Ensembl REST API queried")
  print("Results from Ensembl REST API:")
  print(results_new_ensembl_rest)
  
  # Step 5: Merge biomaRt, MyGene.info, and Ensembl REST results with the old Ensembl data
  merged_results <- merge(results_old, results_new_biomart, by.x = "external_gene_name", by.y = "external_gene_name", all.x = TRUE)
  merged_results <- merge(merged_results, results_new_mygene, by.x = "external_gene_name", by.y = "Gene_Symbol", all.x = TRUE)
  merged_results <- merge(merged_results, results_new_ensembl_rest, by.x = "external_gene_name", by.y = "Gene_Symbol", all.x = TRUE)
  print("Merged Results:")
  print(merged_results)
  
  # Step 6: Rename columns for clarity
  colnames(merged_results) <- c("Gene_Symbol (biomaRt)", "Input (Ensembl v80)", 
                                "New_Ensembl_ID (biomaRt)", "New_Ensembl_ID (MyGene.info)", 
                                "New_Ensembl_ID (Ensembl REST)")
  
  return(merged_results)
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
                 
                 # New Advanced Dropdown Button
                 selectInput("advanced_options", "Advanced Options:",
                             choices = list("None" = "none", "Outdated Ensembl Lookup" = "test1"),
                             selected = "none"),
                 
                 # Conditional panel for Test 1 (Outdated Ensembl ID analysis)
                 conditionalPanel(
                   condition = "input.advanced_options == 'test1'",
                   fileInput("outdated_file_human", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                   textInput("outdated_file_name_human", "Enter file name for download:", "outdated_gene_query_results")
                 ),
                 
                 helpText("Enter gene symbols or Ensembl IDs separated by commas, or upload an Excel file.")
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_human"),  # For Human
                 textOutput("percent_mapped_human"),         # Add this to show the percentage for Human
                 downloadButton("download_data_human", "Download Full Results"),
                 downloadButton("download_filtered_human", "Download Filtered Results"),
                 
                 # Download button for outdated ensembl only visible if Test 1 selected
                 conditionalPanel(
                   condition = "input.advanced_options == 'test1'",
                   downloadButton("download_outdated_human", "Download Outdated Ensembl Results")
                 )
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
                 
                 # New Advanced Dropdown Button for Mouse
                 selectInput("advanced_options_mouse", "Advanced Options:",
                             choices = list("None" = "none", "Outdated Ensembl Lookup" = "test1"),
                             selected = "none"),
                 
                 # Conditional panel for Outdated Ensembl Lookup
                 conditionalPanel(
                   condition = "input.advanced_options_mouse == 'test1'",
                   fileInput("outdated_file_mouse", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                   textInput("outdated_file_name_mouse", "Enter file name for download:", "outdated_gene_query_results")
                 )
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_mouse"),  # For Mouse
                 textOutput("percent_mapped_mouse"),         # Add this to show the percentage for Mouse
                 downloadButton("download_data_mouse", "Download Full Results"),
                 downloadButton("download_filtered_mouse", "Download Filtered Results"),
                 
                 # Download button for outdated ensembl only visible if Test 1 selected
                 conditionalPanel(
                   condition = "input.advanced_options_mouse == 'test1'",
                   downloadButton("download_outdated_mouse", "Download Outdated Ensembl Results")
                 )
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
                 textInput("file_name_rat", "Enter file name for download:", "gene_query_results"),
                 textInput("keyword_rat", "Keyword Search (GO Only):", ""),
                 
                 # New Advanced Dropdown Button for Rat
                 selectInput("advanced_options_rat", "Advanced Options:",
                             choices = list("None" = "none", "Outdated Ensembl Lookup" = "test1"),
                             selected = "none"),
                 
                 # Conditional panel for Outdated Ensembl Lookup
                 conditionalPanel(
                   condition = "input.advanced_options_rat == 'test1'",
                   fileInput("outdated_file_rat", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                   textInput("outdated_file_name_rat", "Enter file name for download:", "outdated_gene_query_results")
                 )
               ),
               mainPanel(
                 DT::dataTableOutput("result_table_rat"),    # For Rat
                 textOutput("percent_mapped_rat"),           # Add this to show the percentage for Rat
                 downloadButton("download_data_rat", "Download Full Results"),
                 downloadButton("download_filtered_rat", "Download Filtered Results"),
                 
                 # Download button for outdated ensembl only visible if Test 1 selected
                 conditionalPanel(
                   condition = "input.advanced_options_rat == 'test1'",
                   downloadButton("download_outdated_rat", "Download Outdated Ensembl Results")
                 )
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
  outdated_results_human <- reactiveVal(data.frame())
  outdated_results_mouse <- reactiveVal(data.frame())
  outdated_results_rat <- reactiveVal(data.frame())
  
  # Process outdated Ensembl IDs from uploaded file for Human
  observeEvent(input$outdated_file_human, {
    req(input$outdated_file_human)
    ensembl_ids <- read_excel(input$outdated_file_human$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])  # Assuming the Ensembl IDs are in the first column
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    # Process the outdated Ensembl IDs
    result_data <- process_outdated_ensembl(ensembl_ids, "human", progress)
    outdated_results_human(result_data)  # Set the reactive value
    
    # Display the final merged results in the Shiny app
    output$result_table_human <- DT::renderDataTable({
      DT::datatable(result_data)
    })
    
    # Calculate the combined percentage directly from displayed table
    output$percent_mapped_human <- renderText({
      table_data <- outdated_results_human()
      total_old_ids <- nrow(table_data)
      
      # Count genes successfully mapped to any database
      mapped_to_any <- sum(!is.na(table_data$`New_Ensembl_ID (biomaRt)`) | 
                             !is.na(table_data$`New_Ensembl_ID (MyGene.info)`) | 
                             !is.na(table_data$`New_Ensembl_ID (Ensembl REST)`))
      
      # Calculate combined percentage
      if (total_old_ids > 0) {
        percent_mapped_to_any <- (mapped_to_any / total_old_ids) * 100
      } else {
        percent_mapped_to_any <- 0
      }
      
      paste0("Percentage of genes successfully mapped to any database: ", round(percent_mapped_to_any, 2), "%")
    })
  })
  
  # Download handler for outdated Ensembl results for Human
  output$download_outdated_human <- downloadHandler(
    filename = function() {
      paste(input$outdated_file_name_human, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(outdated_results_human(), file)
    }
  )
  
  # Process outdated Ensembl IDs from uploaded file for Mouse
  observeEvent(input$outdated_file_mouse, {
    req(input$outdated_file_mouse)
    ensembl_ids <- read_excel(input$outdated_file_mouse$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])  # Assuming the Ensembl IDs are in the first column
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    # Process the outdated Ensembl IDs
    result_data <- process_outdated_ensembl(ensembl_ids, "mouse", progress)
    outdated_results_mouse(result_data)  # Set the reactive value
    
    # Display the final merged results in the Shiny app
    output$result_table_mouse <- DT::renderDataTable({
      DT::datatable(result_data)
    })
    
    # Calculate the combined percentage directly from displayed table
    output$percent_mapped_mouse <- renderText({
      table_data <- outdated_results_mouse()
      total_old_ids <- nrow(table_data)
      
      # Count genes successfully mapped to any database
      mapped_to_any <- sum(!is.na(table_data$`New_Ensembl_ID (biomaRt)`) | 
                             !is.na(table_data$`New_Ensembl_ID (MyGene.info)`) | 
                             !is.na(table_data$`New_Ensembl_ID (Ensembl REST)`))
      
      # Calculate combined percentage
      if (total_old_ids > 0) {
        percent_mapped_to_any <- (mapped_to_any / total_old_ids) * 100
      } else {
        percent_mapped_to_any <- 0
      }
      
      paste0("Percentage of genes successfully mapped to any database: ", round(percent_mapped_to_any, 2), "%")
    })
  })
  
  # Download handler for outdated Ensembl results for Mouse
  output$download_outdated_mouse <- downloadHandler(
    filename = function() {
      paste(input$outdated_file_name_mouse, ".xlsx", sep = "")
    },
    content = function(file) {
      write_xlsx(outdated_results_mouse(), file)
    }
  )
  
  # Process outdated Ensembl IDs from uploaded file for Rat
  observeEvent(input$outdated_file_rat, {
    req(input$outdated_file_rat)
    ensembl_ids <- read_excel(input$outdated_file_rat$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])  # Assuming the Ensembl IDs are in the first column
    
    # Initialize the progress bar
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    
    # Close progress bar when done
    on.exit(progress$close())
    
    # Process the outdated Ensembl IDs
    result_data <- process_outdated_ensembl(ensembl_ids, "rat", progress)
    outdated_results_rat(result_data)  # Set the reactive value
    
    # Display the final merged results in the Shiny app
    output$result_table_rat <- DT::renderDataTable({
      DT::datatable(result_data)
    })
    
    # Calculate the combined percentage directly from displayed table
    output$percent_mapped_rat <- renderText({
      table_data <- outdated_results_rat()
      total_old_ids <- nrow(table_data)
      
      # Count genes successfully mapped to any database
      mapped_to_any <- sum(!is.na(table_data$`New_Ensembl_ID (biomaRt)`) | 
                             !is.na(table_data$`New_Ensembl_ID (MyGene.info)`) | 
                             !is.na(table_data$`New_Ensembl_ID (Ensembl REST)`))
      
      # Calculate combined percentage
      if (total_old_ids > 0) {
        percent_mapped_to_any <- (mapped_to_any / total_old_ids) * 100
      } else {
        percent_mapped_to_any <- 0
      }
      
      paste0("Percentage of genes successfully mapped to any database: ", round(percent_mapped_to_any, 2), "%")
    })
  })
  
  # Download handler for outdated Ensembl results for Rat
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

