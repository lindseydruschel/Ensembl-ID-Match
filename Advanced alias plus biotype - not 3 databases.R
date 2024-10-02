#Databases Queried: 
#biomaRt: Databases Queried - Ensembl Gene Database, Ensembl Variation Database, Ensembl Regulation Database, Gencode (Human, Mouse). 
  #Key Information Provided - Gene IDs, symbols, biotypes, genomic variants, regulatory elements.
#MyGene.info: Databases Queried - NCBI Gene, Ensembl, UCSC Genome Browser, RefSeq, UniProt. 
  #Key Information Provided - Gene annotations, IDs, protein information, functional annotations.
#Ensembl REST API: Databases Queried - Ensembl Gene Database, Ensembl Regulation Database, Ensembl Variation Database, Ensembl Protein Databases. 
  #Key Information Provided - Gene annotations, biotypes, transcript structures, regulatory elements.


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

# Function to query MyGene.info API
get_gene_info <- function(gene_symbol, species) {
  species_name <- switch(species, 
                         "human" = "human", 
                         "mouse" = "mouse", 
                         "rat" = "rat")
  
  url <- paste0("https://mygene.info/v3/query?q=", gene_symbol, "&fields=ensembl.gene,alias,symbol&species=", species_name)
  response <- httr::GET(url)
  content <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))
  
  # Extract Ensembl ID
  ensembl_id <- if (!is.null(content$hits$ensembl$gene)) content$hits$ensembl$gene[1] else NA
  
  # Extract and flatten aliases
  aliases <- if (!is.null(content$hits$alias) && length(content$hits$alias) > 0) {
    paste(unique(unlist(content$hits$alias)), collapse = ", ")
  } else NA
  
  return(data.frame(Ensembl_ID = ensembl_id, Aliases = aliases))
}


# Function to query Ensembl REST API for a new Ensembl ID and biotype
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
        ensembl_id <- parsed_content$id[1]
        
        # Query for biotype using Ensembl ID
        ensembl_details_url <- paste0("https://rest.ensembl.org/lookup/id/", ensembl_id, "?content-type=application/json")
        details_response <- httr::GET(ensembl_details_url)
        details_content <- jsonlite::fromJSON(httr::content(details_response, as = "text", encoding = "UTF-8"))
        
        # Extract biotype
        biotype <- if (!is.null(details_content$biotype)) details_content$biotype else "Unknown"
        
        return(list(Ensembl_ID = ensembl_id, Biotype = biotype))
      } else {
        return(list(Ensembl_ID = NA, Biotype = "Unknown"))
      }
    } else {
      return(list(Ensembl_ID = NA, Biotype = "Unknown"))
    }
  }, error = function(e) {
    return(list(Ensembl_ID = NA, Biotype = "Unknown"))
  })
  
  return(result)
}




# Function to process outdated Ensembl IDs with version and 3 database queries
process_outdated_ensembl <- function(ensembl_ids, species, version, progress = NULL) {
  ensembl_v <- useEnsembl(biomart = "ensembl", version = version, dataset = switch(species, 
                                                                                   "human" = "hsapiens_gene_ensembl", 
                                                                                   "mouse" = "mmusculus_gene_ensembl", 
                                                                                   "rat" = "rnorvegicus_gene_ensembl"))
  
  # Step 1: Query for gene symbols using old Ensembl IDs
  results_old <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                       filters = 'ensembl_gene_id', 
                       values = ensembl_ids, 
                       mart = ensembl_v)
  
  # Step 2: Query the latest Ensembl IDs from biomaRt using gene symbols
  ensembl_latest <- useMart("ensembl", dataset = switch(species, 
                                                        "human" = "hsapiens_gene_ensembl", 
                                                        "mouse" = "mmusculus_gene_ensembl", 
                                                        "rat" = "rnorvegicus_gene_ensembl"))
  gene_symbols <- results_old$external_gene_name  # Extract gene symbols
  
  results_new_biomart <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                               filters = 'external_gene_name', 
                               values = gene_symbols, 
                               mart = ensembl_latest)
  
  # Step 3: Query MyGene.info for new Ensembl IDs using gene symbols
  results_new_mygene <- do.call(rbind, lapply(gene_symbols, function(symbol) {
    result <- get_gene_info(symbol, species)  # Pass species argument here
    return(data.frame(Gene_Symbol = symbol, New_Ensembl_ID_MyGene = result$Ensembl_ID, Aliases = result$Aliases))
  }))
  
  
  # Step 4: Query Ensembl REST API for new Ensembl IDs using gene symbols
  results_new_ensembl_rest <- do.call(rbind, lapply(gene_symbols, function(symbol) {
    new_ensembl_id <- get_new_ensembl_id(symbol, species)
    
    # Return gene symbol, new Ensembl ID, and biotype
    return(data.frame(Gene_Symbol = symbol, 
                      New_Ensembl_ID_EnsemblRest = new_ensembl_id$Ensembl_ID, 
                      Biotype = new_ensembl_id$Biotype))
  }))
  
  # Step 6: Merge biomaRt, MyGene.info, Ensembl REST (with biotype), and other results
  merged_results <- merge(results_old, results_new_biomart, by.x = "external_gene_name", by.y = "external_gene_name", all.x = TRUE)
  merged_results <- merge(merged_results, results_new_mygene, by.x = "external_gene_name", by.y = "Gene_Symbol", all.x = TRUE)
  merged_results <- merge(merged_results, results_new_ensembl_rest, by.x = "external_gene_name", by.y = "Gene_Symbol", all.x = TRUE)
  
  return(merged_results)
}




# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background-color: #ffe6f0;
      }
    "))
  ),
  
  titlePanel("Outdated Ensembl ID Search"),
  
  tabsetPanel(
    id = "species",
    
    # Human Tab
    tabPanel("Human", value = "human",
             sidebarLayout(
               sidebarPanel(
                 selectInput("ensembl_version_human", "Select Ensembl Version:",
                             choices = list(
                              
                               "Ensembl v77 (Oct 2014, GRCh38)" = 77,
                               "Ensembl v80 (May 2015, GRCh38)" = 80,
                               "Ensembl v97 (Jul 2019, GRCh38)" = 97,
                               "Ensembl v100 (Apr 2020, GRCh38)" = 100,
                               "Ensembl v105 (Dec 2021, GRCh38)" = 105,
                               "Ensembl v110 (Jul 2023, GRCh38)" = 110
                             ),
                             selected = 100),
                 
                 fileInput("outdated_file_human", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                 textInput("outdated_file_name_human", "Enter file name for download:", "outdated_gene_query_results")
               ),
               mainPanel(
                 DT::dataTableOutput("outdated_result_table_human"),
                 textOutput("percent_mapped_human"),
                 downloadButton("download_outdated_human", "Download Outdated Ensembl Results")
               )
             )
    ),
    
    # Mouse Tab
    tabPanel("Mouse", value = "mouse",
             sidebarLayout(
               sidebarPanel(
                 selectInput("ensembl_version_mouse", "Select Ensembl Version:",
                             choices = list(
                               
                               "Ensembl v77 (Oct 2014, GRCm38)" = 77,
                               "Ensembl v80 (May 2015, GRCm38)" = 80,
                               "Ensembl v97 (Jul 2019, GRCm38)" = 97,
                               "Ensembl v100 (Apr 2020, GRCm38)" = 100,
                               "Ensembl v105 (Dec 2021, GRCm38)" = 105,
                               "Ensembl v110 (Jul 2023, GRCm39)" = 110
                             ),
                             selected = 100),
                 
                 fileInput("outdated_file_mouse", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                 textInput("outdated_file_name_mouse", "Enter file name for download:", "outdated_gene_query_results")
               ),
               mainPanel(
                 DT::dataTableOutput("outdated_result_table_mouse"),
                 textOutput("percent_mapped_mouse"),
                 downloadButton("download_outdated_mouse", "Download Outdated Ensembl Results")
               )
             )
    ),
    
    # Rat Tab
    tabPanel("Rat", value = "rat",
             sidebarLayout(
               sidebarPanel(
                 selectInput("ensembl_version_rat", "Select Ensembl Version:",
                             choices = list(
                               
                               "Ensembl v77 (Oct 2014, Rnor5.0)" = 77,
                               "Ensembl v80 (May 2015, Rnor6.0)" = 80,
                               "Ensembl v97 (Jul 2019, Rnor6.0)" = 97,
                               "Ensembl v100 (Apr 2020, Rnor6.0)" = 100,
                               "Ensembl v105 (Dec 2021, Rnor6.0)" = 105,
                               "Ensembl v110 (Jul 2023, Mrat7.2)" = 110
                             ),
                             selected = 100),
                 
                 fileInput("outdated_file_rat", "Upload Outdated Ensembl IDs (Excel)", accept = c(".xlsx")),
                 textInput("outdated_file_name_rat", "Enter file name for download:", "outdated_gene_query_results")
               ),
               mainPanel(
                 DT::dataTableOutput("outdated_result_table_rat"),
                 textOutput("percent_mapped_rat"),
                 downloadButton("download_outdated_rat", "Download Outdated Ensembl Results")
               )
             )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Calculate and display the combined percentage of successfully mapped Ensembl IDs
  calculate_percent_mapped <- function(data) {
    total_old_ids <- nrow(data)
    mapped_to_any <- sum(!is.na(data$`New_Ensembl_ID (biomaRt)`) | 
                           !is.na(data$`New_Ensembl_ID (MyGene.info)`) | 
                           !is.na(data$`New_Ensembl_ID (Ensembl REST)`))
    if (total_old_ids > 0) {
      percent_mapped <- (mapped_to_any / total_old_ids) * 100
    } else {
      percent_mapped <- 0
    }
    return(percent_mapped)
  }
  observeEvent(input$outdated_file_human, {
    req(input$outdated_file_human)
    ensembl_ids <- read_excel(input$outdated_file_human$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])
    
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    on.exit(progress$close())
    
    result_data <- process_outdated_ensembl(ensembl_ids, "human", input$ensembl_version_human, progress)  # Pass "human"
    
    # Step 3: Update the output table to include Aliases
    output$outdated_result_table_human <- DT::renderDataTable({
      DT::datatable(result_data[, c("external_gene_name", "New_Ensembl_ID_MyGene", "Aliases", "Biotype")])
    })
    
    output$percent_mapped_human <- renderText({
      paste0("Percentage of genes successfully mapped: ", 
             round(calculate_percent_mapped(result_data), 2), "%")
    })
    
    output$download_outdated_human <- downloadHandler(
      filename = function() {
        paste("outdated_ensembl_results_human", Sys.Date(), ".xlsx", sep = "")
      },
      content = function(file) {
        result_data <- process_outdated_ensembl(ensembl_ids, "human", input$ensembl_version_human, progress)
        writexl::write_xlsx(result_data, file)  # Write to file
      }
    )
  })
  
  
  observeEvent(input$outdated_file_mouse, {
    req(input$outdated_file_mouse)
    ensembl_ids <- read_excel(input$outdated_file_mouse$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])
    
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    on.exit(progress$close())
    
    result_data <- process_outdated_ensembl(ensembl_ids, "mouse", input$ensembl_version_mouse, progress)  # Pass "mouse"
    output$outdated_result_table_mouse <- DT::renderDataTable({
      DT::datatable(result_data)
    })
    
    output$percent_mapped_mouse <- renderText({
      paste0("Percentage of genes successfully mapped: ", 
             round(calculate_percent_mapped(result_data), 2), "%")
    })
    
    output$download_outdated_mouse <- downloadHandler(
      filename = function() {
        paste("outdated_ensembl_results_mouse", Sys.Date(), ".xlsx", sep = "")
      },
      content = function(file) {
        result_data <- process_outdated_ensembl(ensembl_ids, "mouse", input$ensembl_version_mouse, progress)
        writexl::write_xlsx(result_data, file)  # Write to file
      }
    )
    
  })
  
  
  # Rat Processing
  observeEvent(input$outdated_file_rat, {
    req(input$outdated_file_rat)
    ensembl_ids <- read_excel(input$outdated_file_rat$datapath)
    ensembl_ids <- as.vector(ensembl_ids[[1]])
    
    progress <- shiny::Progress$new()
    progress$set(message = "Processing outdated Ensembl IDs...", value = 0)
    on.exit(progress$close())
    
    result_data <- process_outdated_ensembl(ensembl_ids, "rat", input$ensembl_version_rat, progress)
    
    # Step 3: Update the output table to include Aliases for Rat
    output$outdated_result_table_rat <- DT::renderDataTable({
      DT::datatable(result_data[, c("external_gene_name", "New_Ensembl_ID_MyGene", "Aliases", "Biotype")])
    })
    
    output$percent_mapped_rat <- renderText({
      paste0("Percentage of genes successfully mapped: ", 
             round(calculate_percent_mapped(result_data), 2), "%")
    })
    
    output$download_outdated_rat <- downloadHandler(
      filename = function() {
        paste("outdated_ensembl_results_rat", Sys.Date(), ".xlsx", sep = "")
      },
      content = function(file) {
        result_data <- process_outdated_ensembl(ensembl_ids, "rat", input$ensembl_version_rat, progress)
        writexl::write_xlsx(result_data, file)  # Write to file
      }
    )
  })
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
