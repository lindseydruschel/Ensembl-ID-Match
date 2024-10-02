# Load necessary packages function
load_packages <- function(pkgs) {
  sapply(pkgs, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  })
}

packages <- c("shiny", "biomaRt", "httr", "jsonlite", "writexl", "DT", "KEGGREST", "readxl")
load_packages(packages)

# Ensure required packages are installed
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("DT", quietly = TRUE)) {
  install.packages("DT", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("shinyjs", quietly = TRUE)) {
  install.packages("shinyjs", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl", repos = "https://cran.rstudio.com/")
}

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
# Function to retrieve KEGG ID based on pathway name
get_kegg_id_from_name <- function(pathway_name) {
  kegg_data <- KEGGREST::keggList("pathway")  # Get list of all pathways
  kegg_pathways <- names(kegg_data)  # Extract KEGG IDs
  
  # Find KEGG ID for the given pathway name
  for (id in kegg_pathways) {
    if (grepl(pathway_name, kegg_data[[id]], ignore.case = TRUE)) {
      return(sub("pathway:", "", id))  # Return the KEGG ID without "pathway:"
    }
  }
  return(NA)  # Return NA if no match is found
}
# Function to get Ensembl ID based on gene symbol
get_ensembl_id <- function(gene_symbol) {
  result <- ensembl_mapping$Ensembl_ID[ensembl_mapping$Gene_Symbol == gene_symbol]
  if (length(result) > 0) {
    return(result[1])  # Return the first matching Ensembl ID
  } else {
    return(NA)  # Return NA if not found
  }
}

# Function to get GO ID based on GO term
get_go_id <- function(go_term) {
  result <- go_mapping$GO_ID[go_mapping$GO_Term == go_term]
  if (length(result) > 0) {
    return(result[1])  # Return the first matching GO ID
  } else {
    return(NA)  # Return NA if not found
  }
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
  titlePanel("Gene Search Dashboard"),
  
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
                 conditionalPanel(
                   condition = "input.species == 'human'",
                   radioButtons("output_type_human", "Table Output:",
                                choices = list("Plain text (quick)" = "plain", 
                                               "NCBI/KEGG Links (slow)" = "links"),
                                selected = "plain")
                 ),
                 
                 textInput("file_name_human", "Enter file name for download:", "gene_query_results"),
                 textInput("keyword_human", "Keyword Search (GO Only):", ""),
                 
                 helpText("Enter gene symbols or Ensembl IDs separated by commas, or upload an Excel file.")
               )
               ,
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
                 conditionalPanel(
                   condition = "input.species == 'mouse'",
                   radioButtons("output_type_mouse", "Table Output:",
                                choices = list("Plain text (quick)" = "plain", 
                                               "NCBI/KEGG Links (slow)" = "links"),
                                selected = "plain")
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
                 conditionalPanel(
                   condition = "input.species == 'rat'",  # Ensure this is within the Rat tab
                   radioButtons("output_type_rat", "Table Output:",
                                choices = list("Plain text (quick)" = "plain", 
                                               "NCBI/KEGG Links (slow)" = "links"),
                                selected = "plain")
                 ),
                 textInput("file_name_rat", "Enter file name for download:", "gene_query_results"),
                 textInput("keyword_rat", "Keyword Search (GO Only):", ""),
                 helpText("Enter gene symbols or Ensembl IDs separated by commas, or upload an Excel file.")
               ),
               mainPanel(
                 uiOutput("progress_bar_rat"),
                 DT::dataTableOutput("result_table_rat"),
                 downloadButton("download_data_rat", "Download Full Results"),
                 downloadButton("download_filtered_rat", "Download Filtered Results"),
                 h3("Gene Information"),
                 verbatimTextOutput("gene_info_rat"),
                 tags$div(
                   style = "margin-top: 20px;",
                   tags$strong("Note:"),
                   tags$p("Gene symbol links in the rat tab direct to human genes. This is due to the lack of comprehensive gene information available for rat genes in the selected database. We recommend using human gene links for further investigation.")
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
  
  # Observe block for human species
  observe({
    human_results <- process_input("human", input$input_type_human, input$input_source_human, 
                                   input$user_input_human, input$ensembl_input_human, 
                                   input$file_input_human, input$keyword_human)
    results_human(human_results)  # Update the reactive value for human
  })
  
  # Observe block for mouse species
  observe({
    mouse_results <- process_input("mouse", input$input_type_mouse, input$input_source_mouse, 
                                   input$user_input_mouse, input$ensembl_input_mouse, 
                                   input$file_input_mouse, input$keyword_mouse)
    results_mouse(mouse_results)  # Update the reactive value for mouse
  })
  
  # Observe block for rat species
  observe({
    rat_results <- process_input("rat", input$input_type_rat, input$input_source_rat, 
                                 input$user_input_rat, input$ensembl_input_rat, 
                                 input$file_input_rat, input$keyword_rat)
    results_rat(rat_results)  # Update the reactive value for rat
  })
  
  # Reactive values to handle different species results
  
  output$result_table_human <- DT::renderDataTable({
    result_data <- process_input("human", input$input_type_human, input$input_source_human, 
                                 input$user_input_human, input$ensembl_input_human, 
                                 input$file_input_human, input$keyword_human)
    
    # Ensure result_data is a data frame
    if (!is.data.frame(result_data)) {
      result_data <- as.data.frame(result_data)
    }
    
    # Print the structure of result_data for debugging
    print(str(result_data))  # Check the structure before modification
    
    if (input$output_type_human == "plain") {
      # Return plain text version
      result_data$Gene_Symbol <- as.character(result_data$Gene_Symbol)
      result_data$Aliases <- as.character(result_data$Aliases)  # Keep as plain text
      result_data$Ensembl_ID <- as.character(result_data$Ensembl_ID)  # Keep as plain text
      result_data$KEGG_Pathways <- as.character(result_data$KEGG_Pathways)  # Keep as plain text
    } else {  # If user selected links
      # Create clickable links for the Gene_Symbol column linking to NCBI
      if ("Gene_Symbol" %in% colnames(result_data)) {
        result_data$Gene_Symbol <- sapply(as.character(result_data$Gene_Symbol), function(gene) {
          if (!is.na(gene) && nchar(gene) > 0) {
            link <- paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'>", gene, "</a>")
            return(link)  # Return the link as a string
          } else {
            return(gene)  # Return original if NA or empty
          }
        })
      }
      
      # Create clickable links for the Ensembl_ID column linking to Ensembl
      if ("Ensembl_ID" %in% colnames(result_data)) {
        result_data$Ensembl_ID <- sapply(as.character(result_data$Ensembl_ID), function(ensembl_id) {
          if (!is.na(ensembl_id) && nchar(ensembl_id) > 0) {
            link <- paste0("<a href='https://www.ensembl.org/id/", ensembl_id, "' target='_blank'>", ensembl_id, "</a>")
            return(link)  # Return the link as a string
          } else {
            return(ensembl_id)  # Return original if NA or empty
          }
        })
      }
      
      # Keep Aliases as plain text
      if ("Aliases" %in% colnames(result_data)) {
        result_data$Aliases <- as.character(result_data$Aliases)  # Keep as plain text
      }
      
      # Create clickable links for the KEGG_Pathways column linking to KEGG
      if ("KEGG_Pathways" %in% colnames(result_data)) {
        result_data$KEGG_Pathways <- sapply(as.character(result_data$KEGG_Pathways), function(pathways) {
          if (is.na(pathways) || nchar(pathways) == 0) {
            return("")  # Return empty if NA or missing
          }
          pathway_list <- strsplit(pathways, ",")[[1]]
          pathway_list <- trimws(pathway_list)
          
          links <- sapply(pathway_list, function(pathway) {
            if (nchar(pathway) > 0) {
              kegg_id <- get_kegg_id_from_name(pathway)
              if (!is.na(kegg_id)) {
                link <- paste0("<a href='https://www.kegg.jp/pathway/", kegg_id, "' target='_blank'>", pathway, "</a>")
                return(link)
              } else {
                return(pathway)
              }
            } else {
              return(pathway)
            }
          })
          return(paste(links, collapse = ", "))
        })
      }
    }
    
    # Render the DataTable with escape = FALSE to allow clickable HTML links
    DT::datatable(result_data, escape = FALSE, options = list(pageLength = 5))
  })
  
  
  output$result_table_mouse <- DT::renderDataTable({
    result_data <- process_input("mouse", input$input_type_mouse, input$input_source_mouse, 
                                 input$user_input_mouse, input$ensembl_input_mouse, 
                                 input$file_input_mouse, input$keyword_mouse)
    
    # Ensure result_data is a data frame
    if (!is.data.frame(result_data)) {
      result_data <- as.data.frame(result_data)
    }
    
    # Print the structure of result_data for debugging
    print(str(result_data))  # Check the structure before modification
    
    if (input$output_type_mouse == "plain") {
      # Return plain text version
      result_data$Gene_Symbol <- as.character(result_data$Gene_Symbol)
      result_data$Aliases <- as.character(result_data$Aliases)  # Keep as plain text
      result_data$Ensembl_ID <- as.character(result_data$Ensembl_ID)  # Keep as plain text
      result_data$KEGG_Pathways <- as.character(result_data$KEGG_Pathways)  # Keep as plain text
    } else {  # If user selected links
      # Create clickable links for the Gene_Symbol column linking to NCBI
      if ("Gene_Symbol" %in% colnames(result_data)) {
        result_data$Gene_Symbol <- sapply(as.character(result_data$Gene_Symbol), function(gene) {
          if (!is.na(gene) && nchar(gene) > 0) {
            link <- paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'>", gene, "</a>")
            return(link)  # Return the link as a string
          } else {
            return(gene)  # Return original if NA or empty
          }
        })
      }
      
      # Create clickable links for the Ensembl_ID column linking to Ensembl
      if ("Ensembl_ID" %in% colnames(result_data)) {
        result_data$Ensembl_ID <- sapply(as.character(result_data$Ensembl_ID), function(ensembl_id) {
          if (!is.na(ensembl_id) && nchar(ensembl_id) > 0) {
            link <- paste0("<a href='https://www.ensembl.org/id/", ensembl_id, "' target='_blank'>", ensembl_id, "</a>")
            return(link)  # Return the link as a string
          } else {
            return(ensembl_id)  # Return original if NA or empty
          }
        })
      }
      
      # Keep Aliases as plain text
      if ("Aliases" %in% colnames(result_data)) {
        result_data$Aliases <- as.character(result_data$Aliases)  # Keep as plain text
      }
      
      # Create clickable links for the KEGG_Pathways column linking to KEGG
      if ("KEGG_Pathways" %in% colnames(result_data)) {
        result_data$KEGG_Pathways <- sapply(as.character(result_data$KEGG_Pathways), function(pathways) {
          if (is.na(pathways) || nchar(pathways) == 0) {
            return("")  # Return empty if NA or missing
          }
          pathway_list <- strsplit(pathways, ",")[[1]]
          pathway_list <- trimws(pathway_list)
          
          links <- sapply(pathway_list, function(pathway) {
            if (nchar(pathway) > 0) {
              kegg_id <- get_kegg_id_from_name(pathway)
              if (!is.na(kegg_id)) {
                link <- paste0("<a href='https://www.kegg.jp/pathway/", kegg_id, "' target='_blank'>", pathway, "</a>")
                return(link)
              } else {
                return(pathway)
              }
            } else {
              return(pathway)
            }
          })
          return(paste(links, collapse = ", "))
        })
      }
    }
    
    # Render the DataTable with escape = FALSE to allow clickable HTML links
    DT::datatable(result_data, escape = FALSE, options = list(pageLength = 5))
  })
  
  output$result_table_rat <- DT::renderDataTable({
    result_data <- results_rat()  # Use the reactive results for Rat
    
    # Ensure result_data is a data frame
    if (!is.data.frame(result_data)) {
      result_data <- as.data.frame(result_data)
    }
    
    # Print the structure of result_data for debugging
    print(str(result_data))  # Check the structure before modification
    
    if (input$output_type_rat == "plain") {
      # Return plain text version
      result_data$Gene_Symbol <- as.character(result_data$Gene_Symbol)
      result_data$Aliases <- as.character(result_data$Aliases)  # Keep as plain text
      result_data$Ensembl_ID <- as.character(result_data$Ensembl_ID)  # Keep as plain text
      result_data$KEGG_Pathways <- as.character(result_data$KEGG_Pathways)  # Keep as plain text
      
    } else {  # If user selected links
      # Create clickable links for the Gene_Symbol column linking to NCBI
      if ("Gene_Symbol" %in% colnames(result_data)) {
        result_data$Gene_Symbol <- sapply(as.character(result_data$Gene_Symbol), function(gene) {
          if (!is.na(gene) && nchar(gene) > 0) {
            link <- paste0("<a href='https://www.ncbi.nlm.nih.gov/gene/?term=", gene, "' target='_blank'>", gene, "</a>")
            return(link)  # Return the link as a string
          } else {
            return(gene)  # Return original if NA or empty
          }
        })
      }
      
      # Create clickable links for the Ensembl_ID column linking to Ensembl
      if ("Ensembl_ID" %in% colnames(result_data)) {
        result_data$Ensembl_ID <- sapply(as.character(result_data$Ensembl_ID), function(ensembl_id) {
          if (!is.na(ensembl_id) && nchar(ensembl_id) > 0) {
            link <- paste0("<a href='https://www.ensembl.org/id/", ensembl_id, "' target='_blank'>", ensembl_id, "</a>")
            return(link)  # Return the link as a string
          } else {
            return(ensembl_id)  # Return original if NA or empty
          }
        })
      }
      
      # Keep Aliases as plain text
      if ("Aliases" %in% colnames(result_data)) {
        result_data$Aliases <- as.character(result_data$Aliases)  # Keep as plain text
      }
      
      # Create clickable links for the KEGG_Pathways column linking to KEGG
      if ("KEGG_Pathways" %in% colnames(result_data)) {
        result_data$KEGG_Pathways <- sapply(as.character(result_data$KEGG_Pathways), function(pathways) {
          if (is.na(pathways) || nchar(pathways) == 0) {
            return("")  # Return empty if NA or missing
          }
          pathway_list <- strsplit(pathways, ",")[[1]]
          pathway_list <- trimws(pathway_list)
          
          links <- sapply(pathway_list, function(pathway) {
            if (nchar(pathway) > 0) {
              kegg_id <- get_kegg_id_from_name(pathway)
              if (!is.na(kegg_id)) {
                link <- paste0("<a href='https://www.kegg.jp/pathway/", kegg_id, "' target='_blank'>", pathway, "</a>")
                return(link)
              } else {
                return(pathway)
              }
            } else {
              return(pathway)
            }
          })
          return(paste(links, collapse = ", "))
        })
      }
    }
    
    # Render the DataTable with escape = FALSE to allow clickable HTML links
    DT::datatable(result_data, escape = FALSE, options = list(pageLength = 5))
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
      data <- results_rat()  # Retrieve the updated reactive value
      if (nrow(data) == 0) {  # Check if data is empty
        showNotification("No data available for download", type = "error")
        return()
      }
      write_xlsx(data, file)
    }
  )
  
  # Download handlers for filtered results
  output$download_filtered_human <- downloadHandler(
    filename = function() {
      paste(input$file_name_human, "_filtered.xlsx", sep = "")
    },
    content = function(file) {
      filtered_data <- results_human()[input$result_table_human_rows_all, ]
      if (nrow(filtered_data) == 0) {  # Check if filtered data is empty
        showNotification("No filtered data available for download", type = "error")
        return()
      }
      write_xlsx(filtered_data, file)
    }
  )
  
  output$download_filtered_mouse <- downloadHandler(
    filename = function() {
      paste(input$file_name_mouse, "_filtered.xlsx", sep = "")
    },
    content = function(file) {
      filtered_data <- results_mouse()[input$result_table_mouse_rows_all, ]
      if (nrow(filtered_data) == 0) {  # Check if filtered data is empty
        showNotification("No filtered data available for download", type = "error")
        return()
      }
      write_xlsx(filtered_data, file)
    }
  )
  
  
  output$download_filtered_rat <- downloadHandler(
    filename = function() {
      paste(input$file_name_rat, "_filtered.xlsx", sep = "")
    },
    content = function(file) {
      filtered_data <- results_rat()[input$result_table_rat_rows_all, ]
      if (nrow(filtered_data) == 0) {  # Check if filtered data is empty
        showNotification("No filtered data available for download", type = "error")
        return()
      }
      write_xlsx(filtered_data, file)
    }
  )
}
# Run the application
shinyApp(ui = ui, server = server)
