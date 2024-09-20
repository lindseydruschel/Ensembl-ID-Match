# Load necessary libraries
library(biomaRt)
library(readxl)
library(writexl)
library(httr)
library(jsonlite)
#now, it lists the ensembl from 2 diff databases and any aliases. Next step: add new databases (in GPT under ensembl gene lookup)
# and plug alaises back into the databases. 


# Step 1: Read the Excel file containing old Ensembl IDs from the "Given" column
file_path <- "C:\\Users\\druschel\\Documents\\Yulia RNASeq\\127 unmapped Yulia Ensembl.xlsx"  # Update this with the directory
old_ensembl_data <- read_excel(file_path)

# Extract the old Ensembl IDs from the column named "Given"
old_ensembl_ids <- old_ensembl_data$Given  # Adjust column name if different

# Step 2: Connect to Ensembl version 80 (Rnor6.0) to retrieve gene symbols
ensembl_v80 <- useEnsembl(biomart = "ensembl", version = 80, dataset = "rnorvegicus_gene_ensembl")

# Query for gene symbols and descriptions from version 80 using the old Ensembl IDs
results_old <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                     filters = 'ensembl_gene_id', 
                     values = old_ensembl_ids, 
                     mart = ensembl_v80)

# Ensure we keep all old Ensembl IDs in the result
results_old_full <- merge(data.frame(Old_Ensembl_ID = old_ensembl_ids), results_old, by.x = "Old_Ensembl_ID", by.y = "ensembl_gene_id", all.x = TRUE)

# Step 3: Connect to the latest Ensembl version to map to new Ensembl IDs
ensembl_latest <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

# Use the gene symbols from the old IDs to get the new Ensembl IDs
gene_symbols <- results_old$external_gene_name  # Extract gene symbols

# Query for new Ensembl IDs using gene symbols in the latest release
results_new <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                     filters = 'external_gene_name', 
                     values = gene_symbols, 
                     mart = ensembl_latest)

# Merge to keep all old Ensembl IDs even if no new Ensembl ID is found
results_new_full <- merge(results_old_full, results_new, by.x = "external_gene_name", by.y = "external_gene_name", all.x = TRUE)

# Now using mygene.info to query for aliases and Ensembl IDs for the gene symbols

# Step 4: Install and load necessary packages for mygene.info
library(httr)
library(jsonlite)

# Create a combined gene symbols list, adding some known genes to test
gene_symbolsT <- c("TP53", results_old$external_gene_name)

# Function to query mygene.info for Ensembl IDs and aliases using gene symbols
get_gene_info <- function(gene_symbol) {
  # Set species to rat (or another taxon if needed)
  species <- "rat"
  url <- paste0("https://mygene.info/v3/query?q=", gene_symbol, 
                "&fields=ensembl.gene,alias,symbol&species=", species)
  
  # Make the GET request
  response <- GET(url)
  
  # Parse the response
  content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
  
  # Extract the Ensembl ID and alias if available
  ensembl_id <- if (!is.null(content$hits$ensembl$gene)) content$hits$ensembl$gene[1] else NA
  aliases <- if (!is.null(content$hits$alias)) paste(content$hits$alias, collapse = ", ") else NA
  
  # Return the results as a data frame
  return(data.frame(Gene_Symbol = gene_symbol, Ensembl_ID = ensembl_id, Aliases = aliases))
}

# Step 5: Query for all gene symbols to get their Ensembl IDs and aliases from mygene.info
gene_info_results <- do.call(rbind, lapply(gene_symbolsT, get_gene_info))

# Step 6: Combine results from biomaRt and mygene.info, ensuring all old IDs are included
final_results <- merge(results_new_full, gene_info_results, by.x = "external_gene_name", by.y = "Gene_Symbol", all.x = TRUE)

# View the final combined results (with all input genes, even if some have missing data)
print(final_results)

# Step 7: Save the final combined results to an Excel file
output_file <- "Combined_Ensembl_Aliases_Full.xlsx"
write_xlsx(final_results, output_file)

# Success message
cat("Combined results saved to", output_file, "\n")
