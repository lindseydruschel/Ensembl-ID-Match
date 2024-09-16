#README/Notes: This R script reads old Ensembl gene IDs from the "Given" column in an Excel file (currently 
#"named 2000 unmapped ensemble yulia""), maps them to the
#new Ensembl gene IDs using the biomaRt package, and saves the results to a new Excel file. Ensure the input
#file is named "Unmatched Ensembl.xlsx" and contains a column titled "Given" with the old IDs. The output 
#will be saved as "Updated Ensembl.xlsx" in the same directory. Typically not every old ensemble ID will have
#a new ID (not sure why as of rn) but they should all have a gene symbol if you're using the correct reference 
#database. Even within a release (Rnor 6.0 in this case), there are different versions. Not every version is in
#biomaRt, for example I wanted version 84 but that's not available so I used 80 (the closest version under Rnor 6.0).
#To see all available versions use listEnsemblArchives().WARNING: if you rerun the code, Updated Ensembl.xlsx will
#be overwrite the current version, so if you have important into in that sheet save it under a new name before running. 


library(biomaRt)
library(readxl)
library(writexl)


# Read the Excel file containing old Ensembl IDs from the "Given" column
file_path <- "C:\\Users\\druschel\\Documents\\Yulia RNASeq\\2000 unmapped ensemble yulia.xlsx"  # Update this with the directory
old_ensembl_data <- read_excel(file_path)

# Extract the old Ensembl IDs from the column named "Given"
old_ensembl_ids <- old_ensembl_data$Given  # Using the column name "Given"

# Connect to Ensembl version 80 (Rnor6.0)
ensembl_v80 <- useEnsembl(biomart = "ensembl", version = 80, dataset = "rnorvegicus_gene_ensembl")

# Query for gene symbols and descriptions from version 80 using the old Ensembl IDs
results_old <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                     filters = 'ensembl_gene_id', 
                     values = old_ensembl_ids, 
                     mart = ensembl_v80)

# View the result to make sure it works
print(results_old)

# Connect to the latest Ensembl version
ensembl_latest <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

# Use the gene symbols from the old IDs to get the new Ensembl IDs
gene_symbols <- results_old$external_gene_name  # Extract gene symbols

# Query for new Ensembl IDs using gene symbols in the latest release
results_new <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                     filters = 'external_gene_name', 
                     values = gene_symbols, 
                     mart = ensembl_latest)

# View the result
print(results_new)

# Step 4: Merge the old IDs with new Ensembl IDs (ensure no rows are dropped)
merged_results <- merge(results_old, results_new, by = "external_gene_name", all.x = TRUE)

output_file <- "Updated Ensembl.xlsx"
write_xlsx(merged_results, output_file)

# Success message
cat("Merged results saved to", output_file, "\n")

# View the final results
print(merged_results)

# Calculate the percentage mapped to gene symbols and new Ensembl IDs
total_old_ids <- length(old_ensembl_ids)
mapped_to_symbol <- sum(!is.na(results_old$external_gene_name))
mapped_to_new_id <- sum(!is.na(merged_results$ensembl_gene_id.y))

percent_mapped_to_symbol <- (mapped_to_symbol / total_old_ids) * 100
percent_mapped_to_new_id <- (mapped_to_new_id / total_old_ids) * 100

# Display the percentages
cat("Percentage mapped to a gene symbol:", percent_mapped_to_symbol, "%\n")
cat("Percentage mapped to a new Ensembl ID:", percent_mapped_to_new_id, "%\n")

# Step 6: Write the result to a new Excel file
output_file <- "Updated Ensembl.xlsx"
write_xlsx(merged_results, output_file)

# Success message
cat("Results saved to", output_file, "\n")






##STOP HERE! THIS IS EXTRA

## Checks for duplicated symbols from ol Ensembl, don't need these for results
#Test mult mappings
results_new_dup <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                     filters = 'external_gene_name', 
                     values = gene_symbols, 
                     mart = ensembl_latest)

# Check if multiple mappings exist for the same symbol
duplicated_symbols <- results_new_dup[duplicated(results_new$external_gene_name), ]
print(duplicated_symbols)

##Clear everything
rm(list = ls())
# Clear the console
cat("\014")
