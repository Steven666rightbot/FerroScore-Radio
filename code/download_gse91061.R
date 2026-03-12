#!/usr/bin/env Rscript
# Download and analyze GSE91061 using GEOquery

# Install GEOquery if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery")

library(GEOquery)

# Set working directory
setwd("D:/Research/FerroScore-Radio/data/external")

cat("========================================\n")
cat("Downloading GSE91061\n")
cat("========================================\n\n")

# Download GSE91061
tryCatch({
    gse <- getGEO("GSE91061", GSEMatrix = TRUE, getGPL = FALSE)
    
    cat("Download successful!\n")
    cat("Number of samples:", length(gse), "\n")
    
    # Extract expression data
    exprs_data <- exprs(gse[[1]])
    cat("Expression matrix dimensions:", dim(exprs_data), "\n")
    
    # Extract phenotype data
    pheno_data <- pData(gse[[1]])
    cat("Phenotype data dimensions:", dim(pheno_data), "\n")
    cat("Phenotype columns:", names(pheno_data)[1:10], "...\n")
    
    # Check for response information
    cat("\nSearching for response information...\n")
    response_cols <- grep("response|responder|outcome|therapy", names(pheno_data), 
                          ignore.case = TRUE, value = TRUE)
    cat("Potential response columns:", response_cols, "\n")
    
    # Save expression data
    write.csv(exprs_data, "GSE91061_expression.csv")
    cat("\nSaved: GSE91061_expression.csv\n")
    
    # Save phenotype data
    write.csv(pheno_data, "GSE91061_phenotype.csv")
    cat("Saved: GSE91061_phenotype.csv\n")
    
    cat("\n========================================\n")
    cat("Download and extraction complete!\n")
    cat("========================================\n")
    
}, error = function(e) {
    cat("Error:", e$message, "\n")
    cat("\nPossible solutions:\n")
    cat("1. Check internet connection\n")
    cat("2. Try downloading manually from GEO\n")
    cat("3. Use alternative dataset\n")
})
