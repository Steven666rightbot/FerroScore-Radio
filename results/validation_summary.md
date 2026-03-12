
## Available Immunotherapy Datasets

### GSE78220 (Melanoma, Anti-PD-1)
- Samples: 28
- Responders: 15 (CR/PR)
- Non-responders: 13 (PD)
- Status: Response labels extracted
- Missing: Expression matrix (need SRA download)

### GSE91061 (Melanoma, Anti-PD-1)
- Samples: 109
- Pre/On treatment pairs
- Response labels available
- Status: Phenotype downloaded
- Missing: Expression matrix (need SRA download)

### IMvigor210 (Urothelial, Anti-PD-L1)
- Samples: 348
- Status: Need to download via R package

## Next Steps

1. Download expression data from SRA for GSE78220/GSE91061
   - Use SRA Toolkit: prefetch + fasterq-dump
   - Align with STAR/HISAT2
   - Quantify with featureCounts/salmon

2. Or use pre-processed data from:
   - UCSC Xena (if available)
   - GEO2R (for microarray data)
   - Published supplementary files

3. Calculate FerroScore-Immuno for validation samples

4. Evaluate prediction performance (AUC, accuracy)

## Alternative: Use TCGA with Proxy Labels

Since immunotherapy data is limited, use TCGA with:
- High FerroImmuno Score = Predicted responder
- TMB (Tumor Mutation Burden) as additional filter
- Immune infiltration scores
