#!/usr/bin/env python3
"""
Download and analyze GSE91061 using GEOparse
"""

import GEOparse
import pandas as pd
import os

EXTERNAL_DIR = '../data/external'
os.chdir(EXTERNAL_DIR)

print("=" * 60)
print("Downloading GSE91061 using GEOparse")
print("=" * 60)

try:
    # Download GSE91061
    print("\nDownloading... This may take a few minutes.")
    gse = GEOparse.get_GEO(geo="GSE91061", destdir=".", silent=True)
    
    print(f"\nDownload successful!")
    print(f"GSE title: {gse.metadata.get('title', ['N/A'])[0]}")
    
    # Get GSMs (samples)
    gsms = gse.gsms
    print(f"Number of samples: {len(gsms)}")
    
    # Check first sample
    first_gsm = list(gsms.values())[0]
    print(f"\nFirst sample: {first_gsm.name}")
    print(f"Metadata keys: {list(first_gsm.metadata.keys())[:10]}")
    
    # Try to get expression data
    print("\nAttempting to extract expression data...")
    
    # Check if expression data is available
    expr_data = []
    sample_names = []
    
    for gsm_name, gsm in list(gsms.items())[:5]:  # Check first 5 samples
        if hasattr(gsm, 'table') and gsm.table is not None:
            print(f"  {gsm_name}: Has table data")
            if not expr_data:
                # First sample with data
                expr_data = gsm.table
                sample_names.append(gsm_name)
        else:
            print(f"  {gsm_name}: No table data")
    
    if len(sample_names) > 0:
        print(f"\nExpression data found for {len(sample_names)} samples")
        print(f"Expression matrix shape: {expr_data.shape}")
        
        # Save expression data
        expr_data.to_csv('GSE91061_expression.csv', index=False)
        print("Saved: GSE91061_expression.csv")
    else:
        print("\nNo expression table data found in samples.")
        print("This is RNA-seq data. Expression data may need to be downloaded from SRA.")
    
    # Extract and save phenotype data
    print("\nExtracting phenotype data...")
    pheno_data = []
    
    for gsm_name, gsm in gsms.items():
        pheno_row = {'sample_id': gsm_name}
        
        # Extract characteristics
        for key, value in gsm.metadata.items():
            if isinstance(value, list):
                pheno_row[key] = '; '.join(value)
            else:
                pheno_row[key] = value
        
        pheno_data.append(pheno_row)
    
    pheno_df = pd.DataFrame(pheno_data)
    print(f"Phenotype data shape: {pheno_df.shape}")
    print(f"Phenotype columns: {list(pheno_df.columns)[:15]}")
    
    # Save phenotype data
    pheno_df.to_csv('GSE91061_phenotype.csv', index=False)
    print("Saved: GSE91061_phenotype.csv")
    
    # Search for response information
    print("\nSearching for response information...")
    response_keywords = ['response', 'responder', 'outcome', 'therapy', 'treatment']
    
    for col in pheno_df.columns:
        if any(keyword in col.lower() for keyword in response_keywords):
            print(f"  Potential response column: {col}")
            print(f"    Unique values: {pheno_df[col].unique()[:10]}")
    
    print("\n" + "=" * 60)
    print("Download and extraction complete!")
    print("=" * 60)

except Exception as e:
    print(f"\nError: {e}")
    print("\nPossible solutions:")
    print("1. Check internet connection")
    print("2. Dataset may require manual download from GEO")
    print("3. Try alternative dataset")
