#!/usr/bin/env python3
"""
Analyze GSE35452 - Rectal Cancer Radiotherapy Response
"""

import pandas as pd
import numpy as np
import gzip
import os

# Paths
EXTERNAL_DIR = '../data/external'
RESULTS_DIR = '../results'
os.makedirs(f'{RESULTS_DIR}/external', exist_ok=True)

print("=" * 60)
print("GSE35452 Analysis - Rectal Cancer Radiotherapy Response")
print("=" * 60)

# Read GSE35452 data
file_path = f'{EXTERNAL_DIR}/GSE35452_series_matrix.txt.gz'

print(f"\nReading: {file_path}")

# Parse series matrix file
sample_metadata = {}
expression_rows = []
in_data_section = False

with gzip.open(file_path, 'rt', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        
        if not line:
            continue
            
        # Parse metadata
        if line.startswith('!Sample_'):
            parts = line.split('\t')
            key = parts[0].replace('!Sample_', '')
            values = [p.strip('"') for p in parts[1:]]
            sample_metadata[key] = values
            
        # Start of expression data
        elif line.startswith('"ID_REF"'):
            in_data_section = True
            header = line.split('\t')
            sample_ids = [h.strip('"') for h in header[1:]]
            
        elif in_data_section and not line.startswith('!'):
            parts = line.split('\t')
            gene_id = parts[0].strip('"')
            values = [float(p) if p != 'null' else np.nan for p in parts[1:]]
            expression_rows.append([gene_id] + values)

print(f"Samples: {len(sample_ids)}")
print(f"Metadata fields: {list(sample_metadata.keys())}")

# Create expression DataFrame
expr_df = pd.DataFrame(expression_rows, columns=['Gene'] + sample_ids)
expr_df = expr_df.set_index('Gene')
print(f"Expression matrix: {expr_df.shape}")

# Extract sample characteristics
if 'characteristics_ch1' in sample_metadata:
    characteristics = sample_metadata['characteristics_ch1']
    print(f"\nSample characteristics:")
    for i, char in enumerate(characteristics[:5]):
        print(f"  {sample_ids[i]}: {char}")

# Try to identify response groups
response_info = []
for i, char in enumerate(characteristics):
    if 'response' in char.lower() or 'responder' in char.lower():
        response_info.append((sample_ids[i], char))

if response_info:
    print(f"\nFound response information for {len(response_info)} samples")
    for sid, resp in response_info[:5]:
        print(f"  {sid}: {resp}")
else:
    print("\nNo explicit response info found in characteristics")
    print("May need to extract from sample titles")

# Check sample titles for response info
if 'title' in sample_metadata:
    titles = sample_metadata['title']
    print(f"\nSample titles (first 5):")
    for i, title in enumerate(titles[:5]):
        print(f"  {sample_ids[i]}: {title}")

# Extract response labels
print("\n" + "=" * 60)
print("Extracting Response Labels")
print("=" * 60)

response_labels = {}
for i, char in enumerate(characteristics):
    sample_id = sample_ids[i]
    if 'responder' in char.lower() or 'response' in char.lower():
        if 'non-responder' in char.lower() or 'non responder' in char.lower():
            response_labels[sample_id] = 0  # Non-responder
        elif 'responder' in char.lower() and 'non' not in char.lower():
            response_labels[sample_id] = 1  # Responder
        else:
            response_labels[sample_id] = -1  # Unknown

responder_count = sum(1 for v in response_labels.values() if v == 1)
non_responder_count = sum(1 for v in response_labels.values() if v == 0)
unknown_count = sum(1 for v in response_labels.values() if v == -1)

print(f"Responder: {responder_count}")
print(f"Non-responder: {non_responder_count}")
print(f"Unknown: {unknown_count}")

# Load gene sets
print("\n" + "=" * 60)
print("Loading Ferro-Radio Gene Sets")
print("=" * 60)

gene_sets = {
    'Ferroptosis Driver': ['ACSL4', 'LPCAT3', 'ALOX15', 'ALOX5', 'NOX1', 'NOX4', 'P53', 'SAT1', 'CARS1'],
    'Ferroptosis Suppressor': ['GPX4', 'SLC7A11', 'SLC3A2', 'NFE2L2', 'HMOX1', 'FTH1', 'FTL', 'SOD1', 'SOD2', 'GCLC', 'GCLM'],
    'DNA Repair': ['BRCA1', 'BRCA2', 'ATM', 'ATR', 'CHEK1', 'CHEK2', 'TP53', 'RAD51', 'PARP1'],
    'ROS-related': ['NOX1', 'NOX2', 'NOX4', 'SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1']
}

all_genes = []
for genes in gene_sets.values():
    all_genes.extend(genes)
all_genes = list(set(all_genes))

print(f"Total Ferro-Radio genes: {len(all_genes)}")

# Check which genes are available in the dataset
available_genes = [g for g in all_genes if g in expr_df.index]
missing_genes = [g for g in all_genes if g not in expr_df.index]

print(f"Available genes: {len(available_genes)}/{len(all_genes)}")
print(f"Missing genes: {missing_genes}")

if len(available_genes) == 0:
    print("\nWARNING: No Ferro-Radio genes found!")
    print("This may be because the dataset uses different gene ID format.")
    print("Checking first few gene IDs in dataset...")
    print(list(expr_df.index[:10]))
    
    # Try to convert using GPL570 annotation
    print("\n" + "=" * 60)
    print("Attempting to convert Probe ID to Gene Symbol")
    print("=" * 60)
    
    gpl_file = f'{EXTERNAL_DIR}/GPL570-55999.txt'
    if os.path.exists(gpl_file):
        print(f"Loading annotation from {gpl_file}...")
        annotation = pd.read_csv(gpl_file, sep='\t', skiprows=16, low_memory=False)
        probe_to_gene = dict(zip(annotation['ID'], annotation['Gene Symbol']))
        
        # Convert expression matrix index
        expr_df['Gene_Symbol'] = expr_df.index.map(probe_to_gene)
        
        # Remove unmapped probes
        expr_df = expr_df[expr_df['Gene_Symbol'].notna()]
        expr_df = expr_df[expr_df['Gene_Symbol'] != '']
        
        # Handle multiple probes per gene (take mean)
        expr_df = expr_df.groupby('Gene_Symbol').mean()
        
        print(f"After conversion: {expr_df.shape}")
        
        # Check again for Ferro-Radio genes
        available_genes = [g for g in all_genes if g in expr_df.index]
        print(f"Available Ferro-Radio genes after conversion: {len(available_genes)}/{len(all_genes)}")
        
        if len(available_genes) > 0:
            print(f"Available: {available_genes}")
else:
    print(f"\nAvailable genes: {available_genes[:10]}...")

# Calculate FerroScore for all samples
print("\n" + "=" * 60)
print("Calculating FerroScore")
print("=" * 60)

# Extract available genes by category
driver_genes = [g for g in gene_sets['Ferroptosis Driver'] if g in expr_df.index]
suppressor_genes = [g for g in gene_sets['Ferroptosis Suppressor'] if g in expr_df.index]
ddr_genes = [g for g in gene_sets['DNA Repair'] if g in expr_df.index]
ros_genes = [g for g in gene_sets['ROS-related'] if g in expr_df.index]

print(f"Driver genes available: {len(driver_genes)}")
print(f"Suppressor genes available: {len(suppressor_genes)}")
print(f"DDR genes available: {len(ddr_genes)}")

# Calculate scores
if len(driver_genes) > 0:
    driver_score = expr_df.loc[driver_genes].mean()
else:
    driver_score = pd.Series(0.5, index=expr_df.columns)

if len(suppressor_genes) > 0:
    suppressor_score = 1 - expr_df.loc[suppressor_genes].mean()  # Inverse
else:
    suppressor_score = pd.Series(0.5, index=expr_df.columns)

if len(ddr_genes) > 0:
    ddr_score = expr_df.loc[ddr_genes].mean()
else:
    ddr_score = pd.Series(0.5, index=expr_df.columns)

# Combined FerroScore
ferroscore = (driver_score + suppressor_score) / 2

# Normalize to 0-1
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
ferroscore_norm = pd.Series(
    scaler.fit_transform(ferroscore.values.reshape(-1, 1)).flatten(),
    index=ferroscore.index
)
ddr_norm = pd.Series(
    scaler.fit_transform(ddr_score.values.reshape(-1, 1)).flatten(),
    index=ddr_score.index
)

# FerroRadio Score
ferro_radio = 0.7 * ferroscore_norm + 0.3 * (1 - ddr_norm)

# Create results DataFrame
results_df = pd.DataFrame({
    'FerroScore': ferroscore_norm,
    'DDR_Score': ddr_norm,
    'FerroRadio_Score': ferro_radio,
    'Response': [response_labels.get(sid, -1) for sid in expr_df.columns]
})

# Remove samples with unknown response
results_df = results_df[results_df['Response'] != -1]

print(f"\nCalculated scores for {len(results_df)} samples")
print(f"Responder: {(results_df['Response'] == 1).sum()}")
print(f"Non-responder: {(results_df['Response'] == 0).sum()}")

print("\nScore statistics:")
print(results_df.describe())

# Compare Responder vs Non-responder
print("\n" + "=" * 60)
print("Comparing Responder vs Non-responder")
print("=" * 60)

responder_scores = results_df[results_df['Response'] == 1]['FerroRadio_Score']
non_responder_scores = results_df[results_df['Response'] == 0]['FerroRadio_Score']

print(f"\nResponder (n={len(responder_scores)}):")
print(f"  Mean FerroRadio Score: {responder_scores.mean():.3f}")
print(f"  Std: {responder_scores.std():.3f}")

print(f"\nNon-responder (n={len(non_responder_scores)}):")
print(f"  Mean FerroRadio Score: {non_responder_scores.mean():.3f}")
print(f"  Std: {non_responder_scores.std():.3f}")

# T-test
from scipy import stats
t_stat, p_value = stats.ttest_ind(responder_scores, non_responder_scores)
print(f"\nT-test: t={t_stat:.3f}, p={p_value:.4f}")

if p_value < 0.05:
    print("*** SIGNIFICANT DIFFERENCE ***")
else:
    print("No significant difference")

# Train machine learning model
print("\n" + "=" * 60)
print("Training Machine Learning Model")
print("=" * 60)

from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, roc_curve, classification_report

# Features and labels
X = results_df[['FerroScore', 'DDR_Score', 'FerroRadio_Score']]
y = results_df['Response']

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)

print(f"Training set: {X_train.shape}")
print(f"Test set: {X_test.shape}")

# Train models
models = {
    'LogisticRegression': LogisticRegression(random_state=42),
    'RandomForest': RandomForestClassifier(n_estimators=100, random_state=42),
    'SVM': SVC(probability=True, random_state=42)
}

model_results = {}
for name, model in models.items():
    print(f"\nTraining {name}...")
    
    # Cross-validation
    cv_scores = cross_val_score(model, X, y, cv=5, scoring='roc_auc')
    print(f"  CV AUC: {cv_scores.mean():.3f} (+/- {cv_scores.std()*2:.3f})")
    
    # Train on full data
    model.fit(X_train, y_train)
    
    # Test performance
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    
    auc = roc_auc_score(y_test, y_prob)
    print(f"  Test AUC: {auc:.3f}")
    
    model_results[name] = {
        'model': model,
        'cv_auc': cv_scores.mean(),
        'test_auc': auc,
        'y_prob': y_prob
    }

# Best model
best_model_name = max(model_results, key=lambda x: model_results[x]['test_auc'])
best_model = model_results[best_model_name]

print(f"\n*** Best Model: {best_model_name} (Test AUC={best_model['test_auc']:.3f}) ***")

# Save results
results_df.to_csv(f'{RESULTS_DIR}/external/gse35452_results.csv')
print(f"\nSaved: gse35452_results.csv")

print("\n" + "=" * 60)
print("Analysis Complete!")
print("=" * 60)
