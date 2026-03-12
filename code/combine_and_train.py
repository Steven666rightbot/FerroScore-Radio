#!/usr/bin/env python3
"""
Extract GSE91061 response labels and combine with GSE78220
Train final model on combined immunotherapy data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix
import xgboost as xgb
import lightgbm as lgb
import joblib
import os

RESULTS_DIR = '../results'
EXTERNAL_DIR = '../results/external'

print("=" * 60)
print("Combining GSE78220 + GSE91061 and Training Final Model")
print("=" * 60)

# Load GSE78220
print("\n1. Loading GSE78220...")
scores_78220 = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_ferro_immuno_scores.csv', index_col=0)
response_78220 = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_response.csv')

import gzip
sample_metadata = {}
with gzip.open('../data/external/GSE78220_series_matrix_full.txt.gz', 'rt', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('!Sample_title'):
            parts = line.strip().split('\t')
            titles = [p.strip('"') for p in parts[1:]]
            sample_metadata['title'] = titles
        elif line.startswith('!Sample_geo_accession'):
            parts = line.strip().split('\t')
            gsm_ids = [p.strip('"') for p in parts[1:]]
            sample_metadata['gsm'] = gsm_ids

gsm_to_patient = dict(zip(sample_metadata['gsm'], sample_metadata['title']))
response_78220['patient_id'] = response_78220['sample_id'].map(gsm_to_patient)
scores_78220['patient_id'] = scores_78220.index.str.replace('.baseline', '')

data_78220 = scores_78220.merge(response_78220[['patient_id', 'response']], on='patient_id', how='inner')
data_78220['dataset'] = 'GSE78220'
data_78220['sample_id'] = data_78220['patient_id']
print(f"   GSE78220: {len(data_78220)} samples")

# Load GSE91061 phenotype to get response labels
print("\n2. Loading GSE91061 phenotype...")
pheno_91061 = pd.read_csv(f'../data/external/GSE91061_phenotype.csv')
print(f"   Phenotype shape: {pheno_91061.shape}")
print(f"   Columns: {pheno_91061.columns.tolist()}")
print(f"   First few rows:")
print(pheno_91061.head())

# Extract response from 'response' column
# Response values: PD (progressive disease), SD (stable disease), PRCR (partial/complete response), UNK (unknown)
def map_response(resp):
    if pd.isna(resp):
        return -1
    resp = str(resp).upper()
    if resp in ['PRCR', 'CR', 'PR']:  # Responder
        return 1
    elif resp in ['PD', 'SD']:  # Non-responder
        return 0
    else:
        return -1

pheno_91061['response_binary'] = pheno_91061['response'].apply(map_response)

# Filter to pre-treatment samples only and known response
pheno_91061_pre = pheno_91061[pheno_91061['visit (pre or on treatment)'] == 'Pre'].copy()
pheno_91061_pre = pheno_91061_pre[pheno_91061_pre['response_binary'] != -1]

print(f"\n   Pre-treatment samples with known response: {len(pheno_91061_pre)}")
print(f"   Responders: {(pheno_91061_pre['response_binary']==1).sum()}")
print(f"   Non-responders: {(pheno_91061_pre['response_binary']==0).sum()}")

# Load GSE91061 scores
scores_91061 = pd.read_csv(f'{EXTERNAL_DIR}/gse91061_ferro_immuno_scores.csv', index_col=0)
print(f"\n   GSE91061 scores: {len(scores_91061)} samples")

# Match sample IDs
# Phenotype has sample_id like 'GSM2420259', scores have column names like 'Pt1_Pre_AD101148-6'
# Need to extract patient ID from phenotype sample_id

# For now, let's use a simplified approach - assume order matches or use patient number
# Extract patient number from sample titles in series matrix
print("\n3. Matching GSE91061 samples...")

# Load series matrix to get sample titles
with gzip.open('../data/external/GSE91061_series_matrix_full.txt.gz', 'rt', encoding='utf-8', errors='ignore') as f:
    gsm_list = []
    title_list = []
    for line in f:
        if line.startswith('!Sample_geo_accession'):
            gsm_list = [p.strip('"') for p in line.strip().split('\t')[1:]]
        elif line.startswith('!Sample_title'):
            title_list = [p.strip('"') for p in line.strip().split('\t')[1:]]

gsm_to_title = dict(zip(gsm_list, title_list))

# Use title as the matching key (it matches the score index)
pheno_91061_pre['sample_title'] = pheno_91061_pre['sample_id'].map(gsm_to_title)
pheno_91061_pre = pheno_91061_pre.rename(columns={'sample_title': 'sample_key'})

# Score index is like 'Pt1_Pre_AD101148-6'
scores_91061['sample_key'] = scores_91061.index

# Merge on sample_key
data_91061 = scores_91061.merge(pheno_91061_pre[['sample_key', 'response_binary']], on='sample_key', how='inner')
data_91061 = data_91061.rename(columns={'response_binary': 'response'})
data_91061['dataset'] = 'GSE91061'
data_91061['sample_id'] = data_91061['sample_key']
print(f"\n   GSE91061 matched: {len(data_91061)} samples")

# Combine datasets
print("\n4. Combining datasets...")
combined = pd.concat([
    data_78220[['sample_id', 'FerroScore', 'Immune_Score', 'FerroImmuno_Score', 'response', 'dataset']],
    data_91061[['sample_id', 'FerroScore', 'Immune_Score', 'FerroImmuno_Score', 'response', 'dataset']]
], ignore_index=True)

print(f"\n   Combined: {len(combined)} samples")
print(f"   GSE78220: {(combined['dataset']=='GSE78220').sum()}")
print(f"   GSE91061: {(combined['dataset']=='GSE91061').sum()}")
print(f"   Responders: {(combined['response']==1).sum()}")
print(f"   Non-responders: {(combined['response']==0).sum()}")

# Prepare features
features = ['FerroScore', 'Immune_Score', 'FerroImmuno_Score']
X = combined[features].values
y = combined['response'].values

print(f"\n   Features: {features}")
print(f"   X shape: {X.shape}")

# Train model
print("\n5. Training model on combined data...")

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

models = {
    'RandomForest': RandomForestClassifier(n_estimators=100, random_state=42),
    'GradientBoosting': GradientBoostingClassifier(random_state=42),
    'LogisticRegression': LogisticRegression(max_iter=1000, random_state=42),
}

results = {}

for name, model in models.items():
    print(f"\n   Training {name}...")
    
    # LOO CV
    loo = LeaveOneOut()
    y_true_list = []
    y_prob_list = []
    
    for train_idx, test_idx in loo.split(X):
        X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        model.fit(X_train, y_train)
        if hasattr(model, 'predict_proba'):
            y_prob = model.predict_proba(X_test)[:, 1][0]
        else:
            y_prob = model.decision_function(X_test)[0]
        y_pred = model.predict(X_test)[0]
        
        y_true_list.append(y_test[0])
        y_prob_list.append(y_prob)
    
    auc = roc_auc_score(y_true_list, y_prob_list)
    acc = accuracy_score(y_true_list, [1 if p > 0.5 else 0 for p in y_prob_list])
    
    results[name] = {'AUC': auc, 'Accuracy': acc}
    print(f"      LOO AUC: {auc:.3f}, Accuracy: {acc:.3f}")

# Best model
best_name = max(results, key=lambda x: results[x]['AUC'])
print(f"\n   Best model: {best_name} (AUC={results[best_name]['AUC']:.3f})")

# Train final model on all data
print("\n6. Training final model on all data...")
final_model = models[best_name]
final_model.fit(X_scaled, y)

# Save
model_path = f'{RESULTS_DIR}/models/best_model_combined.pkl'
joblib.dump({
    'model': final_model,
    'model_name': best_name,
    'scaler': scaler,
    'features': features,
    'datasets': ['GSE78220', 'GSE91061'],
    'n_samples': len(combined)
}, model_path)
print(f"   [OK] Saved: {model_path}")

# Save combined data
combined.to_csv(f'{RESULTS_DIR}/tables/combined_immunotherapy_data.csv', index=False)
print(f"   [OK] Saved: combined_immunotherapy_data.csv")

# Visualization
print("\n7. Generating visualization...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Model comparison
model_names = list(results.keys())
aucs = [results[m]['AUC'] for m in model_names]
axes[0].barh(model_names, aucs, color='steelblue')
axes[0].set_xlabel('LOO AUC')
axes[0].set_title('Model Performance (Combined Data)')
axes[0].axvline(x=0.5, color='red', linestyle='--')

# Score distribution
responder_scores = combined[combined['response']==1]['FerroImmuno_Score']
non_responder_scores = combined[combined['response']==0]['FerroImmuno_Score']
axes[1].hist(non_responder_scores, bins=10, alpha=0.6, label='Non-responder', color='red')
axes[1].hist(responder_scores, bins=10, alpha=0.6, label='Responder', color='green')
axes[1].set_xlabel('FerroImmuno Score')
axes[1].set_ylabel('Count')
axes[1].set_title('Score Distribution')
axes[1].legend()

plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/figures/combined_model_results.png', dpi=300)
print(f"   [OK] Saved: combined_model_results.png")

# Statistical test
from scipy import stats
t_stat, p_value = stats.ttest_ind(responder_scores, non_responder_scores)

print("\n" + "=" * 60)
print("COMPLETE!")
print("=" * 60)
print(f"""
Final Model Summary:
- Datasets: GSE78220 + GSE91061
- Total samples: {len(combined)}
  - GSE78220: {(combined['dataset']=='GSE78220').sum()}
  - GSE91061: {(combined['dataset']=='GSE91061').sum()}
- Responders: {(combined['response']==1).sum()}
- Non-responders: {(combined['response']==0).sum()}
- Best model: {best_name}
- LOO AUC: {results[best_name]['AUC']:.3f}
- LOO Accuracy: {results[best_name]['Accuracy']:.3f}
- T-test p-value: {p_value:.4f}

Model saved and ready for use!
""")
