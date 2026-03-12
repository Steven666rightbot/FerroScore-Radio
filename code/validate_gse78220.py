#!/usr/bin/env python3
"""
Validate FerroScore-Immuno on GSE78220
Compare predictions with actual anti-PD-1 response
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, confusion_matrix
import joblib
import os

RESULTS_DIR = '../results'
EXTERNAL_DIR = '../results/external'
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)

print("=" * 60)
print("GSE78220 Validation - FerroScore-Immuno vs Anti-PD-1 Response")
print("=" * 60)

# Load model
print("\n1. Loading trained model...")
model_data = joblib.load(f'{RESULTS_DIR}/models/best_model.pkl')
model = model_data['model']
scaler = model_data['scaler']
features = model_data['features']
print(f"   Model: {model_data['model_name']}")
print(f"   Features: {features}")

# Load GSE78220 scores
print("\n2. Loading GSE78220 scores...")
scores_df = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_ferro_immuno_scores.csv', index_col=0)
print(f"   Samples: {len(scores_df)}")
print(f"   Score columns: {list(scores_df.columns)}")

# Load response data
print("\n3. Loading response data...")
response_df = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_response.csv')
print(f"   Response data: {len(response_df)}")
print(f"   Responders: {(response_df['response']==1).sum()}")
print(f"   Non-responders: {(response_df['response']==0).sum()}")

# Map sample IDs
print("\n4. Mapping sample IDs...")
# GSE78220 scores have sample names like 'Pt1.baseline', 'Pt2.baseline', etc.
# Response data has GSM IDs, need to map to patient IDs

# Load original metadata to get mapping
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
        elif line.startswith('!Sample_source_name_ch1'):
            parts = line.strip().split('\t')
            sources = [p.strip('"') for p in parts[1:]]
            sample_metadata['source'] = sources

# Create mapping: GSM -> Patient ID -> Response
gsm_to_patient = {}
if 'gsm' in sample_metadata and 'title' in sample_metadata:
    for gsm, title in zip(sample_metadata['gsm'], sample_metadata['title']):
        gsm_to_patient[gsm] = title

print(f"   GSM to patient mapping: {len(gsm_to_patient)}")

# Map response data
response_df['patient_id'] = response_df['sample_id'].map(gsm_to_patient)
print(f"   Mapped responses: {response_df['patient_id'].notna().sum()}")

# Score sample names are like 'Pt1.baseline', need to match to 'Pt1', 'Pt2', etc.
scores_df['patient_id'] = scores_df.index.str.replace('.baseline', '')

# Merge scores with response
merged_df = scores_df.merge(response_df[['patient_id', 'response']], on='patient_id', how='inner')
print(f"\n   Merged data: {len(merged_df)} samples")
print(f"   Responders: {(merged_df['response']==1).sum()}")
print(f"   Non-responders: {(merged_df['response']==0).sum()}")

if len(merged_df) == 0:
    print("\n   [X] No matching samples found!")
    print("   Sample IDs in scores:", scores_df['patient_id'].tolist()[:5])
    print("   Sample IDs in response:", response_df['patient_id'].dropna().tolist()[:5])
    exit(1)

# Prepare features for prediction
print("\n5. Making predictions...")
X = merged_df[features].values
X_scaled = scaler.transform(X)

# Predict probabilities
y_true = merged_df['response'].values
y_prob = model.predict_proba(X_scaled)[:, 1]
y_pred = model.predict(X_scaled)

# Calculate metrics
auc = roc_auc_score(y_true, y_prob)
acc = accuracy_score(y_true, y_pred)

print(f"\n   AUC: {auc:.3f}")
print(f"   Accuracy: {acc:.3f}")

# Confusion matrix
cm = confusion_matrix(y_true, y_pred)
print(f"\n   Confusion Matrix:")
print(f"   TN={cm[0,0]}, FP={cm[0,1]}")
print(f"   FN={cm[1,0]}, TP={cm[1,1]}")

# Save results
results_df = pd.DataFrame({
    'patient_id': merged_df['patient_id'],
    'FerroScore': merged_df['FerroScore'],
    'Immune_Score': merged_df['Immune_Score'],
    'FerroImmuno_Score': merged_df['FerroImmuno_Score'],
    'actual_response': y_true,
    'predicted_probability': y_prob,
    'predicted_response': y_pred
})
results_df.to_csv(f'{EXTERNAL_DIR}/gse78220_validation_results.csv', index=False)
print(f"\n   [OK] Saved: gse78220_validation_results.csv")

# Visualizations
print("\n6. Generating visualizations...")

# ROC curve
fpr, tpr, _ = roc_curve(y_true, y_prob)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

axes[0].plot(fpr, tpr, label=f'FerroScore-Immuno (AUC={auc:.3f})', linewidth=2)
axes[0].plot([0, 1], [0, 1], 'k--', label='Random')
axes[0].set_xlabel('False Positive Rate')
axes[0].set_ylabel('True Positive Rate')
axes[0].set_title('ROC Curve - GSE78220 Validation')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Score distribution by response
responder_scores = merged_df[merged_df['response']==1]['FerroImmuno_Score']
non_responder_scores = merged_df[merged_df['response']==0]['FerroImmuno_Score']

axes[1].hist(non_responder_scores, bins=10, alpha=0.6, label='Non-responder', color='red')
axes[1].hist(responder_scores, bins=10, alpha=0.6, label='Responder', color='green')
axes[1].set_xlabel('FerroImmuno Score')
axes[1].set_ylabel('Count')
axes[1].set_title('Score Distribution by Response')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/figures/gse78220_validation.png', dpi=300)
print(f"   [OK] Saved: gse78220_validation.png")
plt.close()

# Statistical test
from scipy import stats
t_stat, p_value = stats.ttest_ind(responder_scores, non_responder_scores)
print(f"\n   T-test: t={t_stat:.3f}, p={p_value:.4f}")

print("\n" + "=" * 60)
print("Validation Complete!")
print("=" * 60)
print(f"""
Summary:
- Dataset: GSE78220 (Melanoma Anti-PD-1)
- Samples: {len(merged_df)}
- AUC: {auc:.3f}
- Accuracy: {acc:.3f}
- P-value: {p_value:.4f}

Interpretation:
[{'OK' if auc > 0.6 else 'NG'}] AUC {'> 0.6 (acceptable)' if auc > 0.6 else '< 0.6 (poor)'}
[{'OK' if p_value < 0.05 else 'NG'}] Significant difference {'(p < 0.05)' if p_value < 0.05 else '(p >= 0.05)'}
""")
