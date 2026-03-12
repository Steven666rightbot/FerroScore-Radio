#!/usr/bin/env python3
"""
Retrain FerroScore-Immuno model using GSE78220 as training set
This is the real immunotherapy data (melanoma anti-PD-1)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, precision_score, recall_score, f1_score
import xgboost as xgb
import lightgbm as lgb
import joblib
import os

RESULTS_DIR = '../results'
EXTERNAL_DIR = '../results/external'
os.makedirs(f'{RESULTS_DIR}/models', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)

print("=" * 60)
print("Retraining FerroScore-Immuno Model on GSE78220")
print("=" * 60)

# Load GSE78220 data
print("\n1. Loading GSE78220 data...")
scores_df = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_ferro_immuno_scores.csv', index_col=0)
response_df = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_response.csv')

print(f"   Scores: {len(scores_df)} samples")
print(f"   Responses: {len(response_df)} samples")

# Map sample IDs
scores_df['patient_id'] = scores_df.index.str.replace('.baseline', '')
response_df['patient_id'] = response_df['sample_id'].str.replace('GSM', 'Pt')  # Simplified mapping

# Actually, let's use the correct mapping from before
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

gsm_to_patient = {}
if 'gsm' in sample_metadata and 'title' in sample_metadata:
    for gsm, title in zip(sample_metadata['gsm'], sample_metadata['title']):
        gsm_to_patient[gsm] = title

response_df['patient_id'] = response_df['sample_id'].map(gsm_to_patient)

# Merge data
merged_df = scores_df.merge(response_df[['patient_id', 'response']], on='patient_id', how='inner')
print(f"\n   Merged data: {len(merged_df)} samples")
print(f"   Responders: {(merged_df['response']==1).sum()}")
print(f"   Non-responders: {(merged_df['response']==0).sum()}")

# Prepare features
features = ['FerroScore', 'Immune_Score', 'FerroImmuno_Score']
X = merged_df[features].values
y = merged_df['response'].values

print(f"\n   Features: {features}")
print(f"   X shape: {X.shape}")
print(f"   y distribution: {np.bincount(y)}")

# Split data
print("\n2. Splitting data...")
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42, stratify=y
)
print(f"   Train: {len(X_train)} samples")
print(f"   Test: {len(X_test)} samples")

# Scale features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train multiple models
print("\n3. Training models...")

models = {
    'LogisticRegression': LogisticRegression(max_iter=1000, random_state=42),
    'RandomForest': RandomForestClassifier(n_estimators=100, random_state=42),
    'GradientBoosting': GradientBoostingClassifier(random_state=42),
    'SVM': SVC(probability=True, random_state=42),
    'XGBoost': xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42),
    'LightGBM': lgb.LGBMClassifier(random_state=42, verbose=-1),
    'MLP': MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=500, random_state=42)
}

results = {}
trained_models = {}

for name, model in models.items():
    print(f"\n   Training {name}...")
    
    try:
        model.fit(X_train_scaled, y_train)
        
        y_pred = model.predict(X_test_scaled)
        if hasattr(model, 'predict_proba'):
            y_prob = model.predict_proba(X_test_scaled)[:, 1]
        else:
            y_prob = model.decision_function(X_test_scaled)
        
        auc = roc_auc_score(y_test, y_prob)
        acc = accuracy_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred, zero_division=0)
        rec = recall_score(y_test, y_pred, zero_division=0)
        f1 = f1_score(y_test, y_pred, zero_division=0)
        
        results[name] = {'AUC': auc, 'Accuracy': acc, 'Precision': prec, 'Recall': rec, 'F1': f1}
        trained_models[name] = model
        
        print(f"      AUC: {auc:.3f}, Acc: {acc:.3f}, F1: {f1:.3f}")
        
    except Exception as e:
        print(f"      [ERR] {e}")

# Cross-validation
print("\n4. Cross-validation (5-fold)...")
cv_results = {}
for name, model in models.items():
    try:
        cv_scores = cross_val_score(model, X_train_scaled, y_train, cv=StratifiedKFold(5), scoring='roc_auc')
        cv_results[name] = cv_scores
        print(f"   {name}: {cv_scores.mean():.3f} (+/- {cv_scores.std()*2:.3f})")
    except:
        pass

# Find best model
best_model_name = max(results, key=lambda x: results[x]['AUC'])
best_model = trained_models[best_model_name]
best_auc = results[best_model_name]['AUC']

print(f"\n   Best model: {best_model_name} (AUC={best_auc:.3f})")

# Save best model
print("\n5. Saving best model...")
model_path = f'{RESULTS_DIR}/models/best_model_gse78220.pkl'
joblib.dump({
    'model': best_model,
    'model_name': best_model_name,
    'scaler': scaler,
    'features': features
}, model_path)
print(f"   [OK] Saved: {model_path}")

# Save results
results_df = pd.DataFrame(results).T
results_df.to_csv(f'{RESULTS_DIR}/tables/model_performance_gse78220.csv')
print(f"   [OK] Saved: model_performance_gse78220.csv")

# Visualizations
print("\n6. Generating visualizations...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Model comparison
model_names = list(results.keys())
aucs = [results[m]['AUC'] for m in model_names]
axes[0].barh(model_names, aucs, color='steelblue')
axes[0].set_xlabel('AUC')
axes[0].set_title('Model Performance on GSE78220')
axes[0].axvline(x=0.5, color='red', linestyle='--', label='Random')
axes[0].legend()

# ROC curve for best model
y_prob_best = best_model.predict_proba(X_test_scaled)[:, 1]
fpr, tpr, _ = roc_curve(y_test, y_prob_best)
axes[1].plot(fpr, tpr, label=f'{best_model_name} (AUC={best_auc:.3f})', linewidth=2)
axes[1].plot([0, 1], [0, 1], 'k--', label='Random')
axes[1].set_xlabel('False Positive Rate')
axes[1].set_ylabel('True Positive Rate')
axes[1].set_title('ROC Curve - Best Model')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/figures/model_training_gse78220.png', dpi=300)
print(f"   [OK] Saved: model_training_gse78220.png")
plt.close()

print("\n" + "=" * 60)
print("Retraining Complete!")
print("=" * 60)
print(f"""
Summary:
- Training set: GSE78220 (Melanoma Anti-PD-1)
- Samples: {len(merged_df)}
- Best model: {best_model_name}
- Test AUC: {best_auc:.3f}
- Test Accuracy: {results[best_model_name]['Accuracy']:.3f}

This model is now trained on real immunotherapy data!
""")
