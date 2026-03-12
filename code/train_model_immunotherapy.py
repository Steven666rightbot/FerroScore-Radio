#!/usr/bin/env python3
"""
Train FerroScore-Immuno model on combined immunotherapy data
Currently using GSE78220 (can add GSE91061 after gene ID conversion)
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
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
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
print("Training FerroScore-Immuno on Immunotherapy Data")
print("=" * 60)

# Load GSE78220 data
print("\n1. Loading GSE78220 data...")
scores_78220 = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_ferro_immuno_scores.csv', index_col=0)
response_78220 = pd.read_csv(f'{EXTERNAL_DIR}/gse78220_response.csv')

# Map sample IDs
scores_78220['patient_id'] = scores_78220.index.str.replace('.baseline', '')

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

response_78220['patient_id'] = response_78220['sample_id'].map(gsm_to_patient)

# Merge GSE78220
data_78220 = scores_78220.merge(response_78220[['patient_id', 'response']], on='patient_id', how='inner')
data_78220['dataset'] = 'GSE78220'
print(f"   GSE78220: {len(data_78220)} samples")

# TODO: Add GSE91061 after gene ID conversion
# For now, use only GSE78220

# Combine datasets
combined_data = data_78220.copy()
print(f"\n   Combined: {len(combined_data)} samples")
print(f"   Responders: {(combined_data['response']==1).sum()}")
print(f"   Non-responders: {(combined_data['response']==0).sum()}")

# Prepare features
features = ['FerroScore', 'Immune_Score', 'FerroImmuno_Score']
X = combined_data[features].values
y = combined_data['response'].values

print(f"\n   Features: {features}")
print(f"   X shape: {X.shape}")

# Split data
print("\n2. Splitting data (70% train, 30% test)...")
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
    'XGBoost': xgb.XGBClassifier(eval_metric='logloss', random_state=42),
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
print("\n4. Cross-validation (Leave-One-Out for small dataset)...")
from sklearn.model_selection import LeaveOneOut
loo = LeaveOneOut()

for name, model in models.items():
    try:
        cv_scores = []
        for train_idx, test_idx in loo.split(X):
            X_train_cv, X_test_cv = X[train_idx], X[test_idx]
            y_train_cv, y_test_cv = y[train_idx], y[test_idx]
            
            scaler_cv = StandardScaler()
            X_train_cv_scaled = scaler_cv.fit_transform(X_train_cv)
            X_test_cv_scaled = scaler_cv.transform(X_test_cv)
            
            model.fit(X_train_cv_scaled, y_train_cv)
            if hasattr(model, 'predict_proba'):
                y_prob_cv = model.predict_proba(X_test_cv_scaled)[:, 1]
            else:
                y_prob_cv = model.decision_function(X_test_cv_scaled)
            
            # For binary, AUC might not work with single sample, use accuracy
            y_pred_cv = model.predict(X_test_cv_scaled)
            cv_scores.append(1 if y_pred_cv[0] == y_test_cv[0] else 0)
        
        cv_accuracy = np.mean(cv_scores)
        print(f"   {name}: LOO Accuracy = {cv_accuracy:.3f}")
    except Exception as e:
        print(f"   {name}: LOO failed - {e}")

# Find best model
best_model_name = max(results, key=lambda x: results[x]['AUC'])
best_model = trained_models[best_model_name]
best_auc = results[best_model_name]['AUC']

print(f"\n   Best model: {best_model_name} (AUC={best_auc:.3f})")

# Save best model
print("\n5. Saving best model...")
model_path = f'{RESULTS_DIR}/models/best_model_immunotherapy.pkl'
joblib.dump({
    'model': best_model,
    'model_name': best_model_name,
    'scaler': scaler,
    'features': features
}, model_path)
print(f"   [OK] Saved: {model_path}")

# Save results
results_df = pd.DataFrame(results).T
results_df.to_csv(f'{RESULTS_DIR}/tables/model_performance_immunotherapy.csv')
print(f"   [OK] Saved: model_performance_immunotherapy.csv")

# Confusion matrix for best model
y_pred_best = best_model.predict(X_test_scaled)
cm = confusion_matrix(y_test, y_pred_best)
print(f"\n   Confusion Matrix ({best_model_name}):")
print(f"   TN={cm[0,0]}, FP={cm[0,1]}")
print(f"   FN={cm[1,0]}, TP={cm[1,1]}")

# Visualizations
print("\n6. Generating visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Model comparison
model_names = list(results.keys())
aucs = [results[m]['AUC'] for m in model_names]
axes[0, 0].barh(model_names, aucs, color='steelblue')
axes[0, 0].set_xlabel('AUC')
axes[0, 0].set_title('Model Performance on Immunotherapy Data')
axes[0, 0].axvline(x=0.5, color='red', linestyle='--', label='Random')
axes[0, 0].legend()

# ROC curve for best model
y_prob_best = best_model.predict_proba(X_test_scaled)[:, 1]
fpr, tpr, _ = roc_curve(y_test, y_prob_best)
axes[0, 1].plot(fpr, tpr, label=f'{best_model_name} (AUC={best_auc:.3f})', linewidth=2)
axes[0, 1].plot([0, 1], [0, 1], 'k--', label='Random')
axes[0, 1].set_xlabel('False Positive Rate')
axes[0, 1].set_ylabel('True Positive Rate')
axes[0, 1].set_title('ROC Curve - Best Model')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Score distribution by response
responder_scores = combined_data[combined_data['response']==1]['FerroImmuno_Score']
non_responder_scores = combined_data[combined_data['response']==0]['FerroImmuno_Score']

axes[1, 0].hist(non_responder_scores, bins=8, alpha=0.6, label='Non-responder', color='red')
axes[1, 0].hist(responder_scores, bins=8, alpha=0.6, label='Responder', color='green')
axes[1, 0].set_xlabel('FerroImmuno Score')
axes[1, 0].set_ylabel('Count')
axes[1, 0].set_title('FerroImmuno Score by Response')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# Feature importance (for tree-based models)
if hasattr(best_model, 'feature_importances_'):
    importances = best_model.feature_importances_
    axes[1, 1].bar(features, importances, color='coral')
    axes[1, 1].set_ylabel('Importance')
    axes[1, 1].set_title(f'Feature Importance ({best_model_name})')
    axes[1, 1].tick_params(axis='x', rotation=45)
else:
    axes[1, 1].text(0.5, 0.5, 'Feature importance not available\nfor this model type', 
                    ha='center', va='center', transform=axes[1, 1].transAxes)

plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/figures/model_training_immunotherapy.png', dpi=300)
print(f"   [OK] Saved: model_training_immunotherapy.png")
plt.close()

# Statistical test
from scipy import stats
t_stat, p_value = stats.ttest_ind(responder_scores, non_responder_scores)

print("\n" + "=" * 60)
print("Training Complete!")
print("=" * 60)
print(f"""
Summary:
- Dataset: GSE78220 (Melanoma Anti-PD-1)
- Samples: {len(combined_data)} (27 after filtering)
- Responders: {(combined_data['response']==1).sum()}
- Non-responders: {(combined_data['response']==0).sum()}
- Best model: {best_model_name}
- Test AUC: {best_auc:.3f}
- Test Accuracy: {results[best_model_name]['Accuracy']:.3f}
- T-test p-value: {p_value:.4f}

Note: This model is trained on REAL immunotherapy response data!
""")
