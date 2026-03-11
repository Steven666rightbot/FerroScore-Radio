#!/usr/bin/env python3
"""
Step 4: Machine Learning Model Training
Train multiple ML models to predict radiotherapy response
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.linear_model import LogisticRegression, RidgeClassifier, ElasticNet
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, precision_score, recall_score, f1_score
from sklearn.metrics import classification_report, confusion_matrix
import xgboost as xgb
import lightgbm as lgb
import joblib
import warnings
warnings.filterwarnings('ignore')

PROC_DIR = '../data/processed'
RESULTS_DIR = '../results'
os.makedirs(f'{RESULTS_DIR}/models', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)

def load_data():
    """Load features and labels"""
    print("Loading data...")
    
    # Load scores
    scores_file = f'{RESULTS_DIR}/tables/ferro_radio_scores.csv'
    if not os.path.exists(scores_file):
        print(f"✗ Scores file not found: {scores_file}")
        return None, None, None
    
    scores_df = pd.read_csv(scores_file, index_col=0)
    print(f"  Loaded scores: {scores_df.shape}")
    
    # Load clinical data
    clin_file = f'{PROC_DIR}/tcga_clinical.csv'
    if not os.path.exists(clin_file):
        print(f"✗ Clinical file not found: {clin_file}")
        return None, None, None
    
    clinical = pd.read_csv(clin_file)
    print(f"  Loaded clinical: {clinical.shape}")
    
    # Load survival data
    surv_file = f'{PROC_DIR}/tcga_survival.csv'
    if os.path.exists(surv_file):
        survival = pd.read_csv(surv_file)
        print(f"  Loaded survival: {survival.shape}")
    else:
        survival = None
    
    return scores_df, clinical, survival

def prepare_labels(clinical_df, survival_df=None):
    """
    Prepare labels for training
    
    Strategy:
    1. Patients with radiotherapy and good outcome (long survival) = Positive
    2. Patients with radiotherapy and poor outcome (short survival) = Negative
    3. Use median survival as cutoff
    """
    print("\nPreparing labels...")
    
    # Filter patients with radiotherapy info
    if 'received_radiotherapy' not in clinical_df.columns:
        print("  ✗ No radiotherapy information")
        return None
    
    rt_patients = clinical_df[clinical_df['received_radiotherapy'] == True].copy()
    print(f"  Patients with radiotherapy: {len(rt_patients)}")
    
    if survival_df is None or len(rt_patients) == 0:
        print("  Using FerroRadio score as proxy label")
        # Use median FerroRadio score as cutoff
        return None
    
    # Merge survival data
    if 'sample' in survival_df.columns:
        rt_patients = rt_patients.merge(
            survival_df[['sample', 'OS', 'OS.time']], 
            left_on='sample_id', 
            right_on='sample',
            how='left'
        )
    
    # Filter patients with survival data
    rt_patients = rt_patients.dropna(subset=['OS.time'])
    print(f"  Patients with survival data: {len(rt_patients)}")
    
    if len(rt_patients) < 50:
        print("  Too few samples, using alternative labeling")
        return None
    
    # Define good vs poor outcome
    median_survival = rt_patients['OS.time'].median()
    rt_patients['label'] = (rt_patients['OS.time'] >= median_survival).astype(int)
    
    print(f"  Median survival: {median_survival:.1f} days")
    print(f"  Good outcome (label=1): {(rt_patients['label']==1).sum()}")
    print(f"  Poor outcome (label=0): {(rt_patients['label']==0).sum()}")
    
    return rt_patients[['sample_id', 'label', 'OS', 'OS.time']]

def prepare_features(scores_df, labels_df=None):
    """Prepare feature matrix"""
    print("\nPreparing features...")
    
    # Use FerroScore, DDR_Score, FerroRadio_Score as features
    feature_cols = ['FerroScore', 'DDR_Score', 'FerroRadio_Score']
    available_cols = [c for c in feature_cols if c in scores_df.columns]
    
    if len(available_cols) == 0:
        print("  ✗ No feature columns found")
        return None, None
    
    features = scores_df[available_cols].copy()
    
    # If no labels, use FerroRadio score to create proxy labels
    if labels_df is None:
        print("  Creating proxy labels from FerroRadio score")
        median_score = features['FerroRadio_Score'].median()
        labels = (features['FerroRadio_Score'] >= median_score).astype(int)
        labels.name = 'label'
    else:
        # Merge with labels
        features = features.reset_index()
        features = features.merge(labels_df[['sample_id', 'label']], 
                                  left_on='index', 
                                  right_on='sample_id',
                                  how='inner')
        features = features.set_index('index')
        labels = features['label']
        features = features[available_cols]
    
    # Remove NaN
    valid_idx = features.dropna().index
    features = features.loc[valid_idx]
    labels = labels.loc[valid_idx]
    
    print(f"  Features shape: {features.shape}")
    print(f"  Labels distribution: {labels.value_counts().to_dict()}")
    
    return features, labels

def train_models(X_train, X_test, y_train, y_test):
    """Train and evaluate multiple ML models"""
    print("\n" + "=" * 60)
    print("Training Machine Learning Models")
    print("=" * 60)
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Define models
    models = {
        'LogisticRegression': LogisticRegression(max_iter=1000, random_state=42),
        'Ridge': RidgeClassifier(random_state=42),
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
        print(f"\nTraining {name}...")
        
        try:
            # Train
            model.fit(X_train_scaled, y_train)
            
            # Predict
            y_pred = model.predict(X_test_scaled)
            y_prob = model.predict_proba(X_test_scaled)[:, 1]
            
            # Evaluate
            auc = roc_auc_score(y_test, y_prob)
            acc = accuracy_score(y_test, y_pred)
            prec = precision_score(y_test, y_pred, zero_division=0)
            rec = recall_score(y_test, y_pred, zero_division=0)
            f1 = f1_score(y_test, y_pred, zero_division=0)
            
            results[name] = {
                'AUC': auc,
                'Accuracy': acc,
                'Precision': prec,
                'Recall': rec,
                'F1': f1
            }
            
            trained_models[name] = model
            
            print(f"  AUC: {auc:.3f}, Acc: {acc:.3f}, F1: {f1:.3f}")
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
            continue
    
    # Save scaler
    joblib.dump(scaler, f'{RESULTS_DIR}/models/scaler.pkl')
    
    return results, trained_models, scaler

def cross_validate_models(X, y, models_dict):
    """Perform cross-validation"""
    print("\n" + "=" * 60)
    print("Cross-Validation (5-fold)")
    print("=" * 60)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    cv_results = {}
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    for name, model in models_dict.items():
        try:
            scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='roc_auc')
            cv_results[name] = {
                'mean_auc': scores.mean(),
                'std_auc': scores.std(),
                'scores': scores
            }
            print(f"{name}: {scores.mean():.3f} (+/- {scores.std()*2:.3f})")
        except Exception as e:
            print(f"{name}: ✗ Error - {e}")
    
    return cv_results

def plot_results(results_dict, trained_models, X_test, y_test, scaler):
    """Plot model performance"""
    print("\nGenerating plots...")
    
    # 1. Model comparison bar plot
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    metrics = ['AUC', 'Accuracy', 'Precision', 'Recall', 'F1']
    results_df = pd.DataFrame(results_dict).T
    
    for i, metric in enumerate(metrics):
        ax = axes[i // 3, i % 3]
        results_df[metric].sort_values(ascending=True).plot(kind='barh', ax=ax)
        ax.set_xlabel(metric)
        ax.set_title(f'{metric} Comparison')
        ax.set_xlim(0, 1)
    
    # Remove empty subplot
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/figures/model_comparison.png', dpi=300)
    print("  ✓ Saved: model_comparison.png")
    plt.close()
    
    # 2. ROC curves
    fig, ax = plt.subplots(figsize=(10, 8))
    
    X_test_scaled = scaler.transform(X_test)
    
    for name, model in trained_models.items():
        y_prob = model.predict_proba(X_test_scaled)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, y_prob)
        auc = roc_auc_score(y_test, y_prob)
        ax.plot(fpr, tpr, label=f'{name} (AUC={auc:.3f})')
    
    ax.plot([0, 1], [0, 1], 'k--', label='Random')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curves')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    plt.savefig(f'{RESULTS_DIR}/figures/roc_curves.png', dpi=300)
    print("  ✓ Saved: roc_curves.png")
    plt.close()

def save_best_model(trained_models, results_dict, X_train, y_train, scaler):
    """Save the best performing model"""
    print("\n" + "=" * 60)
    print("Saving Best Model")
    print("=" * 60)
    
    # Find best model by AUC
    best_model_name = max(results_dict, key=lambda x: results_dict[x]['AUC'])
    best_model = trained_models[best_model_name]
    best_auc = results_dict[best_model_name]['AUC']
    
    print(f"Best model: {best_model_name} (AUC={best_auc:.3f})")
    
    # Retrain on full dataset
    X_train_scaled = scaler.fit_transform(X_train)
    best_model.fit(X_train_scaled, y_train)
    
    # Save model
    model_path = f'{RESULTS_DIR}/models/best_model.pkl'
    joblib.dump({
        'model': best_model,
        'model_name': best_model_name,
        'scaler': scaler,
        'features': ['FerroScore', 'DDR_Score', 'FerroRadio_Score']
    }, model_path)
    print(f"  ✓ Saved: {model_path}")
    
    # Save results
    results_df = pd.DataFrame(results_dict).T
    results_df.to_csv(f'{RESULTS_DIR}/tables/model_performance.csv')
    print(f"  ✓ Saved: model_performance.csv")
    
    return best_model_name, best_auc

def main():
    """Main training function"""
    print("=" * 60)
    print("FerroScore-Radio: Machine Learning Model Training")
    print("=" * 60)
    
    # 1. Load data
    scores_df, clinical, survival = load_data()
    if scores_df is None:
        return
    
    # 2. Prepare labels
    labels_df = prepare_labels(clinical, survival)
    
    # 3. Prepare features
    features, labels = prepare_features(scores_df, labels_df)
    if features is None:
        return
    
    # 4. Split data
    X_train, X_test, y_train, y_test = train_test_split(
        features, labels, test_size=0.3, random_state=42, stratify=labels
    )
    print(f"\nTrain set: {X_train.shape}, Test set: {X_test.shape}")
    
    # 5. Train models
    results, trained_models, scaler = train_models(X_train, X_test, y_train, y_test)
    
    # 6. Cross-validation
    cv_results = cross_validate_models(
        pd.concat([X_train, X_test]), 
        pd.concat([y_train, y_test]), 
        trained_models
    )
    
    # 7. Plot results
    plot_results(results, trained_models, X_test, y_test, scaler)
    
    # 8. Save best model
    best_name, best_auc = save_best_model(
        trained_models, results, 
        pd.concat([X_train, X_test]), 
        pd.concat([y_train, y_test]),
        scaler
    )
    
    print("\n" + "=" * 60)
    print("Model Training Complete!")
    print("=" * 60)
    print(f"\nBest Model: {best_name} (AUC={best_auc:.3f})")
    print(f"\nNext step: Run 05_validation.py")

if __name__ == "__main__":
    main()
