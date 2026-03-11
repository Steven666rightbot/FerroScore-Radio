#!/usr/bin/env python3
"""
Step 4: Immuno-Oncology Model Training
Train ML models to predict immunotherapy response
"""

import os
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
import joblib
import warnings
warnings.filterwarnings('ignore')

PROC_DIR = '../data/processed'
RESULTS_DIR = '../results'
os.makedirs(f'{RESULTS_DIR}/models', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)

def load_io_data():
    """Load Ferro-IO scores and immunotherapy response labels"""
    print("Loading Ferro-IO data...")
    
    # Load scores
    scores_file = f'{RESULTS_DIR}/tables/ferro_immu_scores.csv'
    if not os.path.exists(scores_file):
        print(f"Scores file not found: {scores_file}")
        return None, None
    
    scores_df = pd.read_csv(scores_file, index_col=0)
    print(f"  Loaded scores: {scores_df.shape}")
    
    # For now, create mock IO response labels based on FerroImmu Score
    # In real analysis, this would come from immunotherapy cohorts
    print("\n  Note: Using mock IO response labels for demonstration")
    print("  In real analysis, load from GSE78220, IMvigor210, etc.")
    
    # Mock: High FerroImmu Score -> Responder
    median_score = scores_df['FerroImmu_Score'].median()
    scores_df['IO_Response'] = (scores_df['FerroImmu_Score'] >= median_score).astype(int)
    
    responder_count = (scores_df['IO_Response'] == 1).sum()
    non_responder_count = (scores_df['IO_Response'] == 0).sum()
    
    print(f"  Responders: {responder_count}")
    print(f"  Non-responders: {non_responder_count}")
    
    return scores_df, scores_df['IO_Response']

def train_io_models(X_train, X_test, y_train, y_test):
    """Train and evaluate IO prediction models"""
    print("\n" + "=" * 60)
    print("Training IO Prediction Models")
    print("=" * 60)
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Define models
    models = {
        'LogisticRegression': LogisticRegression(max_iter=1000, random_state=42),
        'RandomForest': RandomForestClassifier(n_estimators=100, random_state=42),
        'GradientBoosting': GradientBoostingClassifier(random_state=42),
        'SVM': SVC(probability=True, random_state=42),
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
            print(f"  Error: {e}")
            continue
    
    # Save scaler
    joblib.dump(scaler, f'{RESULTS_DIR}/models/io_scaler.pkl')
    
    return results, trained_models, scaler

def plot_io_results(results_dict, trained_models, X_test, y_test, scaler):
    """Plot IO model performance"""
    print("\nGenerating IO model plots...")
    
    # 1. Model comparison
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    metrics = ['AUC', 'Accuracy', 'Precision', 'Recall', 'F1']
    results_df = pd.DataFrame(results_dict).T
    
    for i, metric in enumerate(metrics):
        ax = axes[i // 3, i % 3]
        results_df[metric].sort_values(ascending=True).plot(kind='barh', ax=ax)
        ax.set_xlabel(metric)
        ax.set_title(f'{metric} Comparison')
        ax.set_xlim(0, 1)
    
    axes[1, 2].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/figures/io_model_comparison.png', dpi=300)
    print("  Saved: io_model_comparison.png")
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
    ax.set_title('ROC Curves - IO Response Prediction')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    plt.savefig(f'{RESULTS_DIR}/figures/io_roc_curves.png', dpi=300)
    print("  Saved: io_roc_curves.png")
    plt.close()

def save_best_io_model(trained_models, results_dict, X_train, y_train, scaler):
    """Save best IO model"""
    print("\n" + "=" * 60)
    print("Saving Best IO Model")
    print("=" * 60)
    
    best_model_name = max(results_dict, key=lambda x: results_dict[x]['AUC'])
    best_model = trained_models[best_model_name]
    best_auc = results_dict[best_model_name]['AUC']
    
    print(f"Best model: {best_model_name} (AUC={best_auc:.3f})")
    
    # Retrain on full dataset
    X_train_scaled = scaler.fit_transform(X_train)
    best_model.fit(X_train_scaled, y_train)
    
    # Save model
    model_path = f'{RESULTS_DIR}/models/best_io_model.pkl'
    joblib.dump({
        'model': best_model,
        'model_name': best_model_name,
        'scaler': scaler,
        'features': ['FerroScore', 'Immune_Score', 'FerroImmu_Score']
    }, model_path)
    print(f"  Saved: {model_path}")
    
    # Save results
    results_df = pd.DataFrame(results_dict).T
    results_df.to_csv(f'{RESULTS_DIR}/tables/io_model_performance.csv')
    print(f"  Saved: io_model_performance.csv")
    
    return best_model_name, best_auc

def main():
    """Main IO model training"""
    print("=" * 60)
    print("FerroScore-IO: Model Training")
    print("=" * 60)
    
    # 1. Load data
    scores_df, labels = load_io_data()
    if scores_df is None:
        return
    
    # 2. Prepare features
    feature_cols = ['FerroScore', 'Immune_Score', 'FerroImmu_Score']
    features = scores_df[feature_cols]
    
    print(f"\nFeatures: {features.shape}")
    print(f"Labels: {labels.value_counts().to_dict()}")
    
    # 3. Split data
    X_train, X_test, y_train, y_test = train_test_split(
        features, labels, test_size=0.3, random_state=42, stratify=labels
    )
    print(f"\nTrain: {X_train.shape}, Test: {X_test.shape}")
    
    # 4. Train models
    results, trained_models, scaler = train_io_models(X_train, X_test, y_train, y_test)
    
    # 5. Plot results
    plot_io_results(results, trained_models, X_test, y_test, scaler)
    
    # 6. Save best model
    best_name, best_auc = save_best_io_model(
        trained_models, results, 
        pd.concat([X_train, X_test]), 
        pd.concat([y_train, y_test]),
        scaler
    )
    
    print("\n" + "=" * 60)
    print("IO Model Training Complete!")
    print("=" * 60)
    print(f"\nBest Model: {best_name} (AUC={best_auc:.3f})")
    print("\nNext step: Run 05_io_validation.py")

if __name__ == "__main__":
    main()
