#!/usr/bin/env python3
"""
FerroScore-IO: Immuno-Oncology Response Predictor
Pan-cancer ferroptosis signature for predicting immune checkpoint inhibitor response
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

PROC_DIR = '../data/processed'
RESULTS_DIR = '../results'
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)

def load_gene_sets():
    """Load Ferro-IO gene sets - optimized for immunotherapy prediction"""
    gene_sets = {
        # Core ferroptosis genes
        'Ferroptosis Driver': [
            'ACSL4', 'LPCAT3', 'ALOX15', 'ALOX5', 'NOX1', 'NOX4', 
            'P53', 'SAT1', 'CARS1', 'ACSL3'
        ],
        'Ferroptosis Suppressor': [
            'GPX4', 'SLC7A11', 'SLC3A2', 'NFE2L2', 'HMOX1', 
            'FTH1', 'FTL', 'SOD1', 'SOD2', 'GCLC', 'GCLM', 'MT1G'
        ],
        # Immune-related genes for IO prediction
        'Immune Checkpoints': [
            'CD274', 'PDCD1', 'CTLA4', 'PDCD1LG2', 'LAG3', 
            'TIGIT', 'HAVCR2', 'IDO1'
        ],
        'Antigen Presentation': [
            'HLA-A', 'HLA-B', 'HLA-C', 'B2M', 'TAP1', 'TAP2'
        ],
        'Immune Infiltration': [
            'CD8A', 'CD8B', 'CD4', 'FOXP3', 'CD68', 'CD163'
        ],
        # Iron metabolism
        'Iron Metabolism': [
            'TFRC', 'SLC40A1', 'FTH1', 'FTL', 'HEPC', 'IRP1', 'IRP2'
        ]
    }
    return gene_sets

def calculate_ferroscore(expr_matrix, driver_genes, suppressor_genes):
    """Calculate FerroScore based on ferroptosis potential"""
    print("Calculating FerroScore...")
    
    n_samples = expr_matrix.shape[1]
    raw_scores = np.zeros(n_samples)
    
    for i, sample in enumerate(expr_matrix.columns):
        sample_expr = expr_matrix[sample]
        
        # Get ranks
        ranks_desc = sample_expr.rank(ascending=False, method='min')
        n_genes = len(ranks_desc)
        
        # Driver Score
        driver_in_data = [g for g in driver_genes if g in ranks_desc.index]
        if driver_in_data:
            driver_ranks = ranks_desc[driver_in_data]
            driver_scores = (n_genes - driver_ranks + 1) / n_genes
            driver_score = driver_scores.mean()
        else:
            driver_score = 0.5
        
        # Suppressor Score (inverse)
        suppressor_in_data = [g for g in suppressor_genes if g in ranks_desc.index]
        if suppressor_in_data:
            suppressor_ranks_asc = sample_expr[suppressor_in_data].rank(ascending=True, method='min')
            suppressor_scores = (n_genes - suppressor_ranks_asc + 1) / n_genes
            suppressor_score = suppressor_scores.mean()
        else:
            suppressor_score = 0.5
        
        raw_scores[i] = driver_score + suppressor_score
        
        if i % 500 == 0:
            print(f"    Processed {i}/{n_samples} samples")
    
    # Normalize
    scaler = MinMaxScaler()
    ferroscore = scaler.fit_transform(raw_scores.reshape(-1, 1)).flatten()
    
    result = pd.Series(ferroscore, index=expr_matrix.columns, name='FerroScore')
    print(f"  FerroScore range: [{result.min():.3f}, {result.max():.3f}]")
    
    return result

def calculate_immune_score(expr_matrix, immune_genes):
    """Calculate Immune Score for IO prediction"""
    print("\nCalculating Immune Score...")
    
    # Use ssGSEA-like approach or simple mean
    immune_in_data = [g for g in immune_genes if g in expr_matrix.index]
    print(f"  Immune genes available: {len(immune_in_data)}/{len(immune_genes)}")
    
    if len(immune_in_data) > 0:
        immune_expr = expr_matrix.loc[immune_in_data].mean()
    else:
        immune_expr = pd.Series(0.5, index=expr_matrix.columns)
    
    # Normalize
    scaler = MinMaxScaler()
    immune_score = scaler.fit_transform(immune_expr.values.reshape(-1, 1)).flatten()
    immune_score = pd.Series(immune_score, index=expr_matrix.columns, name='Immune_Score')
    
    return immune_score

def calculate_ferro_immu_score(ferroscore, immune_score, weight_immune=0.4):
    """
    Combine FerroScore and Immune Score for IO prediction
    
    Logic:
    - High FerroScore = sensitive to ferroptosis
    - High Immune Score = hot tumor, better IO response
    - Combined = Ferro-Immu Score
    """
    print("\nCalculating Ferro-Immu Combined Score...")
    
    # Weighted combination
    combined = (1 - weight_immune) * ferroscore + weight_immune * immune_score
    
    # Normalize
    scaler = MinMaxScaler()
    combined_scaled = scaler.fit_transform(combined.values.reshape(-1, 1)).flatten()
    
    result = pd.Series(combined_scaled, index=ferroscore.index, name='FerroImmu_Score')
    
    print(f"  Weight: FerroScore={1-weight_immune}, Immune={weight_immune}")
    print(f"  FerroImmu Score range: [{result.min():.3f}, {result.max():.3f}]")
    
    return result

def predict_io_response(ferro_immu_score):
    """Predict IO response category"""
    if ferro_immu_score >= 0.6:
        return 'Likely Responder', '#2ecc71'
    elif ferro_immu_score >= 0.4:
        return 'Intermediate', '#f39c12'
    else:
        return 'Likely Non-responder', '#e74c3c'

def visualize_scores(scores_df, cancer_types=None):
    """Visualize score distributions"""
    print("\nGenerating visualizations...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Distribution of scores
    axes[0, 0].hist(scores_df['FerroScore'], bins=50, alpha=0.7, label='FerroScore', color='red')
    axes[0, 0].set_xlabel('FerroScore')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('FerroScore Distribution')
    
    axes[0, 1].hist(scores_df['Immune_Score'], bins=50, alpha=0.7, label='Immune_Score', color='blue')
    axes[0, 1].set_xlabel('Immune Score')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Immune Score Distribution')
    
    axes[1, 0].hist(scores_df['FerroImmu_Score'], bins=50, alpha=0.7, label='FerroImmu', color='green')
    axes[1, 0].set_xlabel('FerroImmu Score')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('FerroImmu Score Distribution')
    
    # Correlation plot
    axes[1, 1].scatter(scores_df['FerroScore'], scores_df['Immune_Score'], alpha=0.3)
    axes[1, 1].set_xlabel('FerroScore')
    axes[1, 1].set_ylabel('Immune Score')
    axes[1, 1].set_title('FerroScore vs Immune Score')
    
    # Add correlation
    corr = scores_df['FerroScore'].corr(scores_df['Immune_Score'])
    axes[1, 1].text(0.05, 0.95, f'r={corr:.3f}', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/figures/ferro_immu_distributions.png', dpi=300)
    print(f"  Saved: ferro_immu_distributions.png")
    plt.close()

def main():
    """Main algorithm function"""
    print("=" * 60)
    print("FerroScore-IO: Immuno-Oncology Algorithm")
    print("=" * 60)
    
    # 1. Load gene sets
    gene_sets = load_gene_sets()
    driver_genes = gene_sets.get('Ferroptosis Driver', [])
    suppressor_genes = gene_sets.get('Ferroptosis Suppressor', [])
    immune_genes = (gene_sets.get('Immune Checkpoints', []) + 
                   gene_sets.get('Antigen Presentation', []) +
                   gene_sets.get('Immune Infiltration', []))
    
    print(f"\nGene sets:")
    print(f"  Driver genes: {len(driver_genes)}")
    print(f"  Suppressor genes: {len(suppressor_genes)}")
    print(f"  Immune genes: {len(immune_genes)}")
    
    # 2. Load expression data
    expr_file = f'{PROC_DIR}/tcga_ferro_io_expression.csv'
    if not os.path.exists(expr_file):
        print(f"\nExpression file not found: {expr_file}")
        print("  Please run preprocessing first")
        return
    
    expr_matrix = pd.read_csv(expr_file, index_col=0)
    print(f"\nLoaded expression: {expr_matrix.shape}")
    
    # 3. Calculate FerroScore
    ferroscore = calculate_ferroscore(expr_matrix, driver_genes, suppressor_genes)
    
    # 4. Calculate Immune Score
    immune_score = calculate_immune_score(expr_matrix, immune_genes)
    
    # 5. Calculate combined score
    ferro_immu_score = calculate_ferro_immu_score(ferroscore, immune_score, weight_immune=0.4)
    
    # 6. Combine results
    scores_df = pd.DataFrame({
        'FerroScore': ferroscore,
        'Immune_Score': immune_score,
        'FerroImmu_Score': ferro_immu_score
    })
    
    # 7. Load clinical data
    clin_file = f'{PROC_DIR}/tcga_clinical.csv'
    if os.path.exists(clin_file):
        clinical = pd.read_csv(clin_file)
        clinical_dict = dict(zip(clinical['sample_id'], clinical['cancer_type']))
        scores_df['cancer_type'] = scores_df.index.map(clinical_dict)
    
    # 8. Save results
    scores_df.to_csv(f'{RESULTS_DIR}/tables/ferro_immu_scores.csv')
    print(f"\nSaved: ferro_immu_scores.csv")
    
    # 9. Visualize
    visualize_scores(scores_df, scores_df.get('cancer_type'))
    
    print("\n" + "=" * 60)
    print("FerroScore-IO Calculation Complete!")
    print("=" * 60)
    print(f"\nNext step: Run 04_io_model_training.py")

if __name__ == "__main__":
    main()
