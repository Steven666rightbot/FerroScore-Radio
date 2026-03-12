#!/usr/bin/env python3
"""
Step 3: FerroScore-Immuno Algorithm
Calculate ferroptosis-immunotherapy sensitivity score
Based on SenScoreR methodology
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns

PROC_DIR = '../data/processed'
RESULTS_DIR = '../results'
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)

def load_gene_sets():
    """Load categorized gene sets"""
    gene_file = '../gene_sets/ferro_immuno_genes.txt'
    
    gene_sets = {}
    current_set = None
    
    with open(gene_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('##'):
                # Extract category name (remove ## and number prefix)
                current_set = line.split('(', 1)[0].strip()
                current_set = current_set.replace('## ', '').replace('##', '').strip()
                # Remove leading numbers like "1. "
                if '. ' in current_set:
                    current_set = current_set.split('. ', 1)[1]
                gene_sets[current_set] = []
            elif line and not line.startswith('#') and current_set:
                gene_sets[current_set].append(line)
    
    return gene_sets

def calculate_ferroscore(expr_matrix, driver_genes, suppressor_genes):
    """
    Calculate FerroScore based on SenScoreR algorithm
    
    Parameters:
    -----------
    expr_matrix : DataFrame (genes x samples)
    driver_genes : list of ferroptosis driver genes
    suppressor_genes : list of ferroptosis suppressor genes
    
    Returns:
    --------
    ferroscore : Series of scores for each sample
    """
    print("Calculating FerroScore...")
    print(f"  Expression matrix: {expr_matrix.shape}")
    
    n_samples = expr_matrix.shape[1]
    raw_scores = np.zeros(n_samples)
    
    for i, sample in enumerate(expr_matrix.columns):
        sample_expr = expr_matrix[sample]
        
        # Get ranks (descending order: highest expression = rank 1)
        ranks_desc = sample_expr.rank(ascending=False, method='min')
        n_genes = len(ranks_desc)
        
        # Calculate Driver Score (high driver expression = high score)
        driver_in_data = [g for g in driver_genes if g in ranks_desc.index]
        if driver_in_data:
            driver_ranks = ranks_desc[driver_in_data]
            # Score: (N - rank + 1) / N, then mean
            driver_scores = (n_genes - driver_ranks + 1) / n_genes
            driver_score = driver_scores.mean()
        else:
            driver_score = 0
        
        # Calculate Suppressor Score (low suppressor expression = high score)
        suppressor_in_data = [g for g in suppressor_genes if g in ranks_desc.index]
        if suppressor_in_data:
            suppressor_ranks = ranks_desc[suppressor_in_data]
            # For suppressors: low expression (high rank number) should give high score
            # So we use ascending ranks directly
            suppressor_ranks_asc = sample_expr[suppressor_in_data].rank(ascending=True, method='min')
            suppressor_scores = (n_genes - suppressor_ranks_asc + 1) / n_genes
            suppressor_score = suppressor_scores.mean()
        else:
            suppressor_score = 0
        
        # Combined raw score
        raw_scores[i] = driver_score + suppressor_score
        
        if i % 500 == 0:
            print(f"    Processed {i}/{n_samples} samples")
    
    # Min-max normalization to [0, 1]
    scaler = MinMaxScaler()
    ferroscore = scaler.fit_transform(raw_scores.reshape(-1, 1)).flatten()
    
    result = pd.Series(ferroscore, index=expr_matrix.columns, name='FerroScore')
    
    print(f"  FerroScore range: [{result.min():.3f}, {result.max():.3f}]")
    print(f"  FerroScore mean: {result.mean():.3f}")
    
    return result

def calculate_immune_score(expr_matrix, immune_genes):
    """
    Calculate Immune Infiltration score
    High immune score = potentially responsive to immunotherapy
    """
    print("\nCalculating Immune Score...")
    
    immune_in_data = [g for g in immune_genes if g in expr_matrix.index]
    print(f"  Immune genes available: {len(immune_in_data)}/{len(immune_genes)}")
    
    # Use ssGSEA-like approach or simple mean
    immune_expr = expr_matrix.loc[immune_in_data].mean()
    
    # Normalize to [0, 1]
    scaler = MinMaxScaler()
    immune_score = scaler.fit_transform(immune_expr.values.reshape(-1, 1)).flatten()
    immune_score = pd.Series(immune_score, index=expr_matrix.columns, name='Immune_Score')
    
    return immune_score

def calculate_ferro_immuno_score(ferroscore, immune_score, weight_immune=0.4):
    """
    Combine FerroScore and Immune Score
    
    Logic:
    - High FerroScore = sensitive to ferroptosis
    - High Immune Score = responsive to immunotherapy
    - Combined: High FerroScore + High Immune = High ImmunoSensitivity
    """
    print("\nCalculating Ferro-Immuno Combined Score...")
    
    # Weighted combination (both high = good for immuno)
    combined = (1 - weight_immune) * ferroscore + weight_immune * immune_score
    
    # Normalize
    scaler = MinMaxScaler()
    combined_scaled = scaler.fit_transform(combined.values.reshape(-1, 1)).flatten()
    
    result = pd.Series(combined_scaled, index=ferroscore.index, name='FerroImmuno_Score')
    
    print(f"  Weight: FerroScore={1-weight_immune}, Immune={weight_immune}")
    print(f"  FerroImmuno Score range: [{result.min():.3f}, {result.max():.3f}]")
    
    return result

def visualize_scores(scores_df, cancer_types=None):
    """Visualize score distributions"""
    print("\nGenerating visualizations...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Distribution of scores
    axes[0, 0].hist(scores_df['FerroScore'], bins=50, alpha=0.7, label='FerroScore')
    axes[0, 0].set_xlabel('FerroScore')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('FerroScore Distribution')
    
    axes[0, 1].hist(scores_df['Immune_Score'], bins=50, alpha=0.7, label='Immune_Score', color='orange')
    axes[0, 1].set_xlabel('Immune Score')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Immune Score Distribution')
    
    axes[1, 0].hist(scores_df['FerroImmuno_Score'], bins=50, alpha=0.7, label='FerroImmuno', color='green')
    axes[1, 0].set_xlabel('FerroImmuno Score')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('FerroImmuno Score Distribution')
    
    # Correlation plot
    axes[1, 1].scatter(scores_df['FerroScore'], scores_df['Immune_Score'], alpha=0.3)
    axes[1, 1].set_xlabel('FerroScore')
    axes[1, 1].set_ylabel('Immune Score')
    axes[1, 1].set_title('FerroScore vs Immune Score')
    
    # Add correlation
    corr = scores_df['FerroScore'].corr(scores_df['Immune_Score'])
    axes[1, 1].text(0.05, 0.95, f'r={corr:.3f}', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/figures/score_distributions.png', dpi=300)
    print(f"  [OK] Saved: score_distributions.png")
    plt.close()
    
    # If cancer types provided, plot by cancer
    if cancer_types is not None:
        scores_df['cancer_type'] = cancer_types
        
        fig, ax = plt.subplots(figsize=(14, 6))
        scores_df.boxplot(column='FerroImmuno_Score', by='cancer_type', ax=ax)
        plt.xticks(rotation=45, ha='right')
        plt.title('FerroImmuno Score by Cancer Type')
        plt.suptitle('')
        plt.tight_layout()
        plt.savefig(f'{RESULTS_DIR}/figures/score_by_cancer.png', dpi=300)
        print(f"  [OK] Saved: score_by_cancer.png")
        plt.close()

def main():
    """Main algorithm function"""
    print("=" * 60)
    print("FerroScore-Immuno: Algorithm Calculation")
    print("=" * 60)
    
    # 1. Load gene sets
    gene_sets = load_gene_sets()
    driver_genes = gene_sets.get('Ferroptosis Driver Genes', [])
    suppressor_genes = gene_sets.get('Ferroptosis Suppressor Genes', [])
    immune_genes = gene_sets.get('T Cell Infiltration Markers', []) + \
                   gene_sets.get('Antigen Presentation Genes', []) + \
                   gene_sets.get('Immune Checkpoint Genes', [])
    
    print(f"\nGene sets:")
    print(f"  Driver genes: {len(driver_genes)}")
    print(f"  Suppressor genes: {len(suppressor_genes)}")
    print(f"  Immune genes: {len(immune_genes)}")
    
    # 2. Load expression data
    expr_file = f'{PROC_DIR}/tcga_ferro_immuno_expression.csv'
    if not os.path.exists(expr_file):
        print(f"\n[X] Expression file not found: {expr_file}")
        print("  Please run 02_data_preprocessing.py first")
        return
    
    expr_matrix = pd.read_csv(expr_file, index_col=0)
    print(f"\nLoaded expression: {expr_matrix.shape}")
    
    # 3. Calculate FerroScore
    ferroscore = calculate_ferroscore(expr_matrix, driver_genes, suppressor_genes)
    
    # 4. Calculate Immune Score
    immune_score = calculate_immune_score(expr_matrix, immune_genes)
    
    # 5. Calculate combined score
    ferro_immuno_score = calculate_ferro_immuno_score(ferroscore, immune_score, weight_immune=0.4)
    
    # 6. Combine results
    scores_df = pd.DataFrame({
        'FerroScore': ferroscore,
        'Immune_Score': immune_score,
        'FerroImmuno_Score': ferro_immuno_score
    })
    
    # 7. Load clinical data for cancer type annotation
    clin_file = f'{PROC_DIR}/tcga_clinical.csv'
    if os.path.exists(clin_file):
        clinical = pd.read_csv(clin_file)
        # Map sample IDs
        clinical_dict = dict(zip(clinical['sample_id'], clinical['cancer_type']))
        scores_df['cancer_type'] = scores_df.index.map(clinical_dict)
    
    # 8. Save results
    scores_df.to_csv(f'{RESULTS_DIR}/tables/ferro_immuno_scores.csv')
    print(f"\n[OK] Saved: ferro_immuno_scores.csv")
    
    # 9. Visualize
    visualize_scores(scores_df, scores_df.get('cancer_type'))
    
    print("\n" + "=" * 60)
    print("FerroScore Calculation Complete!")
    print("=" * 60)
    print(f"\nResults saved to: {RESULTS_DIR}/tables/ferro_immuno_scores.csv")
    print("\nNext step: Run 04_model_training.py")

if __name__ == "__main__":
    main()
