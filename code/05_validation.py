#!/usr/bin/env python3
"""
Step 5: Validation Analysis
Survival analysis, immunotherapy-specific validation, and external validation
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from sklearn.metrics import roc_auc_score, roc_curve
import warnings
warnings.filterwarnings('ignore')

PROC_DIR = '../data/processed'
RESULTS_DIR = '../results'
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)
os.makedirs(f'{RESULTS_DIR}/figures', exist_ok=True)

def load_data():
    """Load all necessary data"""
    print("Loading data...")
    
    # Load scores
    scores_file = f'{RESULTS_DIR}/tables/ferro_radio_scores.csv'
    scores_df = pd.read_csv(scores_file, index_col=0)
    
    # Load clinical
    clin_file = f'{PROC_DIR}/tcga_clinical.csv'
    clinical = pd.read_csv(clin_file)
    
    # Load survival
    surv_file = f'{PROC_DIR}/tcga_survival.csv'
    survival = pd.read_csv(surv_file)
    
    print(f"  Scores: {scores_df.shape}")
    print(f"  Clinical: {clinical.shape}")
    print(f"  Survival: {survival.shape}")
    
    return scores_df, clinical, survival

def merge_data(scores_df, clinical, survival):
    """Merge all data sources"""
    print("\nMerging data...")
    
    # Merge scores with clinical
    merged = scores_df.reset_index().merge(
        clinical, left_on='index', right_on='sample_id', how='inner'
    )
    
    # Merge with survival
    if 'sample' in survival.columns:
        merged = merged.merge(
            survival, left_on='sample_id', right_on='sample', how='left'
        )
    
    print(f"  Merged dataset: {merged.shape}")
    
    return merged

def survival_analysis(merged_df, score_col='FerroRadio_Score'):
    """
    Perform survival analysis
    - Kaplan-Meier curves
    - Log-rank test
    - Cox proportional hazards model
    """
    print("\n" + "=" * 60)
    print("Survival Analysis")
    print("=" * 60)
    
    # Filter patients with survival data
    surv_data = merged_df.dropna(subset=['OS.time', 'OS']).copy()
    
    if len(surv_data) < 100:
        print(f"  ✗ Insufficient survival data: {len(surv_data)} samples")
        return None
    
    print(f"  Samples with survival data: {len(surv_data)}")
    
    # Convert survival time from days to years
    surv_data['OS_years'] = surv_data['OS.time'] / 365.25
    
    # Create high/low score groups
    median_score = surv_data[score_col].median()
    surv_data['score_group'] = (surv_data[score_col] >= median_score).astype(int)
    surv_data['score_label'] = surv_data['score_group'].map({0: 'Low', 1: 'High'})
    
    print(f"  Median {score_col}: {median_score:.3f}")
    print(f"  High score group: {(surv_data['score_group']==1).sum()}")
    print(f"  Low score group: {(surv_data['score_group']==0).sum()}")
    
    # Kaplan-Meier analysis
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    high_mask = surv_data['score_group'] == 1
    low_mask = surv_data['score_group'] == 0
    
    kmf_high.fit(
        surv_data.loc[high_mask, 'OS_years'],
        surv_data.loc[high_mask, 'OS'],
        label='High Score'
    )
    kmf_low.fit(
        surv_data.loc[low_mask, 'OS_years'],
        surv_data.loc[low_mask, 'OS'],
        label='Low Score'
    )
    
    # Log-rank test
    results = logrank_test(
        surv_data.loc[high_mask, 'OS_years'],
        surv_data.loc[low_mask, 'OS_years'],
        surv_data.loc[high_mask, 'OS'],
        surv_data.loc[low_mask, 'OS']
    )
    
    print(f"\n  Log-rank test p-value: {results.p_value:.4e}")
    
    # Plot KM curves
    fig, ax = plt.subplots(figsize=(10, 7))
    kmf_high.plot_survival_function(ax=ax, ci_show=True, color='red')
    kmf_low.plot_survival_function(ax=ax, ci_show=True, color='blue')
    
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Overall Survival Probability')
    ax.set_title(f'Kaplan-Meier Curves by {score_col}\n(p={results.p_value:.4e})')
    ax.legend(loc='lower left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/figures/km_survival_{score_col}.png', dpi=300)
    print(f"  ✓ Saved: km_survival_{score_col}.png")
    plt.close()
    
    # Cox regression
    cox_data = surv_data[['OS_years', 'OS', score_col, 'age_at_initial_pathologic_diagnosis', 'gender']].dropna()
    
    if len(cox_data) > 50:
        cph = CoxPHFitter()
        cox_data['gender'] = (cox_data['gender'] == 'female').astype(int)
        cox_data = cox_data.rename(columns={'age_at_initial_pathologic_diagnosis': 'age'})
        
        try:
            cph.fit(cox_data, duration_col='OS_years', event_col='OS')
            
            print(f"\n  Cox Regression Results:")
            print(f"    {score_col} HR: {np.exp(cph.params_[score_col]):.3f} "
                  f"(95% CI: {np.exp(cph.confidence_intervals_.loc[score_col, 'lower-bound']):.3f}-"
                  f"{np.exp(cph.confidence_intervals_.loc[score_col, 'upper-bound']):.3f})")
            print(f"    p-value: {cph.summary['p'][score_col]:.4e}")
            
            # Save Cox results
            cph.summary.to_csv(f'{RESULTS_DIR}/tables/cox_regression_{score_col}.csv')
            print(f"  ✓ Saved: cox_regression_{score_col}.csv")
            
        except Exception as e:
            print(f"  ✗ Cox regression error: {e}")
    
    return results.p_value

def immunotherapy_validation(merged_df):
    """
    Validate specifically in immunotherapy patients
    """
    print("\n" + "=" * 60)
    print("Immunotherapy-Specific Validation")
    print("=" * 60)
    
    # Filter immunotherapy patients
    if 'received_immunotherapy' not in merged_df.columns:
        print("  ✗ No immunotherapy information")
        return None
    
    immuno_patients = merged_df[merged_df['received_immunotherapy'] == True].copy()
    print(f"  Patients with immunotherapy: {len(immuno_patients)}")
    
    if len(immuno_patients) < 50:
        print("  ✗ Insufficient immunotherapy patients")
        return None
    
    # Check survival data
    immuno_surv = immuno_patients.dropna(subset=['OS.time', 'OS'])
    print(f"  Immuno patients with survival: {len(immuno_surv)}")
    
    if len(immuno_surv) < 30:
        print("  ✗ Insufficient survival data for immuno patients")
        return None
    
    # Survival analysis in immuno patients
    immuno_surv['OS_years'] = immuno_surv['OS.time'] / 365.25
    
    # Stratify by FerroImmuno score
    median_score = immuno_surv['FerroImmuno_Score'].median()
    immuno_surv['score_group'] = (immuno_surv['FerroImmuno_Score'] >= median_score).astype(int)
    
    # KM analysis
    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()
    
    high_mask = immuno_surv['score_group'] == 1
    low_mask = immuno_surv['score_group'] == 0
    
    kmf_high.fit(immuno_surv.loc[high_mask, 'OS_years'], immuno_surv.loc[high_mask, 'OS'], label='High FerroImmuno')
    kmf_low.fit(immuno_surv.loc[low_mask, 'OS_years'], immuno_surv.loc[low_mask, 'OS'], label='Low FerroImmuno')
    
    results = logrank_test(
        immuno_surv.loc[high_mask, 'OS_years'],
        immuno_surv.loc[low_mask, 'OS_years'],
        immuno_surv.loc[high_mask, 'OS'],
        immuno_surv.loc[low_mask, 'OS']
    )
    
    print(f"\n  Immuno patients - Log-rank p-value: {results.p_value:.4e}")
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 7))
    kmf_high.plot_survival_function(ax=ax, ci_show=True, color='darkred')
    kmf_low.plot_survival_function(ax=ax, ci_show=True, color='darkblue')
    
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Overall Survival Probability')
    ax.set_title(f'Immunotherapy Patients: Survival by FerroImmuno Score\n(p={results.p_value:.4e})')
    ax.legend(loc='lower left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/figures/km_immunotherapy_patients.png', dpi=300)
    print(f"  ✓ Saved: km_immunotherapy_patients.png")
    plt.close()
    
    # Compare immuno vs non-immuno
    non_immuno = merged_df[merged_df['received_immunotherapy'] == False]
    print(f"\n  Non-immuno patients: {len(non_immuno)}")
    
    if len(non_immuno) > 50:
        print(f"  Mean FerroImmuno Score:")
        print(f"    Immuno patients: {immuno_patients['FerroImmuno_Score'].mean():.3f}")
        print(f"    Non-immuno patients: {non_immuno['FerroImmuno_Score'].mean():.3f}")
        
        # T-test
        t_stat, p_val = stats.ttest_ind(
            immuno_patients['FerroImmuno_Score'].dropna(),
            non_immuno['FerroImmuno_Score'].dropna()
        )
        print(f"    T-test p-value: {p_val:.4f}")
    
    return results.p_value

def cancer_type_analysis(merged_df):
    """
    Analyze by cancer type
    """
    print("\n" + "=" * 60)
    print("Cancer Type-Specific Analysis")
    print("=" * 60)
    
    if 'cancer_type' not in merged_df.columns:
        print("  ✗ No cancer type information")
        return None
    
    # Count by cancer type
    cancer_counts = merged_df['cancer_type'].value_counts()
    print(f"  Number of cancer types: {len(cancer_counts)}")
    print(f"  Top 5 cancer types:")
    for cancer, count in cancer_counts.head().items():
        print(f"    {cancer}: {count}")
    
    # Analyze FerroRadio score distribution by cancer
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Box plot
    top_cancers = cancer_counts.head(15).index
    plot_data = merged_df[merged_df['cancer_type'].isin(top_cancers)]
    
    sns.boxplot(data=plot_data, x='cancer_type', y='FerroImmuno_Score', ax=axes[0])
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha='right')
    axes[0].set_title('FerroImmuno Score Distribution by Cancer Type (Top 15)')
    axes[0].set_ylabel('FerroImmuno Score')
    
    # Mean score by cancer
    mean_scores = merged_df.groupby('cancer_type')['FerroImmuno_Score'].mean().sort_values(ascending=False)
    mean_scores.head(20).plot(kind='bar', ax=axes[1], color='steelblue')
    axes[1].set_title('Mean FerroImmuno Score by Cancer Type (Top 20)')
    axes[1].set_ylabel('Mean FerroImmuno Score')
    axes[1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/figures/cancer_type_analysis.png', dpi=300)
    print(f"\n  ✓ Saved: cancer_type_analysis.png")
    plt.close()
    
    # Save summary
    summary = merged_df.groupby('cancer_type').agg({
        'FerroImmuno_Score': ['mean', 'std', 'count'],
        'FerroScore': ['mean', 'std'],
        'Immune_Score': ['mean', 'std']
    }).round(3)
    
    summary.to_csv(f'{RESULTS_DIR}/tables/cancer_type_summary.csv')
    print(f"  ✓ Saved: cancer_type_summary.csv")
    
    return mean_scores

def prognostic_value_by_stage(merged_df):
    """
    Analyze prognostic value by tumor stage
    """
    print("\n" + "=" * 60)
    print("Prognostic Value by Tumor Stage")
    print("=" * 60)
    
    if 'stage' not in merged_df.columns:
        print("  ✗ No stage information")
        return None
    
    # Filter patients with stage and survival data
    stage_data = merged_df.dropna(subset=['stage', 'OS.time', 'OS']).copy()
    
    if len(stage_data) < 100:
        print(f"  ✗ Insufficient data: {len(stage_data)}")
        return None
    
    print(f"  Samples with stage data: {len(stage_data)}")
    
    # Simplify stage
    def simplify_stage(stage):
        if pd.isna(stage):
            return 'Unknown'
        stage = str(stage).upper()
        if 'I' in stage and 'IV' not in stage and 'III' not in stage and 'II' not in stage:
            return 'Stage I'
        elif 'II' in stage and 'III' not in stage and 'IV' not in stage:
            return 'Stage II'
        elif 'III' in stage and 'IV' not in stage:
            return 'Stage III'
        elif 'IV' in stage:
            return 'Stage IV'
        else:
            return 'Other'
    
    stage_data['stage_simple'] = stage_data['stage'].apply(simplify_stage)
    
    print(f"\n  Stage distribution:")
    print(stage_data['stage_simple'].value_counts())
    
    # Analyze prognostic value within each stage
    results = []
    for stage in ['Stage I', 'Stage II', 'Stage III', 'Stage IV']:
        stage_subset = stage_data[stage_data['stage_simple'] == stage]
        
        if len(stage_subset) < 30:
            continue
        
        # Stratify by FerroImmuno score
        median_score = stage_subset['FerroImmuno_Score'].median()
        stage_subset['score_group'] = (stage_subset['FerroImmuno_Score'] >= median_score).astype(int)
        
        # Log-rank test
        high_mask = stage_subset['score_group'] == 1
        low_mask = stage_subset['score_group'] == 0
        
        if high_mask.sum() < 10 or low_mask.sum() < 10:
            continue
        
        lr_results = logrank_test(
            stage_subset.loc[high_mask, 'OS.time'] / 365.25,
            stage_subset.loc[low_mask, 'OS.time'] / 365.25,
            stage_subset.loc[high_mask, 'OS'],
            stage_subset.loc[low_mask, 'OS']
        )
        
        results.append({
            'Stage': stage,
            'N': len(stage_subset),
            'p_value': lr_results.p_value,
            'significant': lr_results.p_value < 0.05
        })
        
        print(f"\n  {stage}: N={len(stage_subset)}, p={lr_results.p_value:.4f}")
    
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(f'{RESULTS_DIR}/tables/stage_stratified_analysis.csv', index=False)
        print(f"\n  ✓ Saved: stage_stratified_analysis.csv")
        return results_df
    
    return None

def main():
    """Main validation function"""
    print("=" * 60)
    print("FerroScore-Immuno: Validation Analysis")
    print("=" * 60)
    
    # 1. Load data
    scores_df, clinical, survival = load_data()
    
    # 2. Merge data
    merged_df = merge_data(scores_df, clinical, survival)
    
    # 3. Overall survival analysis
    p_val_os = survival_analysis(merged_df, 'FerroImmuno_Score')
    p_val_ferro = survival_analysis(merged_df, 'FerroScore')
    p_val_immune = survival_analysis(merged_df, 'Immune_Score')
    
    # 4. Immunotherapy-specific validation
    immuno_pval = immunotherapy_validation(merged_df)
    
    # 5. Cancer type analysis
    cancer_scores = cancer_type_analysis(merged_df)
    
    # 6. Stage-stratified analysis
    stage_results = prognostic_value_by_stage(merged_df)
    
    # 7. Summary
    print("\n" + "=" * 60)
    print("Validation Summary")
    print("=" * 60)
    
    summary = {
        'Analysis': [
            'Overall Survival (FerroImmuno)',
            'Overall Survival (FerroScore)',
            'Overall Survival (Immune)',
            'Immunotherapy Patients',
            'Cancer Types Analyzed',
            'Stage-stratified'
        ],
        'P-value': [
            p_val_os if p_val_os else 'N/A',
            p_val_ferro if p_val_ferro else 'N/A',
            p_val_immune if p_val_immune else 'N/A',
            immuno_pval if immuno_pval else 'N/A',
            len(cancer_scores) if cancer_scores is not None else 'N/A',
            'Completed' if stage_results is not None else 'N/A'
        ]
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(f'{RESULTS_DIR}/tables/validation_summary.csv', index=False)
    print(f"\n✓ Saved: validation_summary.csv")
    
    print("\n" + "=" * 60)
    print("Validation Complete!")
    print("=" * 60)
    print("\nNext step: Run 06_visualization.py")

if __name__ == "__main__":
    main()
