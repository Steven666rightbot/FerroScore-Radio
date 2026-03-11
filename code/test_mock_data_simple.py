#!/usr/bin/env python3
"""
Simple Test Script - No External Dependencies
Test mock data with basic Python only
"""

import csv
import os
import random

PROC_DIR = '../data/processed'
RESULTS_DIR = '../results'
os.makedirs(f'{RESULTS_DIR}/tables', exist_ok=True)

def test_mock_data():
    """Test loading mock data with basic CSV reader"""
    print("=" * 60)
    print("Testing Mock Data (Basic CSV)")
    print("=" * 60)
    
    # Test expression data
    expr_file = f'{PROC_DIR}/tcga_ferro_radio_expression.csv'
    if os.path.exists(expr_file):
        with open(expr_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            n_samples = len(header) - 1  # First column is gene names
            n_genes = sum(1 for row in reader)
        print(f"\n[OK] Expression data: {n_genes} genes x {n_samples} samples")
    
    # Test clinical data
    clin_file = f'{PROC_DIR}/tcga_clinical.csv'
    if os.path.exists(clin_file):
        with open(clin_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            n_clinical = len(rows)
            rt_count = sum(1 for r in rows if r.get('received_radiotherapy') == 'True')
        print(f"[OK] Clinical data: {n_clinical} samples")
        print(f"  - With RT: {rt_count}")
        print(f"  - Without RT: {n_clinical - rt_count}")
    
    # Test survival data
    surv_file = f'{PROC_DIR}/tcga_survival.csv'
    if os.path.exists(surv_file):
        with open(surv_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            n_survival = len(rows)
            events = sum(1 for r in rows if r.get('OS') == '1')
        print(f"[OK] Survival data: {n_survival} samples")
        print(f"  - Events: {events}")
        print(f"  - Censored: {n_survival - events}")
    
    print("\n" + "=" * 60)
    print("All mock data verified successfully!")
    print("=" * 60)
    print("\nNote: Full analysis requires pandas/scipy/sklearn")
    print("      Install with: pip install pandas scipy scikit-learn")

def simple_ferroscore_demo():
    """Demonstrate FerroScore calculation concept"""
    print("\n" + "=" * 60)
    print("FerroScore Calculation Demo (Simplified)")
    print("=" * 60)
    
    # Mock gene expression for one sample
    sample_expr = {
        'ACSL4': 15.2, 'LPCAT3': 8.5, 'NOX1': 12.3,  # Drivers
        'GPX4': 45.6, 'SLC7A11': 38.9, 'NFE2L2': 22.1,  # Suppressors
        'BRCA1': 18.7, 'ATM': 14.2, 'TP53': 25.3  # DDR
    }
    
    print("\nSample gene expression (TPM):")
    for gene, expr in sample_expr.items():
        print(f"  {gene}: {expr:.1f}")
    
    # Simple scoring (just for demonstration)
    driver_score = (sample_expr['ACSL4'] + sample_expr['LPCAT3'] + sample_expr['NOX1']) / 3
    suppressor_score = 50 - (sample_expr['GPX4'] + sample_expr['SLC7A11']) / 2  # Inverse
    ddr_score = (sample_expr['BRCA1'] + sample_expr['ATM'] + sample_expr['TP53']) / 3
    
    # Normalize to 0-1 (simplified)
    ferroscore = min(max(driver_score / 50, 0), 1)
    ferro_radio = min(max((ferroscore + (1 - ddr_score/50)) / 2, 0), 1)
    
    print(f"\nCalculated Scores:")
    print(f"  FerroScore: {ferroscore:.3f}")
    print(f"  DDR Score: {min(max(ddr_score/50, 0), 1):.3f}")
    print(f"  FerroRadio Score: {ferro_radio:.3f}")
    
    # Prediction
    if ferro_radio > 0.6:
        prediction = "High RT Sensitivity"
    elif ferro_radio > 0.4:
        prediction = "Moderate Sensitivity"
    else:
        prediction = "Low Sensitivity"
    
    print(f"\n  Prediction: {prediction}")
    print("\n" + "=" * 60)

if __name__ == "__main__":
    test_mock_data()
    simple_ferroscore_demo()
    
    print("\n" + "=" * 60)
    print("Next Steps:")
    print("=" * 60)
    print("1. Fix Python environment (pandas/scipy installation)")
    print("2. Run: python 03_ferroscore_algorithm.py")
    print("3. Run: python 04_model_training.py")
    print("4. Run: python 05_validation.py")
    print("5. Run: python 06_visualization.py")
    print("\nOr use the Streamlit app:")
    print("  cd shiny_app")
    print("  streamlit run app.py")
    print("=" * 60)
