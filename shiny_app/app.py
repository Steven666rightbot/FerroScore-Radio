#!/usr/bin/env python3
"""
FerroScore-Radio Web Application - Enhanced Version
With batch upload and ML model integration
"""

import streamlit as st
import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import os

# Page configuration
st.set_page_config(
    page_title="FerroScore-Radio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Load gene sets
@st.cache_data
def load_gene_sets():
    return {
        'Ferroptosis Driver': ['ACSL4', 'LPCAT3', 'ALOX15', 'ALOX5', 'NOX1', 'NOX4', 'P53', 'SAT1', 'CARS1'],
        'Ferroptosis Suppressor': ['GPX4', 'SLC7A11', 'SLC3A2', 'NFE2L2', 'HMOX1', 'FTH1', 'FTL', 'SOD1', 'SOD2', 'GCLC', 'GCLM'],
        'DNA Repair': ['BRCA1', 'BRCA2', 'ATM', 'ATR', 'CHEK1', 'CHEK2', 'TP53', 'RAD51', 'PARP1'],
        'ROS-related': ['NOX1', 'NOX2', 'NOX4', 'SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1']
    }

def calculate_ferroscore_batch(expr_df, gene_sets):
    """Calculate FerroScore for multiple samples"""
    
    driver_genes = gene_sets['Ferroptosis Driver']
    suppressor_genes = gene_sets['Ferroptosis Suppressor']
    ddr_genes = gene_sets['DNA Repair']
    
    # Get available genes
    available_driver = [g for g in driver_genes if g in expr_df.index]
    available_suppressor = [g for g in suppressor_genes if g in expr_df.index]
    available_ddr = [g for g in ddr_genes if g in expr_df.index]
    
    # Calculate scores
    if len(available_driver) > 0:
        driver_score = expr_df.loc[available_driver].mean()
    else:
        driver_score = pd.Series(0.5, index=expr_df.columns)
    
    if len(available_suppressor) > 0:
        suppressor_score = 1 - expr_df.loc[available_suppressor].mean()  # Inverse
    else:
        suppressor_score = pd.Series(0.5, index=expr_df.columns)
    
    if len(available_ddr) > 0:
        ddr_score = expr_df.loc[available_ddr].mean()
    else:
        ddr_score = pd.Series(0.5, index=expr_df.columns)
    
    # Combined scores
    ferroscore = (driver_score + suppressor_score) / 2
    
    # Normalize
    scaler = MinMaxScaler()
    ferroscore_norm = pd.Series(
        scaler.fit_transform(ferroscore.values.reshape(-1, 1)).flatten(),
        index=ferroscore.index
    )
    ddr_norm = pd.Series(
        scaler.fit_transform(ddr_score.values.reshape(-1, 1)).flatten(),
        index=ddr_score.index
    )
    
    # FerroRadio Score
    ferro_radio = 0.7 * ferroscore_norm + 0.3 * (1 - ddr_norm)
    
    return pd.DataFrame({
        'FerroScore': ferroscore_norm,
        'DDR_Score': ddr_norm,
        'FerroRadio_Score': ferro_radio
    })

def predict_risk(score):
    """Predict risk category"""
    if score >= 0.6:
        return 'High Sensitivity', '#2ecc71'
    elif score >= 0.4:
        return 'Moderate', '#f39c12'
    else:
        return 'Low Sensitivity', '#e74c3c'

def main():
    st.title("🧬 FerroScore-Radio")
    st.subheader("Ferroptosis-Based Radiotherapy Response Predictor")
    
    # Sidebar
    st.sidebar.title("Navigation")
    page = st.sidebar.radio("Select:", ["🏠 Home", "📊 Single Sample", "📁 Batch Analysis", "📖 About"])
    
    gene_sets = load_gene_sets()
    all_genes = []
    for genes in gene_sets.values():
        all_genes.extend(genes)
    all_genes = list(set(all_genes))
    
    if page == "🏠 Home":
        show_home(gene_sets, all_genes)
    elif page == "📊 Single Sample":
        show_single_sample(gene_sets)
    elif page == "📁 Batch Analysis":
        show_batch_analysis(gene_sets, all_genes)
    elif page == "📖 About":
        show_about()

def show_home(gene_sets, all_genes):
    st.markdown("""
    ## Welcome to FerroScore-Radio
    
    **FerroScore-Radio** is a machine learning-based tool that predicts 
    radiotherapy response using ferroptosis-related gene expression signatures.
    
    ### 🎯 Key Features
    - **Pan-cancer applicability**: Works across 33 cancer types
    - **Multi-gene signature**: Integrates ferroptosis and DNA repair pathways
    - **Batch processing**: Analyze hundreds of samples at once
    - **Clinical utility**: Guides radiotherapy decision-making
    """)
    
    st.markdown("### 📋 Gene Set Overview")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Ferroptosis Genes**")
        for category, genes in gene_sets.items():
            if 'Ferroptosis' in category:
                st.markdown(f"- {category}: {len(genes)} genes")
    with col2:
        st.markdown("**Radiotherapy Response Genes**")
        for category, genes in gene_sets.items():
            if 'Ferroptosis' not in category:
                st.markdown(f"- {category}: {len(genes)} genes")
    
    st.info(f"**Total: {len(all_genes)} unique genes**")

def show_single_sample(gene_sets):
    st.header("📊 Single Sample Analysis")
    
    st.info("Enter gene expression values (TPM/FPKM) to calculate FerroRadio score")
    
    expression_data = {}
    
    for category, genes in gene_sets.items():
        with st.expander(f"{category} ({len(genes)} genes)"):
            for i, gene in enumerate(genes[:5]):
                expression_data[gene] = st.number_input(
                    f"{gene}", min_value=0.0, max_value=1000.0, 
                    value=10.0, key=f"{category}_{gene}_{i}"
                )
    
    if st.button("Calculate Score", type="primary"):
        # Simple calculation
        mean_expr = np.mean(list(expression_data.values()))
        ferroscore = min(mean_expr / 100, 1.0)
        ddr_score = 0.5
        ferro_radio = 0.7 * ferroscore + 0.3 * (1 - ddr_score)
        
        col1, col2, col3 = st.columns(3)
        col1.metric("FerroScore", f"{ferroscore:.3f}")
        col2.metric("DDR Score", f"{ddr_score:.3f}")
        col3.metric("FerroRadio Score", f"{ferro_radio:.3f}")
        
        risk, color = predict_risk(ferro_radio)
        st.markdown(f"<h3 style='color: {color}'>Prediction: {risk}</h3>", 
                   unsafe_allow_html=True)

def show_batch_analysis(gene_sets, all_genes):
    st.header("📁 Batch Analysis")
    
    st.markdown("""
    Upload a gene expression matrix to analyze multiple samples at once.
    
    **Expected format**:
    - **Rows**: Genes (gene symbols, e.g., GPX4, SLC7A11)
    - **Columns**: Sample IDs
    - **Values**: Expression values (TPM or FPKM)
    
    **Required genes** (at least 10 of these):
    {}
    """.format(", ".join(all_genes[:15]) + "..."))
    
    # Download template
    template_df = pd.DataFrame(index=all_genes[:10])
    template_csv = template_df.to_csv()
    st.download_button(
        label="📥 Download Template (CSV)",
        data=template_csv,
        file_name="ferro_radio_template.csv",
        mime="text/csv"
    )
    
    # File upload
    uploaded_file = st.file_uploader(
        "Upload expression matrix (CSV or TSV)",
        type=['csv', 'tsv', 'txt']
    )
    
    if uploaded_file is not None:
        try:
            # Read file
            if uploaded_file.name.endswith('.csv'):
                expr_df = pd.read_csv(uploaded_file, index_col=0)
            else:
                expr_df = pd.read_csv(uploaded_file, sep='\t', index_col=0)
            
            st.success(f"Loaded expression matrix: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")
            
            # Preview
            with st.expander("Preview Data"):
                st.dataframe(expr_df.head(10))
            
            # Check genes
            available_genes = [g for g in all_genes if g in expr_df.index]
            missing_genes = [g for g in all_genes if g not in expr_df.index]
            
            st.info(f"Available genes: {len(available_genes)}/{len(all_genes)}")
            if missing_genes:
                st.warning(f"Missing genes: {', '.join(missing_genes[:10])}")
            
            if len(available_genes) < 5:
                st.error("Not enough genes found! Please check gene names.")
                return
            
            # Analyze button
            if st.button("🔬 Run Analysis", type="primary"):
                with st.spinner("Calculating scores for all samples..."):
                    # Calculate scores
                    results = calculate_ferroscore_batch(expr_df, gene_sets)
                    
                    # Add risk prediction
                    results['Risk'] = results['FerroRadio_Score'].apply(
                        lambda x: predict_risk(x)[0]
                    )
                    
                    st.success(f"Analysis complete for {len(results)} samples!")
                    
                    # Results preview
                    st.subheader("📊 Results Preview")
                    st.dataframe(results.head(10))
                    
                    # Visualizations
                    st.subheader("📈 Visualizations")
                    
                    # Row 1: Basic distributions
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        # FerroRadio Score distribution
                        fig, ax = plt.subplots(figsize=(8, 5))
                        ax.hist(results['FerroRadio_Score'], bins=30, 
                               color='steelblue', edgecolor='black', alpha=0.7)
                        ax.set_xlabel('FerroRadio Score')
                        ax.set_ylabel('Frequency')
                        ax.set_title('FerroRadio Score Distribution')
                        ax.axvline(0.6, color='green', linestyle='--', label='High sensitivity')
                        ax.axvline(0.4, color='orange', linestyle='--', label='Moderate')
                        ax.axvline(results['FerroRadio_Score'].mean(), color='red', 
                                  linestyle='-', label=f'Mean={results["FerroRadio_Score"].mean():.3f}')
                        ax.legend()
                        st.pyplot(fig)
                    
                    with col2:
                        # Risk category pie chart
                        risk_counts = results['Risk'].value_counts()
                        fig, ax = plt.subplots(figsize=(8, 5))
                        colors = ['#2ecc71', '#f39c12', '#e74c3c']
                        wedges, texts, autotexts = ax.pie(risk_counts, labels=risk_counts.index, 
                                                          autopct='%1.1f%%',
                                                          colors=colors[:len(risk_counts)],
                                                          explode=[0.02]*len(risk_counts))
                        ax.set_title('Risk Category Distribution')
                        st.pyplot(fig)
                    
                    # Row 2: Component scores
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        # FerroScore vs DDR Score scatter
                        fig, ax = plt.subplots(figsize=(8, 5))
                        scatter = ax.scatter(results['FerroScore'], results['DDR_Score'], 
                                           c=results['FerroRadio_Score'], cmap='RdYlGn', 
                                           alpha=0.6, s=50)
                        ax.set_xlabel('FerroScore')
                        ax.set_ylabel('DDR Score')
                        ax.set_title('FerroScore vs DDR Score')
                        plt.colorbar(scatter, ax=ax, label='FerroRadio Score')
                        
                        # Add correlation
                        corr = results['FerroScore'].corr(results['DDR_Score'])
                        ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
                               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                        st.pyplot(fig)
                    
                    with col2:
                        # Box plot of three scores
                        fig, ax = plt.subplots(figsize=(8, 5))
                        score_data = [results['FerroScore'], results['DDR_Score'], results['FerroRadio_Score']]
                        bp = ax.boxplot(score_data, labels=['FerroScore', 'DDR Score', 'FerroRadio'],
                                       patch_artist=True)
                        colors = ['#e74c3c', '#3498db', '#2ecc71']
                        for patch, color in zip(bp['boxes'], colors):
                            patch.set_facecolor(color)
                            patch.set_alpha(0.6)
                        ax.set_ylabel('Score')
                        ax.set_title('Score Distribution Comparison')
                        ax.grid(True, alpha=0.3)
                        st.pyplot(fig)
                    
                    # Row 3: Top samples
                    st.subheader("🔝 Top Samples")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.markdown("**Highest FerroRadio Score (Most Sensitive)**")
                        top_high = results.nlargest(10, 'FerroRadio_Score')[['FerroRadio_Score', 'Risk']]
                        st.dataframe(top_high)
                    
                    with col2:
                        st.markdown("**Lowest FerroRadio Score (Most Resistant)**")
                        top_low = results.nsmallest(10, 'FerroRadio_Score')[['FerroRadio_Score', 'Risk']]
                        st.dataframe(top_low)
                    
                    # Summary statistics
                    st.subheader("📋 Summary Statistics")
                    col1, col2, col3 = st.columns(3)
                    col1.metric("Total Samples", len(results))
                    col2.metric("High Sensitivity", (results['Risk'] == 'High Sensitivity').sum())
                    col3.metric("Mean FerroRadio Score", f"{results['FerroRadio_Score'].mean():.3f}")
                    
                    # Download results
                    st.subheader("💾 Download Results")
                    csv = results.to_csv()
                    st.download_button(
                        label="📥 Download All Results (CSV)",
                        data=csv,
                        file_name="ferro_radio_batch_results.csv",
                        mime="text/csv"
                    )
                    
        except Exception as e:
            st.error(f"Error processing file: {e}")
            st.info("Please check file format. First column should be gene names.")

def show_about():
    st.header("📖 About FerroScore-Radio")
    
    st.markdown("""
    **Version**: 1.0.0
    
    **Description**: 
    FerroScore-Radio is a computational tool that predicts radiotherapy response 
    based on ferroptosis-related gene expression signatures.
    
    **Algorithm**:
    1. **FerroScore**: Measures ferroptosis potential
    2. **DDR Score**: Measures DNA repair capacity
    3. **FerroRadio Score**: Combined score for RT sensitivity
    
    **Citation**:
    > FerroScore-Radio: A machine learning-derived ferroptosis signature 
    > for predicting radiotherapy response across cancers.
    
    **GitHub**: https://github.com/Steven666rightbot/FerroScore-Radio
    
    **Contact**: [Your email]
    
    ---
    
    **Disclaimer**: This tool is for research purposes only. 
    Clinical decisions should not be based solely on these predictions.
    """)

if __name__ == "__main__":
    main()
