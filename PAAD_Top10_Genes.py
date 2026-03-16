# ============================================================================
# PAAD GENOME-WIDE ANALYSIS - TARGET GENE VISUALIZATION
# Author: Lubanah Younes
# Date: 2026
# Description: Generate survival and SHAP plots for specific target genes
# ============================================================================

# ------------------------------------------------------------
# Import required libraries
# Author: Lubanah Younes
# ------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import os
import shap
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

# ============================================================================
# SECTION 1: DATA LOADING
# This section loads the gene expression and survival data files
# Author: Lubanah Younes
# ============================================================================

print("="*60)
print("🔬 PAAD TARGET GENE ANALYSIS")
print(f"👩‍🔬 Principal Investigator: Lubanah Younes")
print("="*60)

# Define file paths
expr_file = "TCGA.PAAD.sampleMap_HiSeqV2.gz"
surv_file = "survival_PAAD_survival.txt"

# Check if files exist
if not os.path.exists(expr_file) or not os.path.exists(surv_file):
    print("❌ Error: Data files not found!")
    exit()

# Load gene expression data
print("\n📂 Loading expression data...")
expr = pd.read_csv(expr_file, compression='gzip', sep='\t', index_col=0)
print(f"✅ Loaded {expr.shape[0]:,} genes and {expr.shape[1]} samples")

# Load survival data
print("📂 Loading survival data...")
surv = pd.read_csv(surv_file, sep='\t')
print(f"✅ Loaded {surv.shape[0]} patients")

# ============================================================================
# SECTION 2: DATA PREPARATION
# Clean and align survival data with expression data
# Author: Lubanah Younes
# ============================================================================

# Clean survival data
surv_clean = surv.dropna(subset=['OS.time', 'OS']).copy()
print(f"\n✅ Clean survival data: {surv_clean.shape[0]} patients with complete data")

# Align samples
X = expr.T
y = surv_clean.set_index('sample').loc[X.index, 'OS']
valid = y.notna()
X = X[valid]
y = y[valid]
print(f"✅ Aligned data: {X.shape[0]} samples, {X.shape[1]:,} genes")

# ============================================================================
# SECTION 3: TARGET GENE SELECTION
# Define the genes we want to analyze (top novel candidates)
# Author: Lubanah Younes
# ============================================================================

# List of top novel genes from our analysis
target_genes = [
    "PHLDB3",    # p = 0.000001
    "SLURP1",    # p = 0.000002
    "MYEOV",     # p = 0.000002
    "USP20",     # p = 0.000003
    "LOC651250", # p = 0.000005
    "DEF8",      # p = 0.000005
    "EPS8",      # p = 0.000007
    "NCAM1",     # p = 0.000009
    "FAM123A",   # p = 0.000016
    "MYO5B"      # p = 0.000017
]

print(f"\n🎯 Target genes to analyze: {', '.join(target_genes)}")

# ============================================================================
# SECTION 4: SURVIVAL ANALYSIS FOR EACH GENE
# Generate Kaplan-Meier plots and calculate log-rank p-values
# Author: Lubanah Younes
# ============================================================================

# Store results
results = []

for gene in target_genes:
    print(f"\n{'='*50}")
    print(f"🔬 Analyzing gene: {gene}")
    print(f"{'='*50}")
    
    # Check if gene exists in data
    if gene not in expr.index:
        print(f"⚠️ Gene {gene} not found in expression data!")
        continue
    
    # Get gene expression data
    gene_expr = expr.loc[gene]
    median_expr = gene_expr.median()
    
    # Split patients into high and low expression groups
    high_patients = gene_expr > median_expr
    low_patients = ~high_patients
    
    # Get survival data for each group
    clin_high = surv_clean[surv_clean['sample'].isin(gene_expr[high_patients].index)]
    clin_low = surv_clean[surv_clean['sample'].isin(gene_expr[low_patients].index)]
    
    print(f"   High expression group: {len(clin_high)} patients")
    print(f"   Low expression group: {len(clin_low)} patients")
    
    # Perform log-rank test
    results_logrank = logrank_test(
        clin_high['OS.time'], clin_low['OS.time'],
        event_observed_A=clin_high['OS'],
        event_observed_B=clin_low['OS']
    )
    p_value = results_logrank.p_value
    
    # Store result
    results.append({
        'gene': gene,
        'p_value': p_value,
        'high_n': len(clin_high),
        'low_n': len(clin_low)
    })
    
    print(f"   📊 Log-rank p-value: {p_value:.6f}")
    
    # ============================================================================
    # SECTION 4.1: GENERATE KAPLAN-MEIER PLOT
    # Create survival curve visualization
    # Author: Lubanah Younes
    # ============================================================================
    
    plt.figure(figsize=(10, 7))
    
    # Fit Kaplan-Meier estimators
    kmf_high = KaplanMeierFitter().fit(
        clin_high['OS.time'], 
        clin_high['OS'], 
        label=f'High Expression (n={len(clin_high)})'
    )
    
    kmf_low = KaplanMeierFitter().fit(
        clin_low['OS.time'], 
        clin_low['OS'], 
        label=f'Low Expression (n={len(clin_low)})'
    )
    
    # Plot survival curves
    ax = kmf_high.plot_survival_function(color='red', linewidth=2)
    kmf_low.plot_survival_function(ax=ax, color='blue', linewidth=2)
    
    # Add p-value to plot
    if p_value < 0.001:
        p_text = f'p = {p_value:.6f} ***'
    elif p_value < 0.01:
        p_text = f'p = {p_value:.4f} **'
    elif p_value < 0.05:
        p_text = f'p = {p_value:.4f} *'
    else:
        p_text = f'p = {p_value:.4f}'
    
    plt.text(0.6, 0.15, p_text, transform=ax.transAxes,
             fontsize=14, fontweight='bold',
             bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Add significance stars explanation
    if p_value < 0.05:
        plt.text(0.6, 0.08, '* p<0.05, ** p<0.01, *** p<0.001',
                transform=ax.transAxes, fontsize=10, style='italic')
    
    # Customize plot
    plt.title(f'PAAD Survival Analysis: {gene}', fontsize=16, fontweight='bold')
    plt.xlabel('Time (Days)', fontsize=14)
    plt.ylabel('Survival Probability', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='lower left', fontsize=12)
    
    # Add author attribution
    plt.text(0.02, 0.02, f'Analysis by Lubanah Younes, 2026',
             transform=ax.transAxes, fontsize=8, style='italic', alpha=0.7)
    
    # Save plot
    plot_filename = f'PAAD_survival_{gene}.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   ✅ Survival plot saved: {plot_filename}")

# ============================================================================
# SECTION 5: RESULTS SUMMARY
# Compile and display all results
# Author: Lubanah Younes
# ============================================================================

print("\n" + "="*60)
print("📊 FINAL RESULTS SUMMARY")
print("="*60)

# Create results dataframe
results_df = pd.DataFrame(results)
results_df = results_df.sort_values('p_value')

print("\n🔬 Genes ranked by statistical significance:")
print("-"*50)
for idx, row in results_df.iterrows():
    stars = ""
    if row['p_value'] < 0.001:
        stars = "***"
    elif row['p_value'] < 0.01:
        stars = "**"
    elif row['p_value'] < 0.05:
        stars = "*"
    
    print(f"   {row['gene']:<12} p = {row['p_value']:.6f} {stars}")

# Save results to CSV
results_df.to_csv('PAAD_target_genes_results.csv', index=False)
print(f"\n💾 Results saved to: PAAD_target_genes_results.csv")

print("\n" + "="*60)
print("✅ Analysis complete!")
print(f"👩‍🔬 Principal Investigator: Lubanah Younes")
print("📊 Generated plots and results are ready for publication")
print("="*60)