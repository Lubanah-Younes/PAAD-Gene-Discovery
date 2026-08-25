# =============================================================================
# PAAD GENE COMPARISON – PHLDB3 vs KNOWN GENES
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script compares the prognostic performance of PHLDB3
# with known PAAD genes (TP53, KRAS, CDKN2A, SMAD4).
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# STEP 1 – LOAD DATA
# =============================================================================
print('[1] Loading data...')

# Load expression data
expr = pd.read_csv('TCGA.PAAD.sampleMap_HiSeqV2.gz',
                   compression='gzip', sep='\t', index_col=0)
expr = expr.T
expr.index = expr.index.str[:12]
expr = expr[~expr.index.duplicated(keep='first')]
print(f'    Expression: {expr.shape[0]} samples, {expr.shape[1]} genes')

# Load survival data
surv = pd.read_csv('survival_PAAD_survival.txt', sep='\t')
surv['sample_id'] = surv['sample'].str[:12]
surv = surv[~surv['sample_id'].duplicated(keep='first')]
print(f'    Survival: {surv.shape[0]} patients')

# =============================================================================
# STEP 2 – DEFINE GENES TO COMPARE
# =============================================================================
print('[2] Defining genes...')

genes = {
    'PHLDB3': 'Our Gene',
    'TP53': 'Known',
    'KRAS': 'Known',
    'CDKN2A': 'Known',
    'SMAD4': 'Known'
}

# Check which genes exist
available_genes = {}
for gene, label in genes.items():
    if gene in expr.columns:
        available_genes[gene] = label
    else:
        print(f'    ⚠️ {gene} not found')

print(f'    Available genes: {list(available_genes.keys())}')

# =============================================================================
# STEP 3 – CALCULATE HR AND P-VALUE
# =============================================================================
print('[3] Calculating survival metrics...')

results = []

for gene, label in available_genes.items():
    # Get expression
    gene_expr = expr[gene]
    
    # Merge with survival
    merged = pd.DataFrame({
        'expression': gene_expr,
        'status': surv.set_index('sample_id')['OS'],
        'time': surv.set_index('sample_id')['OS.time']
    }).dropna()
    
    # Split by median
    median = merged['expression'].median()
    merged['group'] = np.where(merged['expression'] > median, 'High', 'Low')
    
    # Log-rank test
    high = merged[merged['group'] == 'High']
    low = merged[merged['group'] == 'Low']
    p_val = logrank_test(high['time'], low['time'], 
                         event_observed_A=high['status'],
                         event_observed_B=low['status']).p_value
    
    # Cox regression
    cph = CoxPHFitter()
    cph.fit(merged, duration_col='time', event_col='status', formula='expression')
    
    # Get HR and CI from summary
    summary = cph.summary
    hr = summary.loc['expression', 'exp(coef)']
    ci_lower = summary.loc['expression', 'exp(coef) lower 95%']
    ci_upper = summary.loc['expression', 'exp(coef) upper 95%']
    
    results.append({
        'Gene': gene,
        'Type': label,
        'HR': hr,
        'CI_lower': ci_lower,
        'CI_upper': ci_upper,
        'P_value': p_val,
        'Significant': p_val < 0.05
    })
    
    print(f'    {gene}: HR = {hr:.3f}, p = {p_val:.6f}')

# =============================================================================
# STEP 4 – CREATE RESULTS TABLE
# =============================================================================
print('[4] Creating results table...')

results_df = pd.DataFrame(results)
results_df = results_df.sort_values('P_value')

print('\nComparison Results:')
print(results_df.to_string(index=False))

# Save
results_df.to_csv('PAAD_Gene_Comparison.csv', index=False)
print('\nSaved: PAAD_Gene_Comparison.csv')

# =============================================================================
# STEP 5 – BAR PLOT
# =============================================================================
print('[5] Generating bar plot...')

plt.figure(figsize=(10, 6))
colors = ['red' if row['Gene'] == 'PHLDB3' else 'blue' for _, row in results_df.iterrows()]
plt.barh(results_df['Gene'], -np.log10(results_df['P_value']), color=colors)
plt.xlabel('-Log10(P-value)', fontsize=12)
plt.ylabel('Gene', fontsize=12)
plt.title('PHLDB3 vs Known PAAD Genes: Prognostic Significance', fontsize=14, fontweight='bold')
plt.axvline(x=-np.log10(0.05), color='black', linestyle='--', label='p = 0.05')
plt.legend(['p = 0.05', 'PHLDB3', 'Known Genes'])
plt.tight_layout()
plt.savefig('PAAD_Gene_Comparison_Plot.png', dpi=300)
print('    Plot saved: PAAD_Gene_Comparison_Plot.png')

# =============================================================================
# STEP 6 – SURVIVAL CURVES
# =============================================================================
print('[6] Generating survival curves...')

top_genes = results_df.head(4)['Gene'].tolist()

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for idx, gene in enumerate(top_genes):
    gene_expr = expr[gene]
    merged = pd.DataFrame({
        'expression': gene_expr,
        'status': surv.set_index('sample_id')['OS'],
        'time': surv.set_index('sample_id')['OS.time']
    }).dropna()
    
    median = merged['expression'].median()
    merged['group'] = np.where(merged['expression'] > median, 'High', 'Low')
    
    kmf = KaplanMeierFitter()
    ax = axes[idx]
    
    for group in ['High', 'Low']:
        subset = merged[merged['group'] == group]
        kmf.fit(subset['time'], event_observed=subset['status'], label=group)
        kmf.plot_survival_function(ax=ax)
    
    ax.set_title(f'{gene}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Days')
    ax.set_ylabel('Survival Probability')
    ax.grid(alpha=0.3)
    
    row = results_df[results_df['Gene'] == gene].iloc[0]
    ax.text(0.6, 0.1, f'HR = {row["HR"]:.2f}\np = {row["P_value"]:.6f}', 
            transform=ax.transAxes, fontsize=10)

plt.tight_layout()
plt.savefig('PAAD_Gene_Comparison_KM.png', dpi=300)
print('    Survival curves saved: PAAD_Gene_Comparison_KM.png')

print('\n✅ Gene comparison complete!')