# =============================================================================
# PAAD IMMUNE INFILTRATION ANALYSIS – PHLDB3 (TIMER ALTERNATIVE)
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script estimates immune cell abundance using cBioPortal
# data as a proxy for TIMER2.0 analysis.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# =============================================================================
# STEP 1 – LOAD DATA
# =============================================================================
print('[1] Loading expression data...')
expr = pd.read_csv('TCGA.PAAD.sampleMap_HiSeqV2.gz',
                   compression='gzip', sep='\t', index_col=0)
expr = expr.T
print(f'    Expression data: {expr.shape[0]} samples, {expr.shape[1]} genes.')

# =============================================================================
# STEP 2 – DEFINE IMMUNE MARKER GENES (Canonical TIMER markers)
# =============================================================================
print('[2] Defining immune cell marker genes...')

# Based on TIMER2.0 canonical markers [citation:1]
immune_markers = {
    'B cell': 'CD19',
    'CD4+ T cell': 'CD4',
    'CD8+ T cell': 'CD8A',
    'Neutrophil': 'FCGR3B',
    'Macrophage': 'CD68',
    'Dendritic cell': 'ITGAX'
}

# =============================================================================
# STEP 3 – EXTRACT PHLDB3 AND MARKER EXPRESSION
# =============================================================================
print('[3] Extracting expression data...')

phldb3 = expr['PHLDB3']
phldb3.index = phldb3.index.str[:12]
phldb3 = phldb3[~phldb3.index.duplicated(keep='first')]

results = []

for cell_type, marker in immune_markers.items():
    if marker in expr.columns:
        marker_expr = expr[marker]
        marker_expr.index = marker_expr.index.str[:12]
        marker_expr = marker_expr[~marker_expr.index.duplicated(keep='first')]
        
        # Merge and calculate correlation
        merged = pd.DataFrame({'phldb3': phldb3, 'marker': marker_expr}).dropna()
        
        if len(merged) > 10:
            corr, p_val = spearmanr(merged['phldb3'], merged['marker'])
            results.append({
                'Cell_Type': cell_type,
                'Marker_Gene': marker,
                'Correlation': corr,
                'P_value': p_val,
                'N': len(merged)
            })
            print(f'    {cell_type} (n={len(merged)}): rho = {corr:.3f}, p = {p_val:.6f}')
        else:
            print(f'    ⚠️ {cell_type}: insufficient samples (n={len(merged)})')

# =============================================================================
# STEP 4 – CREATE RESULTS TABLE
# =============================================================================
print('[4] Creating results table...')

results_df = pd.DataFrame(results)
results_df['Significant'] = results_df['P_value'] < 0.05
results_df['Direction'] = results_df['Correlation'].apply(lambda x: 'Positive' if x > 0 else 'Negative')

print('\nResults summary:')
print(results_df.to_string(index=False))

# =============================================================================
# STEP 5 – BAR PLOT
# =============================================================================
print('[5] Generating bar plot...')

plt.figure(figsize=(10, 6))
colors = ['red' if p < 0.05 else 'gray' for p in results_df['P_value']]
plt.barh(results_df['Cell_Type'], results_df['Correlation'], color=colors)
plt.xlabel('Spearman Correlation with PHLDB3 Expression', fontsize=12)
plt.ylabel('Immune Cell Type', fontsize=12)
plt.title('PHLDB3 Association with Immune Cell Infiltration in PAAD', fontsize=14, fontweight='bold')
plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
plt.grid(alpha=0.3)

# Add p-value annotations
for i, row in results_df.iterrows():
    if row['P_value'] < 0.05:
        plt.text(0.02, i, f'p={row["P_value"]:.4f}*', 
                 fontsize=9, verticalalignment='center')

plt.tight_layout()
plt.savefig('PAAD_PHLDB3_Immune_Barplot.png', dpi=300)
print('    Bar plot saved: PAAD_PHLDB3_Immune_Barplot.png')

# =============================================================================
# STEP 6 – SAVE RESULTS
# =============================================================================
print('[6] Saving results...')

results_df.to_csv('PAAD_PHLDB3_Immune_Results.csv', index=False)
print('    Results saved: PAAD_PHLDB3_Immune_Results.csv')

print('\n✅ Immune infiltration analysis complete!')