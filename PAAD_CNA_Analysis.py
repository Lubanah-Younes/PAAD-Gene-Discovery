# =============================================================================
# PAAD CNA AND EXPRESSION ANALYSIS – MULTI-OMICS INTEGRATION
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script integrates copy number alteration (CNA) data
# with gene expression data for PHLDB3 in TCGA-PAAD cohort.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, spearmanr
import glob

# =============================================================================
# STEP 1 – LOAD EXPRESSION DATA
# =============================================================================
print('[1] Loading expression data...')
expr = pd.read_csv('TCGA.PAAD.sampleMap_HiSeqV2.gz',
                   compression='gzip', sep='\t', index_col=0)
expr = expr.T  # samples as rows
print(f'    Expression data: {expr.shape[0]} samples, {expr.shape[1]} genes.')

# =============================================================================
# STEP 2 – LOAD CNA DATA
# =============================================================================
print('[2] Looking for CNA file...')
cct_files = glob.glob('*.cct') + glob.glob('*.cct.gz')
if len(cct_files) == 0:
    print('    ❌ No .cct file found.')
    exit()

cna_file = cct_files[0]
print(f'    Found CNA file: {cna_file}')

if cna_file.endswith('.gz'):
    cna = pd.read_csv(cna_file, compression='gzip', sep='\t', index_col=0)
else:
    cna = pd.read_csv(cna_file, sep='\t', index_col=0)
print(f'    Loaded CNA data: {cna.shape[0]} genes, {cna.shape[1]} samples')

# =============================================================================
# STEP 3 – FIX SAMPLE NAMES
# =============================================================================
print('[3] Fixing sample names...')

# CNA names: convert "TCGA.2J.AAB1" → "TCGA-2J-AAB1"
cna.columns = [col.replace('.', '-') for col in cna.columns]

# Expression names: extract first 12 characters (TCGA-XX-XXXX)
expr.index = [idx[:12] for idx in expr.index]

print(f'    CNA samples (first 5): {cna.columns[:5].tolist()}')
print(f'    Expression samples (first 5): {expr.index[:5].tolist()}')

# =============================================================================
# STEP 4 – REMOVE DUPLICATES
# =============================================================================
print('[4] Removing duplicates...')

# Remove duplicate samples from expression data (keep first occurrence)
expr = expr[~expr.index.duplicated(keep='first')]
print(f'    Expression data after removing duplicates: {expr.shape[0]} samples')

# Remove duplicate samples from CNA data
cna = cna.loc[:, ~cna.columns.duplicated(keep='first')]
print(f'    CNA data after removing duplicates: {cna.shape[1]} samples')

# =============================================================================
# STEP 5 – EXTRACT PHLDB3
# =============================================================================
print('[5] Extracting PHLDB3 data...')

if 'PHLDB3' in cna.index:
    phldb3_cna = cna.loc['PHLDB3']
    print(f'    PHLDB3 found in CNA data ({len(phldb3_cna)} samples)')
else:
    print('    ⚠️ PHLDB3 not found in CNA data.')
    exit()

if 'PHLDB3' in expr.columns:
    phldb3_expr = expr['PHLDB3']
    print(f'    PHLDB3 expression: {phldb3_expr.shape[0]} samples')
else:
    print('    ❌ PHLDB3 not found in expression data.')
    exit()

# =============================================================================
# STEP 6 – MERGE DATA
# =============================================================================
print('[6] Merging data...')

# Convert to Series
phldb3_expr_series = phldb3_expr
phldb3_cna_series = phldb3_cna

# Find common samples
common_samples = set(phldb3_expr_series.index) & set(phldb3_cna_series.index)
print(f'    Common samples: {len(common_samples)}')

if len(common_samples) == 0:
    print('    ❌ No common samples found.')
    exit()

# Create merged DataFrame
merged = pd.DataFrame({
    'expression': phldb3_expr_series.loc[list(common_samples)],
    'cna': phldb3_cna_series.loc[list(common_samples)]
})

print(f'    Merged data: {merged.shape[0]} samples')

# =============================================================================
# STEP 7 – STATISTICAL TESTS
# =============================================================================
print('[7] Statistical tests...')

if merged.shape[0] > 1:
    # Spearman correlation
    corr, p_spearman = spearmanr(merged['expression'], merged['cna'])
    print(f'    Spearman correlation: rho = {corr:.3f}, p = {p_spearman:.6f}')
    
    # Kruskal-Wallis test
    groups = [group['expression'].values for name, group in merged.groupby('cna')]
    if len(groups) > 1:
        stat, p_kruskal = kruskal(*groups)
        print(f'    Kruskal-Wallis p-value: {p_kruskal:.6f}')
else:
    print('    ❌ Not enough data for statistical tests.')

# =============================================================================
# STEP 8 – PLOTS
# =============================================================================
print('[8] Generating plots...')

if merged.shape[0] > 1:
    # Boxplot
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=merged, x='cna', y='expression', palette='Set2')
    plt.xlabel('PHLDB3 Copy Number Alteration', fontsize=12)
    plt.ylabel('PHLDB3 Expression (log2)', fontsize=12)
    plt.title('PHLDB3 Expression by Copy Number Status in PAAD', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('PAAD_PHLDB3_CNA_Boxplot.png', dpi=300)
    print('    Boxplot saved: PAAD_PHLDB3_CNA_Boxplot.png')
    
    # Violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(data=merged, x='cna', y='expression', palette='Set2')
    plt.xlabel('PHLDB3 Copy Number Alteration', fontsize=12)
    plt.ylabel('PHLDB3 Expression (log2)', fontsize=12)
    plt.title('PHLDB3 Expression by CNA Status (Violin Plot)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('PAAD_PHLDB3_CNA_Violin.png', dpi=300)
    print('    Violin plot saved: PAAD_PHLDB3_CNA_Violin.png')

# =============================================================================
# STEP 9 – SAVE RESULTS
# =============================================================================
print('[9] Saving results...')
merged.to_csv('PAAD_PHLDB3_CNA_results.csv')
print('    Results saved: PAAD_PHLDB3_CNA_results.csv')

print('\n✅ Multi‑omics analysis complete!')