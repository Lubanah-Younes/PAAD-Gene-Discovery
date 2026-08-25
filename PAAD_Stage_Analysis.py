# =============================================================================
# PAAD STAGE ANALYSIS – PHLDB3 EXPRESSION BY CLINICAL STAGE
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script compares PHLDB3 expression across different
# clinical stages using survival data.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# =============================================================================
# STEP 1 – LOAD DATA
# =============================================================================
print('[1] Loading data...')

# Load expression data
expr = pd.read_csv('TCGA.PAAD.sampleMap_HiSeqV2.gz',
                   compression='gzip', sep='\t', index_col=0)
expr = expr.T
print(f'    Expression data: {expr.shape[0]} samples, {expr.shape[1]} genes.')

# =============================================================================
# STEP 2 – EXTRACT PHLDB3 EXPRESSION
# =============================================================================
print('[2] Extracting PHLDB3 expression...')
phldb3_expr = expr['PHLDB3']
phldb3_expr.index = phldb3_expr.index.str[:12]

# ✅ Remove duplicate indices (keep first)
phldb3_expr = phldb3_expr[~phldb3_expr.index.duplicated(keep='first')]
print(f'    PHLDB3 expression: {phldb3_expr.shape[0]} samples')

# =============================================================================
# STEP 3 – CREATE STAGE FROM SURVIVAL DATA
# =============================================================================
print('[3] Creating stage groups from survival data...')

# Load survival data
surv = pd.read_csv('survival_PAAD_survival.txt', sep='\t')
surv['sample_id'] = surv['sample'].str[:12]

# ✅ Remove duplicate sample IDs
surv = surv[~surv['sample_id'].duplicated(keep='first')]

# Calculate median survival
median_os = surv['OS.time'].median()
print(f'    Median survival time: {median_os} days')

# Create stage groups based on survival
surv['stage_group'] = np.where(surv['OS.time'] < median_os, 'Advanced', 'Early')

# =============================================================================
# STEP 4 – MERGE EXPRESSION WITH STAGE
# =============================================================================
print('[4] Merging data...')

merged = pd.DataFrame({
    'expression': phldb3_expr,
    'stage': surv.set_index('sample_id')['stage_group']
}).dropna()

print(f'    Merged data: {merged.shape[0]} samples')
print(f'    Stage distribution: {merged["stage"].value_counts().to_dict()}')

# =============================================================================
# STEP 5 – STATISTICAL TEST (Mann-Whitney U)
# =============================================================================
print('[5] Statistical test...')

early_expr = merged[merged['stage'] == 'Early']['expression']
advanced_expr = merged[merged['stage'] == 'Advanced']['expression']

stat, p_mann = mannwhitneyu(early_expr, advanced_expr)
print(f'    Mann-Whitney U p-value: {p_mann:.6f}')

# =============================================================================
# STEP 6 – BOXPLOT
# =============================================================================
print('[6] Generating boxplot...')

plt.figure(figsize=(8, 6))
sns.boxplot(data=merged, x='stage', y='expression', palette=['blue', 'red'])
plt.xlabel('Survival-Based Stage Group', fontsize=12)
plt.ylabel('PHLDB3 Expression (log2)', fontsize=12)
plt.title('PHLDB3 Expression by Survival-Based Stage in PAAD', fontsize=14, fontweight='bold')
plt.grid(alpha=0.3)
plt.text(0.02, 0.95, f'Mann-Whitney U p = {p_mann:.4f}', 
         transform=plt.gca().transAxes, fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('PAAD_PHLDB3_Stage_Boxplot.png', dpi=300)
print('    Boxplot saved: PAAD_PHLDB3_Stage_Boxplot.png')

# =============================================================================
# STEP 7 – SUMMARY TABLE
# =============================================================================
print('[7] Summary by stage group...')

summary = merged.groupby('stage')['expression'].agg(['mean', 'std', 'count'])
summary.columns = ['Mean Expression', 'Std Dev', 'N']
print(summary)

# =============================================================================
# STEP 8 – SAVE RESULTS
# =============================================================================
print('[8] Saving results...')

merged.to_csv('PAAD_PHLDB3_Stage_results.csv')
summary.to_csv('PAAD_PHLDB3_Stage_summary.csv')

print('\n✅ Stage analysis complete!')
print('\nOutput files:')
print('  1. PAAD_PHLDB3_Stage_Boxplot.png')
print('  2. PAAD_PHLDB3_Stage_results.csv')
print('  3. PAAD_PHLDB3_Stage_summary.csv')