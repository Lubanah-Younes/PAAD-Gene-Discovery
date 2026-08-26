# ============================================================================
# PAAD PROTEOMICS ANALYSIS – PHLDB3 (CPTAC)
# ============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: Analyzes PHLDB3 protein levels in CPTAC pancreatic cancer cohort.
# ============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

# ============================================================================
# STEP 1 – LOAD PROTEOME DATA
# ============================================================================
print('[1] Loading proteomics data...')

# Load tumor proteome
tumor = pd.read_csv('CPTAC_Proteome_Tumor.cct', sep='\t', index_col=0)
print(f'    Tumor proteome: {tumor.shape[0]} genes, {tumor.shape[1]} samples')

# Load normal proteome
normal = pd.read_csv('CPTAC_Proteome_Normal.cct', sep='\t', index_col=0)
print(f'    Normal proteome: {normal.shape[0]} genes, {normal.shape[1]} samples')

# ============================================================================
# STEP 2 – EXTRACT PHLDB3 PROTEIN LEVELS
# ============================================================================
print('[2] Extracting PHLDB3 protein levels...')

if 'PHLDB3' in tumor.index:
    phldb3_tumor = tumor.loc['PHLDB3']
    print(f'    PHLDB3 found in tumor samples (n={len(phldb3_tumor)})')
else:
    print('    ❌ PHLDB3 not found in tumor proteome')
    exit()

if 'PHLDB3' in normal.index:
    phldb3_normal = normal.loc['PHLDB3']
    print(f'    PHLDB3 found in normal samples (n={len(phldb3_normal)})')
else:
    print('    ❌ PHLDB3 not found in normal proteome')
    exit()

# ============================================================================
# STEP 3 – STATISTICAL TEST
# ============================================================================
print('[3] Statistical test...')

t_stat, p_val = ttest_ind(phldb3_tumor, phldb3_normal)
print(f'    T-test p-value: {p_val:.6f}')

# ============================================================================
# STEP 4 – BOXPLOT
# ============================================================================
print('[4] Generating boxplot...')

# Prepare data for plotting
import pandas as pd
plot_data = pd.DataFrame({
    'Protein Level': np.concatenate([phldb3_tumor.values, phldb3_normal.values]),
    'Group': ['Tumor']*len(phldb3_tumor) + ['Normal']*len(phldb3_normal)
})

plt.figure(figsize=(8, 6))
sns.boxplot(data=plot_data, x='Group', y='Protein Level', palette=['red', 'blue'])
plt.title('PHLDB3 Protein Levels in CPTAC-PAAD', fontsize=14, fontweight='bold')
plt.ylabel('Protein Intensity (log2)', fontsize=12)
plt.xlabel('', fontsize=12)
plt.text(0.02, 0.95, f'p = {p_val:.6f}', transform=plt.gca().transAxes, fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('PAAD_PHLDB3_Proteomics.png', dpi=300)
print('    Boxplot saved: PAAD_PHLDB3_Proteomics.png')

# ============================================================================
# STEP 5 – SUMMARY
# ============================================================================
print('\n[5] Summary:')
print(f'    Tumor mean: {phldb3_tumor.mean():.3f}')
print(f'    Normal mean: {phldb3_normal.mean():.3f}')
print(f'    Fold change: {phldb3_tumor.mean() / phldb3_normal.mean():.3f}')

print('\n✅ Proteomics analysis complete!')