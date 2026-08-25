# =============================================================================
# PAAD STAGE-SPECIFIC ROC ANALYSIS – PHLDB3 (FIXED)
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script compares PHLDB3 ROC performance between
# Early and Advanced stage groups.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

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

# Load survival data
surv = pd.read_csv('survival_PAAD_survival.txt', sep='\t')
surv['sample_id'] = surv['sample'].str[:12]
surv = surv[~surv['sample_id'].duplicated(keep='first')]

# =============================================================================
# STEP 2 – CREATE STAGE GROUPS
# =============================================================================
print('[2] Creating stage groups...')

median_os = surv['OS.time'].median()
surv['stage_group'] = np.where(surv['OS.time'] < median_os, 'Advanced', 'Early')

# Print group sizes
print(f'    Early group: {sum(surv["stage_group"] == "Early")} patients')
print(f'    Advanced group: {sum(surv["stage_group"] == "Advanced")} patients')

# =============================================================================
# STEP 3 – EXTRACT PHLDB3 AND MERGE
# =============================================================================
print('[3] Extracting PHLDB3 expression...')

phldb3 = expr['PHLDB3']
merged = pd.DataFrame({
    'expression': phldb3,
    'status': surv.set_index('sample_id')['OS'],
    'time': surv.set_index('sample_id')['OS.time'],
    'stage': surv.set_index('sample_id')['stage_group']
}).dropna()

print(f'    Merged data: {merged.shape[0]} samples')

# =============================================================================
# STEP 4 – CHECK GROUP SIZES AFTER MERGE
# =============================================================================
print('[4] Group sizes after merge:')
for group in ['Early', 'Advanced']:
    size = sum(merged['stage'] == group)
    print(f'    {group}: {size} samples')

# =============================================================================
# STEP 5 – ROC FOR EACH GROUP (1-year survival)
# =============================================================================
print('[5] Generating ROC curves...')

time_point = 365  # 1 year

plt.figure(figsize=(8, 6))

# Check if we have enough samples in each group
groups_with_data = []
for group in ['Early', 'Advanced']:
    subset = merged[merged['stage'] == group]
    if len(subset) >= 10:
        groups_with_data.append(group)
        y_true = ((subset['time'] > time_point) & (subset['status'] == 0)).astype(int)
        y_pred = subset['expression']
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f'{group} (n={len(subset)}, AUC = {roc_auc:.3f})')
        print(f'    {group}: AUC = {roc_auc:.3f}')
    else:
        print(f'    ⚠️ {group} group too small (n={len(subset)}) – skipping')

if len(groups_with_data) == 0:
    print('    ❌ Not enough data for any group.')
    exit()

plt.plot([0, 1], [0, 1], 'k--', linewidth=1)
plt.xlabel('False Positive Rate', fontsize=12)
plt.ylabel('True Positive Rate', fontsize=12)
plt.title('PHLDB3 ROC: Early vs Advanced (1-year survival)', fontsize=14, fontweight='bold')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('PAAD_PHLDB3_Stage_ROC.png', dpi=300)
print('    ROC plot saved: PAAD_PHLDB3_Stage_ROC.png')

# =============================================================================
# STEP 6 – ALSO CHECK DIFFERENT TIME POINTS
# =============================================================================
print('[6] Checking other time points...')

for t in [365, 1095, 1825]:
    print(f'\n    {t/365:.0f}-year survival:')
    for group in ['Early', 'Advanced']:
        subset = merged[merged['stage'] == group]
        if len(subset) >= 10:
            y_true = ((subset['time'] > t) & (subset['status'] == 0)).astype(int)
            y_pred = subset['expression']
            fpr, tpr, _ = roc_curve(y_true, y_pred)
            roc_auc = auc(fpr, tpr)
            print(f'        {group}: AUC = {roc_auc:.3f}')

print('\n✅ Stage-specific ROC analysis complete!')