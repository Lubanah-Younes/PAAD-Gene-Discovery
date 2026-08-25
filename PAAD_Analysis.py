# =============================================================================
# PAAD GENE DISCOVERY – SURVIVAL AND RISK ANALYSIS (FINAL VERSION)
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script performs:
#   1. Time‑dependent ROC curve for PHLDB3 (1, 3, 5 years)
#   2. Multivariate Cox regression using available top genes
#   3. Kaplan‑Meier plot for High‑risk vs Low‑risk groups
# All results are saved as PNG and CSV files.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from sklearn.metrics import roc_curve, auc

# =============================================================================
# STEP 1 – LOAD SURVIVAL DATA
# =============================================================================
surv = pd.read_csv('survival_PAAD_survival.txt', sep='\t')
surv = surv[['sample', 'OS', 'OS.time']].dropna()
surv.columns = ['sample', 'status', 'time']
print(f'[1] Loaded {len(surv)} patients.')

# =============================================================================
# STEP 2 – LOAD GENE EXPRESSION DATA
# =============================================================================
expr = pd.read_csv('TCGA.PAAD.sampleMap_HiSeqV2.gz',
                   compression='gzip', sep='\t', index_col=0)
expr = expr.T
print(f'[2] Loaded {expr.shape[0]} samples, {expr.shape[1]} genes.')

# =============================================================================
# STEP 3 – EXTRACT PHLDB3 AND MERGE
# =============================================================================
phldb3 = expr['PHLDB3']
data = surv.merge(phldb3, left_on='sample', right_index=True)
print(f'[3] Final dataset: {data.shape[0]} patients.')

# =============================================================================
# STEP 4 – COX REGRESSION FOR PHLDB3 (UNIVARIATE)
# =============================================================================
cph = CoxPHFitter()
cph.fit(data, duration_col='time', event_col='status', formula='PHLDB3')
print('\n[4] Cox model (PHLDB3):')
print(cph.summary)

# =============================================================================
# STEP 5 – TIME‑DEPENDENT ROC CURVE (1, 3, 5 YEARS) - FIXED
# =============================================================================
data['risk_score'] = np.exp(cph.params_['PHLDB3'] * data['PHLDB3'])

time_points = [365, 1095, 1825]
plt.figure(figsize=(9, 7))

for t in time_points:
    # ✅ المرضى اللي ماتوا خلال الفترة الزمنية
    y_true = ((data['time'] <= t) & (data['status'] == 1)).astype(int)
    y_pred = data['risk_score']
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, label=f'{int(t/365)}‑year AUC = {roc_auc:.3f}')

plt.plot([0, 1], [0, 1], 'k--', linewidth=1)
plt.xlabel('False Positive Rate', fontsize=12)
plt.ylabel('True Positive Rate', fontsize=12)
plt.title('Time‑dependent ROC for PHLDB3', fontsize=14, fontweight='bold')
plt.legend(loc='lower right')
plt.tight_layout()
plt.savefig('PAAD_PHLDB3_ROC.png', dpi=300)
print('[5] ROC curve saved: PAAD_PHLDB3_ROC.png')

# =============================================================================
# STEP 6 – MULTIVARIATE COX MODEL (AVAILABLE GENES ONLY)
# =============================================================================
top_genes = [
    'PHLDB3', 'SLURP1', 'MYEOV', 'USP20', 'LOC651250',
    'DEF8', 'EPS8', 'NCAM1', 'FAM123A', 'MYO5B'
]

# ✅ Check which genes actually exist in the expression data
available = []
for gene in top_genes:
    if gene in expr.columns:
        available.append(gene)
    else:
        print(f'⚠️ Gene "{gene}" not found – skipping.')

print(f'\n[6] Genes available for multivariate model: {available}')

# ✅ Proceed only if we have at least 2 genes
if len(available) >= 2:
    # ✅ Merge all available genes into data
    for gene in available:
        data[gene] = expr[gene]
    
    # ✅ Build formula from available genes
    formula = ' + '.join(available)
    print(f'Formula: {formula}')
    
    # ✅ Fit multivariate Cox model
    cph_multi = CoxPHFitter()
    cph_multi.fit(data, duration_col='time', event_col='status', formula=formula)
    print('\nMultivariate Cox model (summary):')
    print(cph_multi.summary)

    # ✅ Compute risk score and split by median
    data['multi_risk'] = np.exp(cph_multi.predict_partial_hazard(data))
    median_risk = data['multi_risk'].median()
    data['risk_group'] = np.where(data['multi_risk'] > median_risk,
                                  'High Risk', 'Low Risk')

    # ✅ Kaplan‑Meier plot for risk groups
    plt.figure(figsize=(9, 7))
    kmf = KaplanMeierFitter()
    for group in ['High Risk', 'Low Risk']:
        subset = data[data['risk_group'] == group]
        if len(subset) > 0:
            kmf.fit(subset['time'], event_observed=subset['status'],
                    label=f'{group} (n={len(subset)})')
            kmf.plot_survival_function()
    plt.title('Survival by Multivariate Risk Score (Top Genes)',
              fontsize=14, fontweight='bold')
    plt.xlabel('Time (Days)', fontsize=12)
    plt.ylabel('Survival Probability', fontsize=12)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('PAAD_Risk_Score_KM.png', dpi=300)
    print('[7] KM plot saved: PAAD_Risk_Score_KM.png')
else:
    print('❌ Not enough genes available for multivariate model – skipping.')

# =============================================================================
# STEP 7 – SAVE FULL RESULTS
# =============================================================================
data.to_csv('PAAD_analysis_results.csv', index=False)
print('\n[8] Results saved: PAAD_analysis_results.csv')
print('\n✅ Analysis complete – all outputs generated.')