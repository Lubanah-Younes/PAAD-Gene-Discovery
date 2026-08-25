# =============================================================================
# PAAD PATHWAY ENRICHMENT ANALYSIS – PHLDB3 INTERACTOME
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script performs pathway enrichment analysis using Enrichr
# for the PHLDB3 interactome genes.
# =============================================================================

import pandas as pd
import requests
import time

# =============================================================================
# STEP 1 – DEFINE GENE LIST
# =============================================================================
print('[1] Defining gene list...')

genes = [
    'LYPD3', 'TEX101', 'PHLDB3', 'SPACA4', 'PATE3', 'LY6L',
    'SLURP1', 'PATE1', 'PATE2', 'PINLYP', 'LY6G6D', 'ZNF575',
    'CRHBP', 'ETHE1', 'LYSMD1', 'SCNM1', 'AGR3', 'UCN2',
    'MDM2', 'RAB30'
]

print(f'    Total genes: {len(genes)}')

# =============================================================================
# STEP 2 – QUERY ENRICHR API
# =============================================================================
print('[2] Querying Enrichr API...')

gene_string = '\n'.join(genes)
ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'

files = {'file': ('genes.txt', gene_string, 'text/plain')}
data = {'description': 'PHLDB3_interactome'}

response = requests.post(ENRICHR_URL, files=files, data=data)
if response.status_code == 200:
    data = response.json()
    user_list_id = data.get('userListId')
    print(f'    Gene list submitted. ID: {user_list_id}')
else:
    print(f'    ❌ Error: {response.status_code}')
    exit()

# =============================================================================
# STEP 3 – FETCH ENRICHMENT RESULTS
# =============================================================================
print('[3] Fetching enrichment results...')

libraries = ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022']
results = {}

for lib in libraries:
    url = f'https://maayanlab.cloud/Enrichr/enrich?userListId={user_list_id}&backgroundType={lib}'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if lib in data and data[lib]:
            results[lib] = data[lib]
            print(f'    {lib}: {len(data[lib])} terms found')
        else:
            print(f'    ⚠️ {lib}: no data')
    else:
        print(f'    ❌ Error fetching {lib}')
    time.sleep(1)

# =============================================================================
# STEP 4 – PROCESS RESULTS (FIXED)
# =============================================================================
print('[4] Processing results...')

top_terms = {}

for lib, terms in results.items():
    if terms:
        # Enrichr returns 9 columns: Rank, Term, P-value, Z-score, Combined Score,
        # Genes, Adjusted P-value, Old P-value, Old Adjusted P-value
        df = pd.DataFrame(terms, columns=[
            'Rank', 'Term', 'P-value', 'Z-score', 'Combined Score',
            'Genes', 'Adjusted P-value', 'Old P-value', 'Old Adjusted P-value'
        ])
        df_sig = df[df['Adjusted P-value'] < 0.05].head(10)
        top_terms[lib] = df_sig
        print(f'\n    Top terms for {lib}:')
        if len(df_sig) > 0:
            for _, row in df_sig.iterrows():
                print(f'        {row["Term"]} (p = {row["Adjusted P-value"]:.4f})')
        else:
            print('        No significant terms found.')

# =============================================================================
# STEP 5 – SAVE RESULTS
# =============================================================================
print('[5] Saving results...')

for lib, df in top_terms.items():
    if not df.empty:
        filename = f'PAAD_PPI_{lib}_Enrichment.csv'
        df.to_csv(filename, index=False)
        print(f'    Saved: {filename}')

print('\n✅ Pathway enrichment analysis complete!')