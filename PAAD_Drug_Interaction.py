# =============================================================================
# PAAD DRUG-GENE INTERACTION ANALYSIS – PHLDB3 INTERACTOME
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: This script queries DGIdb for drugs targeting the
# PHLDB3 interactome genes using the updated API.
# =============================================================================

import pandas as pd
import requests
import time
import json

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
# STEP 2 – QUERY DGIdb API (UPDATED)
# =============================================================================
print('[2] Querying DGIdb API...')

DGIdb_URL = 'https://dgidb.org/api/v2/interactions.json'
drug_results = []

for gene in genes:
    params = {'genes': gene}
    try:
        response = requests.get(DGIdb_URL, params=params, timeout=30)
        print(f'    Checking {gene}: status {response.status_code}')
        
        if response.status_code == 200:
            try:
                data = response.json()
                if data.get('interactions'):
                    for interaction in data['interactions']:
                        drug_results.append({
                            'Gene': gene,
                            'Drug': interaction.get('drugName', 'Unknown'),
                            'Interaction Type': interaction.get('interactionType', 'Unknown'),
                            'Source': interaction.get('source', 'Unknown')
                        })
                        print(f'        Found: {gene} ↔ {interaction.get("drugName", "Unknown")}')
            except json.JSONDecodeError:
                print(f'        ⚠️ {gene}: Invalid JSON response')
        elif response.status_code == 404:
            print(f'        ⚠️ {gene}: No interactions found')
        else:
            print(f'        ⚠️ {gene}: Error {response.status_code}')
    except requests.exceptions.Timeout:
        print(f'        ⚠️ {gene}: Request timeout')
    except Exception as e:
        print(f'        ⚠️ {gene}: {str(e)[:50]}')
    
    time.sleep(0.5)  # Be polite to the API

print(f'\n[3] Found {len(drug_results)} drug-gene interactions')

# =============================================================================
# STEP 3 – SAVE RESULTS
# =============================================================================
print('[4] Saving results...')

if drug_results:
    drug_df = pd.DataFrame(drug_results)
    drug_df.to_csv('PAAD_PPI_Drug_Interactions.csv', index=False)
    print('    Saved: PAAD_PPI_Drug_Interactions.csv')
    
    print('\n    Top drugs:')
    print(drug_df['Drug'].value_counts().head(10))
    
    print('\n    Summary by gene:')
    print(drug_df['Gene'].value_counts())
else:
    print('    No drug interactions found for these genes.')
    print('    This is normal – many genes are not yet targeted by known drugs.')

print('\n✅ Drug-gene interaction analysis complete!')