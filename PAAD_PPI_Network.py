# =============================================================================
# PAAD PPI NETWORK ANALYSIS – PHLDB3
# =============================================================================
# Author: Lubanah Younes
# Date: August 2026
# Description: Builds a protein-protein interaction (PPI) network for PHLDB3
# using STRING data, identifies hub genes, and generates network plots.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns

# =============================================================================
# STEP 1 – LOAD STRING DATA
# =============================================================================
print('[1] Loading STRING data...')

string_file = 'STRING_PHLDB3_interactions.tsv'

try:
    string_data = pd.read_csv(string_file, sep='\t')
    print(f'    Loaded {string_data.shape[0]} interactions.')
except FileNotFoundError:
    print(f'    ❌ File "{string_file}" not found.')
    exit()

print(f'    Columns: {string_data.columns.tolist()}')

# =============================================================================
# STEP 2 – BUILD THE NETWORK
# =============================================================================
print('[2] Building network...')

# Identify column names (they vary slightly)
node1_col = 'node1' if 'node1' in string_data.columns else '#node1'
node2_col = 'node2' if 'node2' in string_data.columns else 'node2'
score_col = 'combined_score'

# Create graph
G = nx.Graph()
for _, row in string_data.iterrows():
    score = row[score_col]
    if score > 0.4:  # Medium confidence threshold
        G.add_edge(row[node1_col], row[node2_col], weight=score)

print(f'    Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

# =============================================================================
# STEP 3 – IDENTIFY HUB GENES
# =============================================================================
print('[3] Identifying hub genes...')

# Degree centrality
degree_centrality = nx.degree_centrality(G)
hub_genes = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)

print('    Top 10 hub genes:')
for gene, score in hub_genes[:10]:
    print(f'        {gene}: {score:.3f}')

# Save hub genes
hub_df = pd.DataFrame(hub_genes, columns=['Gene', 'Centrality'])
hub_df.to_csv('PAAD_PHLDB3_Hub_Genes.csv', index=False)
print('    Hub genes saved: PAAD_PHLDB3_Hub_Genes.csv')

# =============================================================================
# STEP 4 – NETWORK VISUALIZATION
# =============================================================================
print('[4] Visualizing network...')

plt.figure(figsize=(14, 12))

# Layout
pos = nx.spring_layout(G, k=0.8, iterations=100)

# Node sizes based on degree
node_sizes = [300 + degree_centrality[n] * 2000 for n in G.nodes()]

# Draw
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color='skyblue', alpha=0.8)
nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color='gray')
nx.draw_networkx_labels(G, pos, font_size=9)

plt.title('PPI Network for PHLDB3 (PAAD)', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('PAAD_PHLDB3_PPI_Network.png', dpi=300)
print('    Network plot saved: PAAD_PHLDB3_PPI_Network.png')

# =============================================================================
# STEP 5 – SUMMARY TABLE
# =============================================================================
print('[5] Summary statistics...')

summary = {
    'Total Nodes': G.number_of_nodes(),
    'Total Edges': G.number_of_edges(),
    'Average Degree': np.mean([d for n, d in G.degree()]),
    'Network Density': nx.density(G),
    'Top Hub Gene': hub_genes[0][0],
    'Top Hub Centrality': f'{hub_genes[0][1]:.3f}'
}

summary_df = pd.DataFrame([summary])
summary_df.to_csv('PAAD_PPI_Summary.csv', index=False)
print('    Summary saved: PAAD_PPI_Summary.csv')

print('\n✅ PPI network analysis complete!')
print('\nOutput files:')
print('  1. PAAD_PHLDB3_Hub_Genes.csv')
print('  2. PAAD_PHLDB3_PPI_Network.png')
print('  3. PAAD_PPI_Summary.csv')