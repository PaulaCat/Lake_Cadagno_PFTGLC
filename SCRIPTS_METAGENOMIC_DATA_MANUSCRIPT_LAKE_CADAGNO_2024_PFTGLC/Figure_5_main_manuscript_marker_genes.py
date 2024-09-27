
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset
df = pd.read_csv('MASTER_TABLE_METAGENOMIC_DATA_MANUSCRIPT_LAKE_CADAGNO_2024_PFTGLC.csv')

# Define the COG categories and elements with their respective colors
cog_categories = {
    'AcsA': 'Carbon',
    'CdhC': 'Carbon',
    'RuBisCO_Type_I': 'Carbon',
    'RuBisCO_Type_III': 'Carbon',
    'McrA': 'Carbon',
    'DsrA': 'Sulfur',
    'DsrB': 'Sulfur',
    'AprA': 'Sulfur',
    'AprB': 'Sulfur',
    'NarG': 'Nitrogen',
    'NarH': 'Nitrogen',
    'NifD': 'Nitrogen',
    'NifH': 'Nitrogen',
    'NifK': 'Nitrogen',
    'NfrA': 'Nitrogen',
    'NfrH': 'Nitrogen',
    'NirD': 'Nitrogen',
    'NorB': 'Nitrogen',
    'NorC': 'Nitrogen',
    'HAO': 'Nitrogen'
}

# Define sample columns
cm_columns = ['3_cm', '40_cm', '153_cm', '187_cm', '213_cm', '233_cm', '283_cm', '382_cm', '532_cm', '566_cm', '582_cm', '693_cm', '738_cm']

# Filter rows where the column 'Functional_annotation_short_blast_or_tree_placement' contains marker genes
valid_marker_df = df[df['Functional_annotation_short_blast_or_tree_placement'].isin(cog_categories.keys())]

# Filter rows where the column 'Functional_annotation_short_blast_or_tree_placement' is RpS3
rps3_df = df[df['Functional_annotation_short_blast_or_tree_placement'] == 'RpS3']

# Calculate the total abundance of RpS3 per sample
total_rps3_per_sample = rps3_df[cm_columns].astype(float).sum()

# Calculate the abundance of each marker gene per sample and normalize by the total RpS3 abundance
normalized_abundances = pd.DataFrame(columns=cm_columns)

for gene in cog_categories.keys():
    gene_df = valid_marker_df[valid_marker_df['Functional_annotation_short_blast_or_tree_placement'] == gene]
    gene_abundance = gene_df[cm_columns].astype(float).sum()
    normalized_abundances.loc[gene] = gene_abundance / total_rps3_per_sample

# Round values to 6 significant figures
normalized_abundances = normalized_abundances.round(6)

# Save the normalized abundances to a CSV file
normalized_abundances.to_csv('Figure_5_summarized_data_normalized_abundances_marker_genes.csv')

# Prepare the bubble plot
fig, ax = plt.subplots(figsize=(12, 8))

# Define the color map for elements
color_map = {'Carbon': '#fbbd65', 'Sulfur': '#d88c8c', 'Nitrogen': '#b8d39b'}

# Prepare marker genes in the required order
marker_genes_order = ['AcsA', 'CdhC', 'RuBisCO_Type_I', 'RuBisCO_Type_III', 'McrA', 'DsrA', 'DsrB', 'AprA', 'AprB', 
                      'NarG', 'NarH', 'NifD', 'NifH', 'NifK', 'NfrA', 'NfrH', 'NirD', 'NorB', 'NorC', 'HAO']

# Plot bubbles for each marker gene and sample
for i, gene in enumerate(marker_genes_order):
    for j, sample in enumerate(cm_columns):
        abundance = normalized_abundances.loc[gene, sample]
        if not np.isnan(abundance) and abundance > 0:
            ax.scatter(i, j, s=abundance * 800, color=color_map[cog_categories[gene]], alpha=0.6, edgecolor='black')

# Adjust x and y axis
ax.set_xticks(np.arange(len(marker_genes_order)))
ax.set_xticklabels(marker_genes_order, rotation=90)

# Move the x-axis ticks and labels to the top
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 

# Adjust y-axis with sample names
ax.set_yticks(np.arange(len(cm_columns)))
ax.set_yticklabels(cm_columns)

# Invert the y-axis to match the desired sample order
ax.invert_yaxis()

# Adjust bubble size legend
bubble_sizes = [0.05, 0.1, 0.2, 0.3, 0.6, 0.8]
bubble_legend = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=np.sqrt(size * 800), alpha=0.6) for size in bubble_sizes]
bubble_labels = [str(size) for size in bubble_sizes]

# Add the bubble size legend
bubble_legend = ax.legend(bubble_legend, bubble_labels, title="Normalized marker gene abundance", bbox_to_anchor=(1.05, 1), loc='upper left')

# Add color legend for elements
color_legend = ax.legend([plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[element], markersize=10, alpha=0.6) for element in color_map], 
                         color_map.keys(), title="Element", bbox_to_anchor=(1.05, 0.6), loc='upper left')

# Add bubble legend back to the plot
ax.add_artist(bubble_legend)

# Save the plot to a PDF file
plt.savefig('Figure_5_manuscript__marker_genes.pdf', bbox_inches='tight')

plt.show()
