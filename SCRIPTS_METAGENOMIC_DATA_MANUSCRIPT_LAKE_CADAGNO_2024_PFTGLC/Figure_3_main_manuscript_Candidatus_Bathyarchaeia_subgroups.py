import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset
df = pd.read_csv('MASTER_TABLE_METAGENOMIC_DATA_MANUSCRIPT_LAKE_CADAGNO_2024_PFTGLC.csv', low_memory=False)

# Define sample columns
cm_columns = ['3_cm', '40_cm', '153_cm', '187_cm', '213_cm', '233_cm', '283_cm', '382_cm', '532_cm', '566_cm', '582_cm', '693_cm', '738_cm']

# Filter rows where the column 'Assigned_ID_Candidatus_Bathyarchaeia_cluster' is not NaN
valid_bathy_df = df[df['Assigned_ID_Candidatus_Bathyarchaeia_cluster'].notna()]

# Filter rows where the column 'Functional_annotation_short_blast_or_tree_placement' is RpS3
rps3_df = df[df['Functional_annotation_short_blast_or_tree_placement'] == 'RpS3']

# Calculate the total abundance of RpS3 per sample
total_rps3_per_sample = rps3_df[cm_columns].astype(float).sum()

# Calculate the abundance of each bathyarchaeia cluster per sample and normalize by the total RpS3 abundance
normalized_abundances = pd.DataFrame(columns=cm_columns)

for cluster in valid_bathy_df['Assigned_ID_Candidatus_Bathyarchaeia_cluster'].unique():
    cluster_df = valid_bathy_df[valid_bathy_df['Assigned_ID_Candidatus_Bathyarchaeia_cluster'] == cluster]
    cluster_abundance = cluster_df[cm_columns].astype(float).sum()
    normalized_abundances.loc[cluster] = cluster_abundance / total_rps3_per_sample

# Save the normalized abundances to a CSV file
normalized_abundances.to_csv('bathyarchaeia_cluster_summary_abundances.csv')

# Define the color map based on 'Taxonomic_assignment_Zhou_Hou_GDBTk_Ca_Bathyarchaeia'
color_map = {
    'Wuzhiqiibiales (Subgroup 15 - Ca. Bathyarchaeota B23)': '#ab8cc6',
    'Xuanwarculales (Subgroup 17 - Ca. Bathyarchaeota RGB_16_48_13)': '#99d96a',
    'Houtuarculales (Subgroup 18 - Ca. Bathyarchaeota 13_38_9)': '#824fcb',
    'Baizomonadales (Subgroup 6 - Ca. Bathyarchaeota 13_46_16b)': '#3d59ae',
    'Baizomonadales (Subgroup 13 - Ca. Bathyarchaeota B26-1)': '#e5ce53'
}

# Prepare the bubble plot
fig, ax = plt.subplots(figsize=(12, 8))

# Prepare clusters in the required order
cluster_order = ['1B', '1C', '1D', '2A', '2B', '3A', '3B', '3C', '3D', '4A', '4B', '4D', '4E', '4F', '4G', '4H', '5A', '5B', '5C', '5D', '5E', '5F', '5G', '5H']

# Plot bubbles for each cluster and sample
for i, cluster in enumerate(cluster_order):
    if cluster in normalized_abundances.index:
        for j, sample in enumerate(cm_columns):
            abundance = normalized_abundances.loc[cluster, sample]
            if not np.isnan(abundance) and abundance > 0:
                # Retrieve taxonomic group and handle missing or invalid values
                taxonomic_group = valid_bathy_df.loc[valid_bathy_df['Assigned_ID_Candidatus_Bathyarchaeia_cluster'] == cluster, 'Taxonomic_assignment_Zhou_Hou_GDBTk_Ca_Bathyarchaeia'].values[0]
                if pd.isna(taxonomic_group):
                    taxonomic_group = 'Unknown'
                if taxonomic_group in color_map:
                    # Adjust bubble size here
                    ax.scatter(i, j, s=abundance * 5000, color=color_map[taxonomic_group], alpha=0.6, edgecolor='black')
                else:
                    print(f"Warning: Taxonomic group '{taxonomic_group}' not found in color_map.")

# Adjust x and y axis
ax.set_xticks(np.arange(len(cluster_order)))
ax.set_xticklabels(cluster_order, rotation=90)

# Move the x-axis ticks and labels to the top
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 

# Adjust y-axis with sample names
ax.set_yticks(np.arange(len(cm_columns)))
ax.set_yticklabels(cm_columns)

# Invert the y-axis to match the desired sample order
ax.invert_yaxis()

# Adjust bubble size legend
bubble_sizes = [0.001, 0.01, 0.05, 0.1, 0.2, 0.43]
bubble_legend = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=np.sqrt(size * 5000), alpha=0.6) for size in bubble_sizes]
bubble_labels = [str(size) for size in bubble_sizes]

# Add the bubble size legend
bubble_legend = ax.legend(bubble_legend, bubble_labels, title="Normalized abundance", bbox_to_anchor=(1.05, 1), loc='upper left')

# Add color legend for taxonomic groups
color_legend = ax.legend([plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[element], markersize=10, alpha=0.6) for element in color_map], 
                         color_map.keys(), title="Phylogenetic group", bbox_to_anchor=(1.05, 0.6), loc='upper left')

# Add bubble legend back to the plot
ax.add_artist(bubble_legend)

# Save the plot to a PDF file
plt.savefig('Figure_3_manuscript_Candidatus_Bathyarchaeia_sugroups.pdf', bbox_inches='tight')

plt.show()
