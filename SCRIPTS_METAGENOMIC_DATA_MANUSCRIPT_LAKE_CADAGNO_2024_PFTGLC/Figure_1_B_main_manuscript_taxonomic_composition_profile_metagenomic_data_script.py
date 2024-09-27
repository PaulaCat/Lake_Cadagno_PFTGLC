import pandas as pd
import matplotlib.pyplot as plt

# Define the sample order and reverse it to display in reverse order on the y-axis
sample_order = ["3_cm", "40_cm", "153_cm", "187_cm", "213_cm", "233_cm", "283_cm", "382_cm", "532_cm", "566_cm", "582_cm", "693_cm", "738_cm"]
sample_order.reverse()

# Define the taxonomic groups with their respective colors in the order you specified
group_colors = {
    "Acidobacteria": "#80fed4",
    "Actinobacteria": "#238b21",
    "Candidatus Aminicenantes": "#528a8b",
    "Deltaproteobacteria": "#fe0000",
    "Other Proteobacteria": "#fe69b4",
    "Anaerolineae": "#8b0100",
    "Dehalococcoidia": "#bfff3d",
    "Other Chloroflexi": "#2f4f4f",
    "Candidatus Atribacteria": "#473c8b",
    "Bacteroidetes": "#00ffff",
    "Candidatus Bipolaricaulota": "#feff00",
    "Firmicutes": "#ffebcd",
    "Nitrospirae": "#000000",
    "Planctomycetes": "#ff8c00",
    "Candidatus Aenigmarchaeota": "#e5ce15",
    "Euryarchaeota": "#ba2452",
    "Candidatus Woesearchaeota": "#9400d3",
    "Candidatus Bathyarchaeia": "#181a6f",
    "Unclassified archaea": "#ffc0cb",
    "Others": "#bebebe"  # Default color for 'Others'
}

# Load the dataset
file_path = "MASTER_TABLE_METAGENOMIC_DATA_MANUSCRIPT_LAKE_CADAGNO_2024_PFTGLC.csv"
df = pd.read_csv(file_path, low_memory=False)

# Filter rows where 'Functional_annotation_short_blast_or_tree_placement' contains 'rps3' (case-insensitive)
df_filtered = df[df['Functional_annotation_short_blast_or_tree_placement'].str.contains('rps3', case=False, na=False)]

# Convert abundance columns to numeric, coercing errors to NaN
abundance_columns = sample_order
df_filtered[abundance_columns] = df_filtered[abundance_columns].apply(pd.to_numeric, errors='coerce')

# Replace NaN values with 0 in the abundance columns
df_filtered[abundance_columns] = df_filtered[abundance_columns].fillna(0)

# Sum abundance values for each taxonomic group
df_grouped = df_filtered.groupby('BLAST_taxonomic_annotation')[abundance_columns].sum()

# Create a new dataframe to store the groups from the color list and an "Others" category
df_selected = df_grouped.loc[df_grouped.index.intersection(group_colors.keys())]
df_others = df_grouped.loc[~df_grouped.index.isin(group_colors.keys())].sum().to_frame().T
df_others.index = ['Others']

# Combine the selected groups and "Others"
df_combined = pd.concat([df_selected, df_others])

# Sort the rows (taxa) based on the provided group order
df_combined = df_combined.reindex(group_colors.keys())

# Normalize the data by the total sum of each sample (column-wise normalization)
df_combined_normalized = df_combined.div(df_combined.sum(axis=0), axis=1)

# Get the colors in the correct order for the selected groups and "Others"
colors_in_plot = [group_colors[group] for group in df_combined_normalized.index]

# Plotting the stacked bar plot
fig, ax = plt.subplots(figsize=(10, 7))
df_combined_normalized.T.plot(kind='barh', stacked=True, color=colors_in_plot, ax=ax)

# Labeling the axes
ax.set_xlabel("Relative Abundance")
ax.set_ylabel("Sample Depth (cm)")
ax.set_title("Figure 1B (Main Manuscript)-Relative abundances of RpS3 protein sequences  with depth")

# Move the legend outside the plot
ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), title='Taxonomic Groups')

# Save the plot as PNG
plt.savefig("Figure_1_b_rps3_normalized_abundance_by_taxonomic_groups.pdf", bbox_inches='tight')

# Show plot
plt.show()
