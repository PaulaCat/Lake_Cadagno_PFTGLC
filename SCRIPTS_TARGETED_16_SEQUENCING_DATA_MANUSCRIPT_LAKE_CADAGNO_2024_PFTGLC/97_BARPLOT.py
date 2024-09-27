import pandas as pd
import re
import matplotlib.pyplot as plt

# Load the input OTU table
input_file_path = 'REAL0.97__complete_otu_table_qiime_Cara_Classification_97_identity.csv'
otu_df = pd.read_csv(input_file_path)

# Define the list of sample columns and ensure they follow the specified order
sample_order = ['3_cm', '40_cm', '153_cm', '187_cm', '213_cm', '233_cm', '283_cm', '382_cm', '532_cm', '566_cm', '582_cm', '693_cm', '738_cm']
sample_order = [s.replace("_cm", "cm") for s in sample_order]  # Adjust format to match the previous usage

# Strip whitespace from columns in case of extra spaces
otu_df.columns = otu_df.columns.str.strip()

# List of taxonomic groups and their corresponding matching strings
taxonomic_groups = {
    'Candidatus Bathyarchaeia': ['crenarchaeota'],
    'Acidobacteria': ['acidobacteria'],
    'Actinobacteria': ['actinobacteria'],  # Added Actinobacteria
    'Aminicenantes': ['aminicenantes'],
    'Alphaproteobacteria': ['alphaproteobacteria'],
    'Deltaproteobacteria': ['deltaproteobacteria'],
    'Gammaproteobacteria': ['gammaproteobacteria'],
    'Anaerolineae': ['anaerolineae'],
    'Dehalococcoidia': ['dehalococcoidia'],
    'Other Chloroflexi': ['chloroflexi'],  # Special case handled separately
    'Bacteroidetes': ['bacteroidetes'],
    'Caldiserica': ['caldiserica'],
    'Atribacteria': ['atribacteria'],
    'Cyanobacteria/Chloroplast': ['cyanobacteria', 'chloroplast'],
    'Firmicutes': ['firmicutes'],
    'Nitrospirae': ['nitrospirae'],
    'Planctomycetes': ['planctomycetes'],
    'Rhodothermaeota': ['rhodothermaeota'],
    'Spirochaetes': ['spirochaetes'],
    'Verrucomicrobia': ['verrucomicrobia'],
    'Woesearchaeota': ['woesearchaeota'],
    'Euryarchaeota': ['euryarchaeota'],
    'Candidatus Pacearchaeota': ['pacearchaeota']
}

# Custom colors for each taxonomic group using hex codes (reuse colors from the previous plot)
taxonomic_colors = {
    'Candidatus Bathyarchaeia': '#191970',
    'Acidobacteria': '#7fffd4',
    'Actinobacteria': '#228b22',  # Example color for Actinobacteria
    'Aminicenantes': '#528b8b',
    'Alphaproteobacteria': '#dda0dd',
    'Deltaproteobacteria': '#ff0000',
    'Gammaproteobacteria': '#ffc1c1',
    'Anaerolineae': '#8b0000',
    'Dehalococcoidia': '#c0ff3e',
    'Other Chloroflexi': '#2f4f4f',
    'Bacteroidetes': '#00ffff',
    'Caldiserica': '#c1ffc1',
    'Atribacteria': '#ff8c00',
    'Cyanobacteria/Chloroplast': '#ffff00',
    'Firmicutes': '#ffebcd',
    'Nitrospirae': '#000000',
    'Planctomycetes': '#473d8b',
    'Rhodothermaeota': '#20b2aa',
    'Spirochaetes': '#ff6a6a',
    'Verrucomicrobia': '#c6e2ff',
    'Woesearchaeota': '#9400d3',
    'Euryarchaeota': '#ba2552',
    'Candidatus Pacearchaeota': '#00ff00',  # Example color for Candidatus Pacearchaeota
    'Others': '#808080'
}

# Columns with sample data, ensuring they are in the correct order
sample_columns = [col for col in otu_df.columns if '_cm' in col]
sample_columns_ordered = sorted(sample_columns, key=lambda x: sample_order.index(x.replace("_cm", "cm")))

# Initialize a dictionary to store summed abundances for each group and sample
abundance_dict = {group: {sample: 0 for sample in sample_order} for group in taxonomic_groups}
abundance_dict['Others'] = {sample: 0 for sample in sample_order}

# Iterate over each row in the OTU table
for _, row in otu_df.iterrows():
    row_assigned = False  # Flag to check if the row has been assigned to a group
    
    # Check if the row matches any of the taxonomic groups
    for group, keywords in taxonomic_groups.items():
        for keyword in keywords:
            # Case-insensitive check using re.IGNORECASE
            pattern = re.compile(r'\b' + re.escape(keyword) + r'\b', re.IGNORECASE)
            if any(pattern.search(str(row[col])) for col in ['taxonomy', 'a', 'b', 'c', 'd']):
                # Special handling for "Other Chloroflexi"
                if group == 'Other Chloroflexi':
                    # Check for "Anaerolineae" or "Dehalococcoidia"
                    if any(re.search(r'\banaerolineae\b', str(row[col]), re.IGNORECASE) or 
                           re.search(r'\bdehalococcoidia\b', str(row[col]), re.IGNORECASE) 
                           for col in ['taxonomy', 'a', 'b', 'c', 'd']):
                        continue  # Skip this group if it has "Anaerolineae" or "Dehalococcoidia"
                
                # Add the abundance of this row to the respective group and sample columns
                for sample in sample_columns_ordered:
                    sample_name = sample.replace("_cm", "cm")
                    abundance_dict[group][sample_name] += row[sample]
                
                row_assigned = True
                break  # No need to check further keywords for this group
        if row_assigned:
            break  # No need to check further groups
    
    # If the row wasn't assigned to any group, add it to "Others"
    if not row_assigned:
        for sample in sample_columns_ordered:
            sample_name = sample.replace("_cm", "cm")
            abundance_dict['Others'][sample_name] += row[sample]

# Convert the abundance dictionary to a DataFrame
abundance_df = pd.DataFrame(abundance_dict).transpose()

# Reorder the DataFrame columns to match the specified order
abundance_df = abundance_df[sample_order]

# Save the aggregated abundance table to a CSV file
output_csv_path = 'OTU_97_aggregated_abundance.csv'
abundance_df.to_csv(output_csv_path)

# Normalize the abundance data to create a relative abundance plot
relative_abundance_df = abundance_df.div(abundance_df.sum(axis=0), axis=1)

# Plot settings
fig, ax = plt.subplots(figsize=(8, 10))

# Plotting the stacked bar chart
color_list = [taxonomic_colors[group] for group in relative_abundance_df.index]
relative_abundance_df.T.plot(kind='barh', stacked=True, ax=ax, color=color_list, edgecolor='none')

# Customize the plot
ax.set_xlabel('Relative Abundance')
ax.set_ylabel('Sample Depth (cm)')

# Set the correct y-tick positions and labels
ax.set_yticks(range(len(sample_order)))  # Set the y-ticks to match the number of samples
ax.set_yticklabels(sample_order)  # Set y-tick labels to sample_order

ax.invert_yaxis()  # To match the order in the image (top-down)
ax.set_title('Taxonomic Group Relative Abundance by Sample Depth')

# Create a legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, title='Taxonomic group', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Save the plot to a PDF
output_pdf_path = 'OTU_97_relative_abundance_plot.pdf'
plt.tight_layout()
plt.savefig(output_pdf_path, format='pdf')

print(f"Aggregated abundance table saved to {output_csv_path}")
print(f"Plot saved to {output_pdf_path}")
