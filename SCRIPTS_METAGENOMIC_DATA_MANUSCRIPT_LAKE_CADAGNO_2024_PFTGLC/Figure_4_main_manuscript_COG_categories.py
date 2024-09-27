import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset
df = pd.read_csv('MASTER_TABLE_METAGENOMIC_DATA_MANUSCRIPT_LAKE_CADAGNO_2024_PFTGLC.csv')

# Define the COG categories and their respective colors
cog_categories = {
    'E': ('Amino acid transport and metabolism', '#dbc750'),
    'G': ('Carbohydrate transport and metabolism', '#807a38'),
    'H': ('Coenzyme transport and metabolism', '#346566'),
    'V': ('Defense mechanisms', '#ffcd4b'),
    'C': ('Energy production and conversion', '#a14c08'),
    'P': ('Inorganic ion transport and metabolism', '#2ca7bc'),
    'I': ('Lipid transport and metabolism', '#fd9e55'),
    'F': ('Nucleotide transport and metabolism', '#ff6f0c'),
    'S': ('Poorly characterized', '#98c979'),
    'Q': ('Secondary metabolites biosynthesis, transport and catabolism', '#228e44')
}

supercategories = {
    'CELLULAR PROCESSES AND SIGNALING (D,Y,T,M,N,Z,W,U,O)': {
        'cogs': ['D', 'Y', 'T', 'M', 'N', 'Z', 'W', 'U', 'O'],
        'color': '#7fa69c'
    },
    'INFORMATION STORAGE AND PROCESSING (J,A,K,L,B)': {
        'cogs': ['J', 'A', 'K', 'L', 'B'],
        'color': '#76c9de'
    }
}

# Define sample columns
cm_columns = ['566_cm', '3_cm', '187_cm', '382_cm', '40_cm', '283_cm', '233_cm', 
              '153_cm', '532_cm', '738_cm', '582_cm', '693_cm', '213_cm']

# Filter rows where COG_category is not NA or '-'
valid_df = df[(df['COG_category'].notna()) & (df['COG_category'] != '-')]

# Initialize a DataFrame for holding the abundances by category and supercategory
normalized_abundances = pd.DataFrame(columns=cm_columns)

# Initialize abundances for each COG category and supercategory
for _, row in valid_df.iterrows():
    cog_categories_in_row = row['COG_category']
    for letter in cog_categories_in_row:
        # If the letter is a valid COG category
        if letter in cog_categories:
            if letter not in normalized_abundances.index:
                normalized_abundances.loc[letter] = 0
            normalized_abundances.loc[letter] += row[cm_columns].astype(float)
        # Check if letter belongs to any supercategory
        for supercat, data in supercategories.items():
            if letter in data['cogs']:
                if supercat not in normalized_abundances.index:
                    normalized_abundances.loc[supercat] = 0
                normalized_abundances.loc[supercat] += row[cm_columns].astype(float)

# Calculate total abundance per sample
total_abundance_per_sample = normalized_abundances.sum(axis=0)

# Normalize the abundances for each category and supercategory
normalized_abundance_df = normalized_abundances.div(total_abundance_per_sample, axis=1)

# Output the sum of abundances per sample per category to a CSV file before normalization
abundances_sum_df = normalized_abundances.sum(axis=1)
abundances_sum_df.to_csv('Figure_4_abundances_sum_per_sample_per_C0G_category.csv', header=["Total Abundance"])

# Check if the relative abundance sums to 1 for each sample
for sample in cm_columns:
    # Use round() to handle float values safely
    total_relative = np.round(normalized_abundance_df[sample].sum(), 6)
    
    if total_relative != 1:
        print(f"Warning: Total relative abundance for {sample} does not add up to 1. It's {total_relative}.")

# Define the correct order of samples for plotting
desired_order = ['3_cm', '40_cm', '153_cm', '187_cm', '213_cm', '233_cm', '283_cm', 
                 '382_cm', '532_cm', '566_cm', '582_cm', '693_cm', '738_cm']

# Ensure the DataFrame contains only the desired columns
available_columns = [col for col in desired_order if col in normalized_abundance_df.columns]

# Reorder the DataFrame to match the desired order
normalized_abundance_df = normalized_abundance_df[available_columns]

# Define the color list for the plot based on the category and supercategory
colors = []
legend_labels = []

for category in normalized_abundance_df.index:
    if category in cog_categories:
        colors.append(cog_categories[category][1])
        legend_labels.append(cog_categories[category][0])
    else:
        for supercat, data in supercategories.items():
            if category == supercat:
                colors.append(data['color'])
                legend_labels.append(supercat)

# Plot the normalized stacked bar plot with the correct sample order and colors
ax = normalized_abundance_df.T.plot(kind='barh', stacked=True, figsize=(10, 8), color=colors)

plt.xlabel('Relative Abundance')
plt.ylabel('Sample Depth (cm)')
plt.title('Normalized Abundance of COG Categories and Supercategories by Sample')
plt.legend(legend_labels, title='COG Categories', bbox_to_anchor=(1.05, 1), loc='upper left')

# Invert the y-axis to match the desired sample order
plt.gca().invert_yaxis()

# Save the plot to a PDF
plt.savefig('Figure_4_main_manuscript_COG_categories.pdf', bbox_inches='tight')

plt.show()
