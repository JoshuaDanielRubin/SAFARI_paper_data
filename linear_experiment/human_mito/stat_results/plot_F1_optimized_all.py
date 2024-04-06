import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np

def optimize_k_w(group):
    optimal_row = group.loc[group['sensitivity'].idxmax()]
    # Print the selected values of w and k
    print(f"Selected w: {optimal_row['k']}, k: {optimal_row['w']}")
    return optimal_row[['k', 'w']]

def main():
    file_path = sys.argv[1]
    data = pd.read_csv(file_path)

    # Find the optimal k and w for each tool
    optimal_params = data.groupby('tool').apply(optimize_k_w).reset_index()

    # Merge the original data with the optimal parameters to filter it
    optimal_data = pd.merge(data, optimal_params, on=['tool', 'k', 'w'])

    # Calculate median F1 scores
    median_f1_scores = optimal_data.groupby(['tool', 'damage_level'])['f1'].median().reset_index()

    # Define custom order for damage levels
    damage_order = ['None', 'Single-stranded', 'Mid', 'High']
    median_f1_scores['damage_level'] = pd.Categorical(median_f1_scores['damage_level'], categories=damage_order, ordered=True)

    # Extract tool order from the 'High' damage level
    high_damage_order = median_f1_scores[median_f1_scores['damage_level'] == 'High'].sort_values('f1', ascending=False)['tool'].tolist()

    # Ensure the tool column is Categorical with the order determined from the 'High' damage level
    median_f1_scores['tool'] = pd.Categorical(median_f1_scores['tool'], categories=high_damage_order, ordered=True)

    # Sort the DataFrame by the 'damage_level' to ensure the plot follows the custom order
    median_f1_scores.sort_values('damage_level', inplace=True)

    corrected_palette = {
    'vg giraffe': 'orange', 
    'SAFARI': 'green', 
    'BBMAP': '#1f7872',  # Teal
    'SHRiMP': '#6996b3',  # Light Blue
    'BWA-MEM': '#5e3a8c',  # Deep Violet
    'BWA ALN': '#ff99cc',  # Soft Magenta
    'BWA ALN (anc)': '#36454f',  # Charcoal Grey
    'Bowtie2': '#a45a52'  # Rust
    }

    # Create a barplot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=median_f1_scores, x='damage_level', y='f1', hue='tool', palette=corrected_palette)
    plt.title('Median F1 Score (Optimized Parameters, ' + sys.argv[3] + ")")
    plt.xlabel('Damage Level')
    plt.ylabel('Median F1 Score')
    plt.legend(title='Tool', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.ylim(0.90, 1)  # Adjust based on your needs
    plt.tight_layout()

    plt.savefig(sys.argv[2])

if __name__ == "__main__":
    main()

