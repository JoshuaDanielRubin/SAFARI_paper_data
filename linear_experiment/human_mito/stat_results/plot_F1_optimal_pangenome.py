import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():
    file_path = sys.argv[1]
    data = pd.read_csv(file_path)

    # Function to identify the optimal k and w for each tool
    def optimize_k_w(group):
        return group.loc[group['f1'].idxmax()][['k', 'w']]

    # Find the optimal k and w for each tool
    optimal_params = data.groupby('tool').apply(optimize_k_w).reset_index()

    # Merge the original data with the optimal parameters to filter it
    optimal_data = pd.merge(data, optimal_params, on=['tool', 'k', 'w'])

    # Subset the dataframe to only include specific tools
    subset_df = optimal_data[optimal_data['tool'].str.strip().str.lower().isin(['safari', 'vg giraffe'])]

    # Calculate median F1 scores
    median_f1_scores = subset_df.groupby(['tool', 'damage_level'])['f1'].median().reset_index()

    # Define custom order for damage levels
    damage_order = ['None', 'Single-stranded', 'Mid', 'High']
    median_f1_scores['damage_level'] = pd.Categorical(median_f1_scores['damage_level'], categories=damage_order, ordered=True)

    # Sort the DataFrame by the 'damage_level' to ensure the plot follows the custom order
    median_f1_scores.sort_values('damage_level', inplace=True)

    # Correcting the case sensitivity issue for the palette
    corrected_palette = {'vg giraffe': 'orange', 'SAFARI': 'green'}

    # Create a barplot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=median_f1_scores, x='damage_level', y='f1', hue='tool', palette=corrected_palette)
    plt.title('Median F1 Score by Pangenome Tool \n Stratified by Damage Level (Optimized Parameters)')
    plt.xlabel('Damage Level')
    plt.ylabel('Median F1 Score')
    plt.legend(title='Tool')
    plt.ylim(0.95, 1)  # Restrict y-axis to start at 0.7
    plt.tight_layout()

    plt.savefig(sys.argv[2])

if __name__ == "__main__":
    main()

