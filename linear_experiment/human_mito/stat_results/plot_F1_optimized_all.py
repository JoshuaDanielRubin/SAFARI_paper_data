import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():
    file_path = sys.argv[1]
    data = pd.read_csv(file_path)

    def optimize_k_w(group):
        optimal_row = group.loc[group['sensitivity'].idxmax()]
        # Print the selected values of w and k
        print(f"Selected w: {optimal_row['w']}, k: {optimal_row['k']}")
        return optimal_row[['k', 'w']]

    # Find the optimal k and w for each tool
    optimal_params = data.groupby('tool').apply(optimize_k_w).reset_index()

    # Merge the original data with the optimal parameters to filter it
    optimal_data = pd.merge(data, optimal_params, on=['tool', 'k', 'w'])

    # Calculate median F1 scores
    median_f1_scores = optimal_data.groupby(['tool', 'damage_level'])['f1'].mean().reset_index()

    # Define custom order for damage levels
    damage_order = ['None', 'Single-stranded', 'Mid', 'High']
    median_f1_scores['damage_level'] = pd.Categorical(median_f1_scores['damage_level'], categories=damage_order, ordered=True)

    # Sort the DataFrame by the 'damage_level' to ensure the plot follows the custom order
    median_f1_scores.sort_values('damage_level', inplace=True)

    print(median_f1_scores)

    # Correcting the case sensitivity issue for the palette
    corrected_palette = {'vg giraffe': 'orange', 'SAFARI': 'green', 'BBMAP': 'blue', 'SHRiMP': 'red', 'BWA-MEM':'purple', \
                         'BWA ALN': 'pink', 'BWA ALN (anc)':'black', 'Bowtie2': 'brown'}

    # Create a barplot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=median_f1_scores, x='damage_level', y='f1', hue='tool', palette=corrected_palette)
    plt.title('Median F1 Score by Pangenome Tool \n Stratified by Damage Level (Optimized Parameters, ' + sys.argv[3] + ")")
    plt.xlabel('Damage Level')
    plt.ylabel('Median F1 Score')
    plt.legend(title='Tool')
    plt.ylim(0.85, 1)  # Restrict y-axis to start at 0.7
    plt.tight_layout()

    plt.savefig(sys.argv[2])

if __name__ == "__main__":
    main()

