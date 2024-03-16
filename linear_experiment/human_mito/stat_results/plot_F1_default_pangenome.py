import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():
    # Load the data
    file_path = sys.argv[1]
    data = pd.read_csv(file_path)

    # Subset the dataframe
    subset_df = data[(data['k'] == 29) & (data['w'] == 11) & (data['tool'].str.strip().str.lower().isin(['safari', 'vg giraffe']))]

    # Calculate median F1 scores
    median_f1_scores = subset_df.groupby(['tool', 'damage_level'])['f1'].median().reset_index()

    # Define custom order for damage levels
    damage_order = ['None', 'Single-stranded', 'Mid', 'High']
    median_f1_scores['damage_level'] = pd.Categorical(median_f1_scores['damage_level'], categories=damage_order, ordered=True)

    # Sort the DataFrame by the 'damage_level' to ensure the plot follows the custom order
    median_f1_scores.sort_values('damage_level', inplace=True)
    print(median_f1_scores)

    # Correcting the case sensitivity issue for the palette
    corrected_palette = {'vg giraffe': 'orange', 'SAFARI': 'green'}  # Ensure correct case for consistency

    # Create a barplot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=median_f1_scores, x='damage_level', y='f1', hue='tool', palette=corrected_palette)
    plt.title('Median F1 Score by Pangenome Tool \n Stratified by Damage Matrix (Default Parameters k=29, w=11, ' + sys.argv[3] + ")")
    plt.xlabel('Damage Level')
    plt.ylabel('Median F1 Score')
    plt.legend(title='Tool')
    plt.ylim(0.6, 1)  # Restrict y-axis to start at 0.7
    plt.tight_layout()

    plt.savefig(sys.argv[2])

if __name__ == "__main__":
    main()

