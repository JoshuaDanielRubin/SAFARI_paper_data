import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():
    # Load the data
    file_path = sys.argv[1]
    data = pd.read_csv(file_path)

    # Subset the dataframe
    subset_df = data[(data['k'] == 29) & (data['w'] == 11) & (data['tool'].str.lower().isin(['SAFARI', 'vg giraffe']))]

    # Calculate median F1 scores
    median_f1_scores = subset_df.groupby(['tool', 'damage_level'])['f1'].median().reset_index()

    # Correcting the case sensitivity issue for the palette
    corrected_palette = {'vg giraffe': 'orange', 'SAFARI': 'green'}

    # Create a barplot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=median_f1_scores, x='damage_level', y='f1', hue='tool', palette=corrected_palette)
    plt.title('Median F1 Score by Pangenome Tool \n Stratified by Damage Matrix (Default Parameters k=29, w=11)')
    plt.xlabel('Damage Level')
    plt.ylabel('Median F1 Score')
    plt.legend(title='Tool')
    plt.ylim(0.7, 1)  # Restrict y-axis to start at 0.7
    plt.tight_layout()

    plt.savefig(sys.argv[2])

if __name__ == "__main__":
    main()

