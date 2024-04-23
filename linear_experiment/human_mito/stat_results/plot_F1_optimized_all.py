import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np

def optimize_k_w(group):
    optimal_row = group.loc[group['sensitivity'].idxmax()]
    print(f"Selected k: {optimal_row['k']}, w: {optimal_row['w']}")
    return optimal_row[['k', 'w']]

def main():
    file_path = sys.argv[1]
    data = pd.read_csv(file_path)

    optimal_params = data.groupby('tool').apply(optimize_k_w).reset_index()
    optimal_data = pd.merge(data, optimal_params, on=['tool', 'k', 'w'])
    median_f1_scores = optimal_data.groupby(['tool', 'damage_level'])['f1'].median().reset_index()

    damage_order = ['None', 'Single-stranded', 'Mid', 'High']
    median_f1_scores['damage_level'] = pd.Categorical(median_f1_scores['damage_level'], categories=damage_order, ordered=True)

    # Determine the order of tools based on their performance at 'High' damage level
    high_damage_order = median_f1_scores[median_f1_scores['damage_level'] == 'High'].sort_values('f1', ascending=False)['tool'].tolist()

    corrected_palette = {
        'vg giraffe': 'orange', 
        'SAFARI': 'green', 
        'BBMAP': '#1f7872',
        'SHRiMP': '#6996b3',
        'BWA-MEM': '#5e3a8c',
        'BWA ALN': '#ff99cc',
        'BWA ALN (anc)': '#36454f',
        'Bowtie2': '#a45a52'
    }

    plt.figure(figsize=(10, 6))
    sns.barplot(data=median_f1_scores, x='damage_level', y='f1', hue='tool', palette=corrected_palette, hue_order=high_damage_order)
    plt.title(r'Median $F_1$ Score (Optimized Parameters, ' + sys.argv[3] + ")")
    plt.xlabel('Damage Level')
    plt.ylabel(r'Median $F_1$ Score')
    plt.legend(title='Tool', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.ylim(0.90, 1)
    plt.tight_layout()

    plt.savefig(sys.argv[2])

if __name__ == "__main__":
    main()

