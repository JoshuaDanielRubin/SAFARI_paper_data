import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():
    # Load the data directly from the input file
    data = pd.read_csv(sys.argv[1])

    # Calculate median F1 scores for all tools and damage levels
    median_f1_scores = data.groupby(['tool', 'damage_level'])['f1'].median().reset_index()

    # Define custom order for damage levels and tools
    damage_order = ['None', 'Single-stranded', 'Mid', 'High']
    tools_order = ['SAFARI', 'vg giraffe'] + [tool for tool in median_f1_scores['tool'].unique() if tool not in ['SAFARI', 'vg giraffe']]

    median_f1_scores['damage_level'] = pd.Categorical(median_f1_scores['damage_level'], categories=damage_order, ordered=True)
    median_f1_scores['tool'] = pd.Categorical(median_f1_scores['tool'], categories=tools_order, ordered=True)

    # Find the minimum F1 score to adjust y-axis
    min_f1_score = median_f1_scores['f1'].min()
    lower_limit = round(min_f1_score - 0.05, 2)

    # Define custom colors
    remaining_tools = [tool for tool in tools_order if tool not in ['SAFARI', 'vg giraffe']]
    other_colors = sns.color_palette('bright', len(remaining_tools))
    custom_palette = {'SAFARI': 'green', 'vg giraffe': 'orange'}
    for tool, color in zip(remaining_tools, other_colors):
        custom_palette[tool] = color

    # Create the plot with custom settings
    plt.figure(figsize=(16, 9))
    sns.barplot(data=median_f1_scores, x='damage_level', y='f1', hue='tool', palette=custom_palette)
    plt.title('Median F1 Score by Tool Stratified by Damage Level', fontsize=16)
    plt.xlabel('Damage Level', fontsize=14)
    plt.ylabel('Median F1 Score', fontsize=14)
    plt.legend(title='Tool', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    plt.ylim(lower_limit, 1)
    plt.tight_layout()

    # Save the plot
    plt.savefig(sys.argv[2])

if __name__ == '__main__':
    main()

