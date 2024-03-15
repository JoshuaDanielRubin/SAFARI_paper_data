import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_differences(data, column_name, y_label, title, file_name, legend_loc):
    medians = data.groupby(['k', 'w', 'tool'])[column_name].median().reset_index()
    pivot_medians = medians.pivot(index=['k', 'w'], columns='tool', values=column_name).reset_index()
    pivot_medians['k'] = pivot_medians['k'].astype(int)
    pivot_medians['w'] = pivot_medians['w'].astype(int)
    pivot_medians['diff'] = pivot_medians['SAFARI'] - pivot_medians['vg giraffe']
    pivot_medians_sorted_by_diff = pivot_medians.sort_values(by='diff', ascending=False)
    plt.figure(figsize=(14, 10))
    formatted_labels = pivot_medians_sorted_by_diff['k'].astype(str) + ',' + pivot_medians_sorted_by_diff['w'].astype(str)
    plt.bar(formatted_labels, pivot_medians_sorted_by_diff['diff'], color=np.where(pivot_medians_sorted_by_diff['diff'] > 0, 'green', 'orange'))
    plt.axhline(0, color='grey', linewidth=0.8)
    plt.xlabel('(k,w) Values', fontsize=14, fontweight='bold')
    plt.ylabel(y_label, fontsize=14, fontweight='bold')
    plt.title(title, fontsize=16, fontweight='bold', pad=20)
    plt.xticks(ticks=np.arange(len(formatted_labels)), labels=formatted_labels, rotation=45, ha="right", fontsize=12, fontweight='bold')
    plt.ylim([pivot_medians_sorted_by_diff['diff'].min() * 1.1, pivot_medians_sorted_by_diff['diff'].max() * 1.1])
    legend_elements = [plt.Line2D([0], [0], color='green', lw=4, label='SAFARI'), plt.Line2D([0], [0], color='orange', lw=4, label='vg giraffe')]
    plt.legend(handles=legend_elements, title='Tool', title_fontsize='20', fontsize=18, loc=legend_loc, frameon=True, shadow=True)
    plt.tight_layout()
    plt_path = file_name
    plt.savefig(plt_path)
    plt.close()

# Load the dataset
data = pd.read_csv(sys.argv[1])

# Filter data for samples where the tool is either SAFARI or vg giraffe and damage level is 'High'
filtered_data = data[(data['tool'].isin(['SAFARI', 'vg giraffe'])) & (data['damage_level'] == 'High')]

# Generate both plots
plot_differences(filtered_data, 'sensitivity', 'Difference in Median Sensitivity\n(SAFARI - vg giraffe)', 'Difference in Median Sensitivity Between SAFARI and vg giraffe\non High-damage Samples', sys.argv[2], 'upper right')
plot_differences(filtered_data, 'specificity', 'Difference in Median Specificity\n(SAFARI - vg giraffe)', 'Difference in Median Specificity Between SAFARI and vg giraffe\non High-damage Samples', sys.argv[3], 'upper left')
