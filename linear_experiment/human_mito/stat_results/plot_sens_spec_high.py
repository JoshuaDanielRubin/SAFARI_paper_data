
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Load the dataset directly from the provided file path
file_path = sys.argv[1]  # Update this path as necessary
data = pd.read_csv(file_path)

# Filter data for samples where the tool is either SAFARI or vg giraffe and damage level is 'High'
filtered_data = data[(data['tool'].isin(['SAFARI', 'vg giraffe'])) & (data['damage_level'] == 'High')]

# Compute median sensitivity for each (k,w) pair and tool
medians = filtered_data.groupby(['k', 'w', 'tool'])['sensitivity'].median().reset_index()

# Pivot the data to get SAFARI and vg giraffe side by side for each (k,w) pair
pivot_medians = medians.pivot(index=['k', 'w'], columns='tool', values='sensitivity').reset_index()

# Convert k and w from float to int for cleaner axis labels
pivot_medians['k'] = pivot_medians['k'].astype(int)
pivot_medians['w'] = pivot_medians['w'].astype(int)

# Calculate the difference in median sensitivity (SAFARI - vg giraffe)
pivot_medians['sensitivity_diff'] = pivot_medians['SAFARI'] - pivot_medians['vg giraffe']
pivot_medians_sorted_by_diff = pivot_medians.sort_values(by='sensitivity_diff', ascending=False)

# Plotting
plt.figure(figsize=(14, 10))
formatted_labels = pivot_medians_sorted_by_diff['k'].astype(str) + ',' + pivot_medians_sorted_by_diff['w'].astype(str)
bars = plt.bar(formatted_labels, pivot_medians_sorted_by_diff['sensitivity_diff'], color=np.where(pivot_medians_sorted_by_diff['sensitivity_diff'] > 0, 'green', 'orange'))

plt.axhline(0, color='grey', linewidth=0.8)
plt.xlabel('(k,w) Values', fontsize=14, fontweight='bold')
plt.ylabel('Difference in Median Sensitivity\n(SAFARI - vg giraffe)', fontsize=14, fontweight='bold')
plt.title('Difference in Median Sensitivity Between SAFARI and vg giraffe\non High-damage Samples ('+ sys.argv[3] + ")", fontsize=16, fontweight='bold', pad=20)
plt.xticks(ticks=np.arange(len(formatted_labels)), labels=formatted_labels, rotation=45, ha="right", fontsize=12, fontweight='bold')

# Adjusting legend size
legend_elements = [plt.Line2D([0], [0], color='green', lw=4, label='SAFARI'),
                   plt.Line2D([0], [0], color='orange', lw=4, label='vg giraffe')]
plt.legend(handles=legend_elements, title='Tool', title_fontsize='20', fontsize=18, loc='upper right', frameon=True, shadow=True)

plt.tight_layout()
plt.savefig(sys.argv[2])

