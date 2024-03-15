import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def calculate_f1(data):
    # Placeholder for F1 calculation; replace with appropriate calculation method
    data['F1'] = 2 * (data['sensitivity'] * data['specificity']) / (data['sensitivity'] + data['specificity'])
    return data

def plot_best_f1(data_path):
    # Load data
    data = pd.read_csv(data_path)
    
    # Assuming the input data already contains 'sensitivity', 'specificity', 'k', 'w', and 'tool' columns
    # Filter for 'High' damage level samples and calculate F1 scores
    filtered_data = data[(data['damage_level'] == 'High')].copy()
    filtered_data = calculate_f1(filtered_data)

    # Group by tool, k, and w, then calculate median F1 score
    medians_f1 = filtered_data.groupby(['tool', 'k', 'w'])['F1'].median().reset_index()

    # Find the best (k,w) based on median F1 for each tool
    best_f1 = medians_f1.loc[medians_f1.groupby('tool')['F1'].idxmax()]

    # Filter data for the best (k,w) values for each tool
    best_f1_data = filtered_data.merge(best_f1[['tool', 'k', 'w']], on=['tool', 'k', 'w'])

    # Plotting
    plt.figure(figsize=(10, 8))
    for tool in best_f1_data['tool'].unique():
        tool_data = best_f1_data[best_f1_data['tool'] == tool]
        plt.scatter(tool_data['specificity'], tool_data['sensitivity'], label=f'{tool}', alpha=0.7)

    plt.title('Sensitivity vs Specificity for Best (k,w) Parameters Based on F1 Score', fontsize=16, fontweight='bold')
    plt.xlabel('Specificity', fontsize=14, fontweight='bold')
    plt.ylabel('Sensitivity', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(sys.argv[2])
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py <data_path> <output path>")
        sys.exit(1)
    data_path = sys.argv[1]
    plot_best_f1(data_path)

