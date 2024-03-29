import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import random

def calculate_sensitivity(data):
    # No additional calculation needed for sensitivity, assuming it's already a column in the data
    return data

def plot_best_sensitivity(data_path):
    # Load data
    data = pd.read_csv(data_path)

    # Filtering data for 'None' damage level and calculating sensitivity
    filtered_data = calculate_sensitivity(data[(data['damage_level'] == 'None')].copy())
    # Get medians of sensitivity grouped by tool, k, and w
    medians_sensitivity = filtered_data.groupby(['tool', 'k', 'w'])['sensitivity'].median().reset_index()
    # Identify the best (k,w) combination for each tool based on median sensitivity
    best_sensitivity = medians_sensitivity.loc[medians_sensitivity.groupby('tool')['sensitivity'].idxmax()]
    best_sensitivity_data = filtered_data.merge(best_sensitivity[['tool', 'k', 'w']], on=['tool', 'k', 'w'])

    # Plotting
    plt.figure(figsize=(10, 8))
    
    # Predefined colors for specific tools
    predefined_colors = {'SAFARI': 'green', 'vg giraffe': 'orange'}
    
    # Get a list of all tools to ensure unique colors for those not predefined
    tools = best_sensitivity_data['tool'].unique()
    color_palette = plt.cm.get_cmap('hsv', len(tools)) # hsv colormap, adjust as needed
    
    # Assign colors ensuring unique colors for each tool
    tool_colors = {}
    for i, tool in enumerate(tools):
        if tool in predefined_colors:
            tool_colors[tool] = predefined_colors[tool]
        else:
            while color_palette(i) in predefined_colors.values(): # Ensure unique color assignment
                i += 1
            tool_colors[tool] = color_palette(i)
    
    # Modify the plotting loop
    for tool in tools:
        tool_data = best_sensitivity_data[best_sensitivity_data['tool'] == tool]
        print(tool_data)
        if tool == 'vg giraffe':
            # Adding slight jitter to 'vg giraffe' data points
            jittered_sensitivity = tool_data['sensitivity']
            jittered_specificity = tool_data['specificity']
            plt.scatter(jittered_sensitivity, jittered_specificity, label=f'{tool}', color=tool_colors[tool], alpha=0.7, edgecolor='black', linewidth=0.5)
        else:
            plt.scatter(tool_data['sensitivity'], tool_data['specificity'], label=f'{tool}', color=tool_colors[tool], alpha=0.7, marker='^')

    plt.title('Sensitivity vs Specificity for Best (k,w) \n Parameters Based on Median Sensitivity (' + sys.argv[3] + ", No Damage)", fontsize=16, fontweight='bold')
    plt.xlabel('Sensitivity', fontsize=14, fontweight='bold')
    plt.ylabel('Specificity', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(sys.argv[2])
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <data_path> <output path> <Dist>")
        sys.exit(1)
    data_path = sys.argv[1]
    plot_best_sensitivity(data_path)

