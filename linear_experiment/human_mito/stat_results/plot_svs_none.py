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

    # Placeholder for the rest of your existing code before plotting...
    filtered_data = calculate_f1(data[(data['damage_level'] == 'None')].copy())
    medians_f1 = filtered_data.groupby(['tool', 'k', 'w'])['F1'].median().reset_index()
    best_f1 = medians_f1.loc[medians_f1.groupby('tool')['F1'].idxmax()]
    best_f1_data = filtered_data.merge(best_f1[['tool', 'k', 'w']], on=['tool', 'k', 'w'])

    # Plotting
    plt.figure(figsize=(10, 8))
    
    # Predefined colors for specific tools
    predefined_colors = {'SAFARI': 'green', 'vg giraffe': 'orange'}
    
    # Get a list of all tools to ensure unique colors for those not predefined
    tools = best_f1_data['tool'].unique()
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
    
    # Plot each tool with its assigned color
    for tool in tools:
        tool_data = best_f1_data[best_f1_data['tool'] == tool]
        plt.scatter(tool_data['specificity'], tool_data['sensitivity'], label=f'{tool}', color=tool_colors[tool], alpha=0.7)

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

