import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_sensitivity_k29_w11(data_path, output_path, dist):
    # Load data
    data = pd.read_csv(data_path)

    # Filtering data for 'None' damage level and specific k=29 and w=11
    filtered_data = data[(data['damage_level'] == 'None') & (data['k'] == 29) & (data['w'] == 11)]

    # Plotting
    plt.figure(figsize=(10, 8))
    predefined_colors = {'SAFARI': 'green', 'vg giraffe': 'orange'}
    tools = filtered_data['tool'].unique()
    color_palette = plt.cm.get_cmap('hsv', len(tools))

    tool_colors = {}
    for i, tool in enumerate(tools):
        tool_colors[tool] = predefined_colors.get(tool, color_palette(i))

    for tool in tools:
        tool_data = filtered_data[filtered_data['tool'] == tool]
        if tool == 'vg giraffe':
            plt.scatter(tool_data['sensitivity'], tool_data['specificity'], label=tool, color=tool_colors[tool], alpha=1.0, edgecolor='black', linewidth=0.5, s=80)
        else:
            plt.scatter(tool_data['sensitivity'], tool_data['specificity'], label=tool, color=tool_colors[tool], alpha=0.3, marker='^', s=80)

    plt.title(f'Sensitivity vs Specificity (k=29, w=11, No Damage, {dist})', fontsize=16, fontweight='bold')
    plt.xlabel('Sensitivity', fontsize=14, fontweight='bold')
    plt.ylabel('Specificity', fontsize=14, fontweight='bold')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <data_path> <output_path> <Dist>")
        sys.exit(1)
    data_path = sys.argv[1]
    output_path = sys.argv[2]
    dist = sys.argv[3]
    plot_sensitivity_k29_w11(data_path, output_path, dist)

