import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Adjust font size and remove bold weight for the global text properties
plt.rcParams.update({'font.size': 15, 'font.weight': 'normal'})

# Configure numpy and pandas display precision
np.set_printoptions(precision=15)
pd.set_option('display.precision', 15)

def plot_rmse_by_damage_level(input_file_path):
    # Load the data from the specified input file
    dataset = pd.read_csv(input_file_path)

    # Ensure the damage levels are ordered according to your specifications
    dataset['damage_level'] = pd.Categorical(dataset['damage_level'],
                                             categories=["None", "Single-stranded", "Mid", "High"],
                                             ordered=True)

    # Group the data by fragment length distribution, damage level, and tool, then calculate median RMSE
    median_rmse_data = dataset.groupby(['fragment_len_dist', 'damage_level', 'tool'])['RMSE'].median().reset_index()

    # Convert the grouped data into a pivot table for easier plotting
    pivot_rmse = median_rmse_data.pivot_table(index=['fragment_len_dist', 'damage_level'], columns='tool', values='RMSE').reset_index()

    # Melt the pivot table to long-form for seaborn plotting
    long_form_data = pivot_rmse.melt(id_vars=['fragment_len_dist', 'damage_level'], var_name='Tool', value_name='RMSE')

    print(long_form_data)

    # Define the order of tools and their corresponding colors
    tool_order = ['SAFARI', 'vg giraffe']  # Specify the preferred order for these tools
    color_map = {'SAFARI': 'green', 'vg giraffe': 'orange'}  # Assign colors to the specified tools

    # Append additional tools to the order list and assign colors
    additional_tools = sorted(set(long_form_data['Tool']) - set(tool_order))
    tool_order += additional_tools
    # Define a color palette for any additional tools
    additional_colors = sns.color_palette("husl", len(additional_tools))
    for tool, color in zip(additional_tools, additional_colors):
        color_map[tool] = color

    # Create the plots for each site
    fig, axes = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
    for i, site in enumerate(["Chagyrskaya", "Vindija"]):
        site_specific_data = long_form_data[long_form_data['fragment_len_dist'] == site]
        sns.barplot(x='damage_level', y='RMSE', hue='Tool', data=site_specific_data, ax=axes[i], ci=None, hue_order=tool_order, palette=color_map)
        axes[i].set_title(site)
        axes[i].set_ylabel('Median RMSE')

    axes[-1].set_xlabel('Damage Level')

    plt.tight_layout()
    plt.savefig("RMSE_all_tools_comparison.png")

    # Modify the legend to not use bold font
    for ax in axes:
        legend = ax.get_legend()
        for text in legend.get_texts():
            text.set_fontweight('normal')

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python <script_name>.py <input_file_path>")
    else:
        file_path = sys.argv[1]
        plot_rmse_by_damage_level(file_path)

