import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set global text properties
plt.rcParams.update({'font.size': 15,  # Adjust font size as needed
                     'font.weight': 'bold'})  # Make text bold

# Setting options for numpy and pandas
np.set_printoptions(precision=15)
pd.set_option('display.precision', 15)

def generate_percentage_decrease_plots(input_file_path):
    # Load the dataset
    data = pd.read_csv(input_file_path)

    # Reordering the damage levels
    data['damage_level'] = pd.Categorical(data['damage_level'],
                                          categories=["None", "Single-stranded", "Mid", "High"],
                                          ordered=True)

    # Calculate the median RMSE across samples by tool, for each 'fragment_len_dist' and 'damage_level'
    median_rmse = data.groupby(['fragment_len_dist', 'damage_level', 'tool'])['RMSE'].median().reset_index()

    # Pivot the table to have tools as columns
    pivot_data = median_rmse.pivot_table(index=['fragment_len_dist', 'damage_level'], columns='tool', values='RMSE').reset_index()

    # Melt the pivot_data for plotting
    melted_data = pivot_data.melt(id_vars=['fragment_len_dist', 'damage_level'], var_name='Tool', value_name='RMSE')

    # Order and palette customization
    tools_order = ['SAFARI', 'vg giraffe']  # Add other tools to the list as needed
    tools_palette = {'SAFARI': 'green', 'vg giraffe': 'orange'}  # Add distinct colors for other tools
    
    # You can dynamically adjust the order and palette if the tools vary
    # For example, add any tools not already included to the order list and generate colors
    additional_tools = sorted(set(melted_data['Tool']) - set(tools_order))
    tools_order += additional_tools
    # Use a colormap or define a list of colors for additional tools
    cmap = plt.get_cmap("tab10")  # This is an example; adjust as needed
    for i, tool in enumerate(additional_tools):
        tools_palette[tool] = cmap(i % cmap.N)

    # Creating the plots
    fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
    
    # Filtering data for each site and plotting with customized order and palette
    for i, site in enumerate(["Chagyrskaya", "Vindija"]):
        site_data = melted_data[melted_data['fragment_len_dist'] == site]
        sns.barplot(x='damage_level', y='RMSE', hue='Tool', data=site_data, ax=axs[i], ci=None, hue_order=tools_order, palette=tools_palette)
        axs[i].set_title(site)
        axs[i].set_ylabel('Median RMSE')
        axs[i].set_xlabel('Damage Level')

    plt.tight_layout()
    plt.savefig("RMSE_all.png")

if __name__ == '__main__':
    # Directly running the script requires the input file path as an argument
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script_name.py input_file_path")
    else:
        input_file_path = sys.argv[1]
        generate_percentage_decrease_plots(input_file_path)

