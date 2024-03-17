
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def generate_percentage_decrease_plots(input_file_path):
    # Load the dataset
    data = pd.read_csv(input_file_path)

    # Reordering the damage levels
    data['damage_level'] = pd.Categorical(data['damage_level'], 
                                          categories=["None", "Single-stranded", "Mid", "High"],
                                          ordered=True)

    # Filter the dataset based on the criteria
    filtered_data = data[~((data['tool'].isin(['SAFARI', 'vg giraffe'])) & ~((data['k'] == '29') & (data['w'] == '11')))]

    # Calculate the median RMSE across samples by tool, for each 'fragment_len_dist' and 'damage_level'
    median_rmse = filtered_data.groupby(['fragment_len_dist', 'damage_level', 'tool'])['RMSE'].median().reset_index()

    # Pivot the table to have tools as columns
    pivot_data = median_rmse.pivot_table(index=['fragment_len_dist', 'damage_level'], columns='tool', values='RMSE').reset_index()

    # Calculate the percentage decrease from vg giraffe to SAFARI
    pivot_data['percentage_decrease'] = ((pivot_data['vg giraffe'] - pivot_data['SAFARI']) / pivot_data['vg giraffe']) * 100

    # Melt the pivot_data for plotting
    melted_data = pivot_data.melt(id_vars=['fragment_len_dist', 'damage_level'], value_vars=['percentage_decrease'])

    print(melted_data)

    # Set the aesthetic style of the plots, including making the text larger and bolder
    sns.set_style("whitegrid")
    sns.set_context("talk", font_scale=1)

    # Create the bar plot
    plt.figure(figsize=(12, 8))
    sns.barplot(x='damage_level', y='value', hue='fragment_len_dist', data=melted_data, palette="coolwarm")
    plt.title('Percentage Decrease in Median RMSE from vg giraffe to SAFARI\nacross Different Damage Levels', fontsize=16, fontweight='bold', loc='left')
    plt.xlabel('Damage Level', fontsize=14, fontweight='bold')
    plt.ylabel('Percentage Decrease in Median RMSE', fontsize=14, fontweight='bold')
    plt.ylim(-0.01, 0.065)
    plt.legend(title='Fragment Length Dist.', fontsize='medium', title_fontsize='14')

    plt.tight_layout()
    plt.savefig(sys.argv[2])
    plt.close()

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print("Usage: python generate_percentage_decrease_plots.py <input_file_path>")
    else:
        input_file_path = sys.argv[1]
        generate_percentage_decrease_plots(input_file_path)

