import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def load_data(file_path):
    return pd.read_csv(file_path)

def calculate_percent_change(data):
    giraffe_scores = data[data['Tool'] == 'giraffe'].set_index('Damage')
    
    def percent_change(row, column):
        giraffe_score = giraffe_scores.loc[row['Damage'], column]
        return ((row[column] - giraffe_score) / giraffe_score) * 100
    
    data['% Change All Reads'] = data.apply(percent_change, args=('F1 Score (All Reads)',), axis=1)
    data['% Change Mito Reads'] = data.apply(percent_change, args=('F1 Score (Mitochondrial Reads)',), axis=1)
    return data

def plot_data_custom_bold_large_corrected(data):
    data_filtered = data[data['Tool'] != 'giraffe']
    color_palette = ['#6495ED', '#FF69B4', '#BA55D3', '#20B2AA', '#87CEFA', '#32CD32', '#FFD700']
    custom_order = {'safari': 'SAFARI', 'aln': 'BWA-ALN', 'aln_anc': 'BWA-ALN (anc)', 'mem': 'BWA-MEM', 'shrimp': 'SHRiMP', 'bowtie2': 'Bowtie2', 'bb': 'BBMAP'}
    damages_custom_labels = {'dnone': 'None', 'dsingle': 'Single', 'ddmid': 'Mid', 'ddhigh': 'High'}
    damages = sorted(data['Damage'].unique(), key=lambda x: ['dnone', 'dsingle', 'ddmid', 'ddhigh'].index(x))
    
    data_filtered['Tool'] = data_filtered['Tool'].map(custom_order)
    tools_order = ['SAFARI', 'BWA-ALN', 'BWA-ALN (anc)', 'BWA-MEM', 'SHRiMP', 'Bowtie2', 'BBMAP']
    x = np.arange(len(tools_order))
    width = 0.2
    
    fig, axs = plt.subplots(2, 1, figsize=(20, 18))
    
    for i, dmg in enumerate(damages):
        mito_scores = data_filtered[data_filtered['Damage'] == dmg].groupby('Tool')['% Change Mito Reads'].mean().reindex(tools_order)
        axs[0].bar(x - width*1.5 + i*width, mito_scores, width, label=damages_custom_labels[dmg], color=color_palette[i % len(color_palette)])
    
    for i, dmg in enumerate(damages):
        all_scores = data_filtered[data_filtered['Damage'] == dmg].groupby('Tool')['% Change All Reads'].mean().reindex(tools_order)
        axs[1].bar(x - width*1.5 + i*width, all_scores, width, label=damages_custom_labels[dmg], color=color_palette[i % len(color_palette)])
    
    axs[0].set_xticks([])  # Removing x-ticks entirely for the first subplot
    
    title_fontsize = 24
    label_fontsize = 22
    tick_labelsize = 20
    legend_fontsize = 18
    
    axs[0].set_ylabel('Percent Change in F1 Score', fontsize=label_fontsize, fontweight='bold')
    axs[0].set_title('Percent Change in F1 Score for Mitochondrial Reads from vg giraffe', fontsize=title_fontsize, fontweight='bold')
    
    axs[1].set_ylabel('Percent Change in F1 Score', fontsize=label_fontsize, fontweight='bold')
    axs[1].set_title('Percent Change in F1 Score for All Reads from vg giraffe', fontsize=title_fontsize, fontweight='bold')
    axs[1].tick_params(axis='x', labelsize=tick_labelsize, labelrotation=45)
    axs[1].tick_params(axis='y', labelsize=tick_labelsize)
    
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(tools_order, rotation=45, fontsize=tick_labelsize, fontweight='bold')
    
    axs[0].legend(title='Damage', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=legend_fontsize, title_fontsize='large')
    axs[1].legend(title='Damage', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=legend_fontsize, title_fontsize='large')
    
    plt.tight_layout()
    plt.savefig('linear.png')  # Saving the plot as 'linear.png'
    plt.show()

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_data_file>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    data = load_data(file_path)
    data = calculate_percent_change(data)
    plot_data_custom_bold_large_corrected(data)

if __name__ == "__main__":
    main()

