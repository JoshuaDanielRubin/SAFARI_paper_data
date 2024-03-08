def load_and_prepare_data_adjusted(linear_stats_path, k10_w2_path):
    linear_stats_df = pd.read_csv(linear_stats_path)
    k10_w2_df = pd.read_csv(k10_w2_path)
    combined_df = pd.concat([linear_stats_df, k10_w2_df])
    
    # Adjusting tool names
    tool_name_adjustments = {
        'aln': 'BWA ALN',
        'aln (anc)': 'BWA ALN (anc)',
        'bbmap': 'BBMap',
        'bowtie2': 'Bowtie2',
        'bwa-mem': 'BWA-MEM',
        'shrimp': 'SHRiMP',
        'giraffe': 'vg giraffe',
        'safari': 'SAFARI'
    }
    combined_df['Tool_Display'] = combined_df['Tool'].apply(lambda x: tool_name_adjustments.get(x.lower(), x))
    
    # Correcting the damage level adjustments
    damage_level_adjustments_corrected = {
        'ddlow': 'Single',
        'ddmid': 'Mid',
        'ddhigh': 'High',
        'none': 'None',
        'dnone': 'None',
        'dsingle': 'Single'
    }
    combined_df['Damage_Level'] = combined_df['Damage'].apply(lambda x: damage_level_adjustments_corrected.get(x, x))
    
    return combined_df

def plot_data_final(combined_df):
    metrics_of_interest_corrected = ['Mito reads mapped', 'Bacteria reads mapped', 'NuMT reads mapped']
    subplot_titles = [
        'Mitochondrial Reads Correctly Mapped',
        'Bacterial Reads Spuriously Mapped',
        'NuMT Reads Spuriously Mapped'
    ]
    damage_order = ['None', 'Single', 'Mid', 'High']
    
    fig, axes = plt.subplots(3, 1, figsize=(14, 22))
    plt.rc('font', size=12, weight='bold', style='italic')
    
    for i, metric in enumerate(metrics_of_interest_corrected):
        sns.barplot(data=combined_df[combined_df['Metric'] == metric], x='Tool_Display', y='Value', hue='Damage_Level',
                    hue_order=damage_order, ax=axes[i])
        axes[i].set_title(subplot_titles[i], fontsize=16, fontweight='bold')
        axes[i].set_ylabel('Count', fontsize=14, fontweight='bold')
        axes[i].set_xlabel('Tool', fontsize=14, fontweight='bold')
        axes[i].tick_params(axis='x', labelrotation=45)
        axes[i].tick_params(axis='both', labelsize=12)

        # Adjust y-axis limit for the first plot to avoid legend overlap and disable scientific notation
        if i == 0:
            ymin, ymax = combined_df[combined_df['Metric'] == metric]['Value'].min(), combined_df[combined_df['Metric'] == metric]['Value'].max()
            delta = (ymax - ymin) * 0.1
            ymax_adjusted = ymax + delta + (ymax * 0.05)  # Increase max limit a bit more
            axes[i].set_ylim(ymin - delta, ymax_adjusted)
            axes[i].legend(title='Damage Level', title_fontsize='13', fontsize='12', loc='upper left')
            axes[i].ticklabel_format(style='plain', axis='y')  # Disable scientific notation
        else:
            ymin, ymax = combined_df[combined_df['Metric'] == metric]['Value'].min(), combined_df[combined_df['Metric'] == metric]['Value'].max()
            delta = (ymax - ymin) * 0.1
            axes[i].set_ylim(ymin - delta, ymax + delta)
            axes[i].legend(title='Damage Level', title_fontsize='13', fontsize='12', loc='upper left')

    plt.tight_layout()
    return plt

# Re-load and prepare data with adjusted functions
combined_df_final = load_and_prepare_data_adjusted(linear_stats_path, k10_w2_path)

# Generate and display the final adjusted plot
plot_final = plot_data_final(combined_df_final)
plot_final.show()

