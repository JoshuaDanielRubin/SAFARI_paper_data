import pandas as pd
import matplotlib.pyplot as plt

def calculate_precision_recall(df):
    precision = df['TP'] / (df['TP'] + df['FP'])
    recall = df['TP'] / (df['TP'] + df['FN'])
    return precision.mean(), recall.mean()

def main(file_path):
    # Load data
    data = pd.read_csv(file_path)

    # Filter data for 'giraffe' and 'SAFARI'
    giraffe_data = data[data['Aligner_Name'].str.lower() == 'giraffe']
    safari_data = data[data['Aligner_Name'].str.upper() == 'SAFARI']

    # Damage types
    damage_types = data['Damage_Type'].unique()

    # Prepare plot
    plt.figure(figsize=(10, 8))

    # Markers for different aligners
    markers = {'giraffe': 'o', 'SAFARI': 's'}

    # Plotting
    for damage_type in damage_types:
        # Filter data by damage type
        giraffe_dt_data = giraffe_data[giraffe_data['Damage_Type'] == damage_type]
        safari_dt_data = safari_data[safari_data['Damage_Type'] == damage_type]

        # Calculate precision and recall
        giraffe_precision, giraffe_recall = calculate_precision_recall(giraffe_dt_data)
        safari_precision, safari_recall = calculate_precision_recall(safari_dt_data)

        # Plot points for giraffe and SAFARI data
        plt.scatter(giraffe_recall, giraffe_precision, marker=markers['giraffe'], label=f'giraffe - {damage_type}')
        plt.scatter(safari_recall, safari_precision, marker=markers['SAFARI'], label=f'SAFARI - {damage_type}')

    # Adding labels and title
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision vs Recall: giraffe vs SAFARI (By Damage Type)', fontsize=18)
    plt.legend(loc='best')
    plt.grid(True)

    # Save the plot with high resolution
    plt.savefig('precision_recall_comparison.png', dpi=300)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print("Usage: python this_script.py <path_to_alignment_stats.csv>")
    else:
        main(sys.argv[1])

