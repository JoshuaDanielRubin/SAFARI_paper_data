import pandas as pd
import matplotlib.pyplot as plt

# Function to calculate specificity and sensitivity
def calculate_specificity_sensitivity(df):
    sensitivity = df['TP'] / (df['TP'] + df['FN'])  # Same as Recall
    specificity = df['TN'] / (df['TN'] + df['FP'])
    return specificity.mean(), sensitivity.mean()

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

    # Markers and colors for different aligners and damage types
    markers = {'giraffe': 'o', 'SAFARI': 's'}
    colors = ['blue', 'green', 'red', 'purple']

    # Plotting
    for i, damage_type in enumerate(damage_types):
        # Filter data by damage type
        giraffe_dt_data = giraffe_data[giraffe_data['Damage_Type'] == damage_type]
        safari_dt_data = safari_data[safari_data['Damage_Type'] == damage_type]

        # Calculate specificity and sensitivity
        giraffe_specificity, giraffe_sensitivity = calculate_specificity_sensitivity(giraffe_dt_data)
        safari_specificity, safari_sensitivity = calculate_specificity_sensitivity(safari_dt_data)

        # Print the values
        print(f"{damage_type} - giraffe: Specificity={giraffe_specificity:.12f}, Sensitivity={giraffe_sensitivity:.12f}")
        print(f"{damage_type} - SAFARI: Specificity={safari_specificity:.12f}, Sensitivity={safari_sensitivity:.12f}")

        # Plot points for giraffe and SAFARI data with different colors
        plt.scatter(giraffe_sensitivity, giraffe_specificity, color=colors[i], marker=markers['giraffe'], label=f'giraffe - {damage_type}')
        plt.scatter(safari_sensitivity, safari_specificity, color=colors[i], marker=markers['SAFARI'], label=f'SAFARI - {damage_type}')

    # Adding labels and title
    plt.xlabel('Sensitivity (True Positive Rate)')
    plt.ylabel('Specificity (True Negative Rate)')
    plt.title('Specificity vs Sensitivity: giraffe vs SAFARI (By Damage Type)', fontsize=18)
    plt.legend(loc='best')
    plt.grid(True)

    # Save the plot
    output_file_path = file_path.replace('.csv', '_specificity_sensitivity_plot.png')
    plt.savefig(output_file_path, dpi=300)
    plt.close()

    print(f"Plot saved as {output_file_path}")

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print("Usage: python this_script.py <path_to_alignment_stats.csv>")
    else:
        main(sys.argv[1])

