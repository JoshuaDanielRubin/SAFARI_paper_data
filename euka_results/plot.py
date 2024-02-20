import pandas as pd
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap

def read_and_parse_files(directory_path):
    """
    Reads and parses files in the specified directory that match the pattern.
    
    Parameters:
    - directory_path: str, path to the directory containing the files.
    
    Returns:
    - DataFrame with all the parsed data.
    """
    all_data = []

    # Loop through each file in the directory
    for file_name in os.listdir(directory_path):
        if file_name.endswith(".tsv") and "_corrected" in file_name and "detected" in file_name:
            threshold = float(file_name.split("_")[-2])
            file_path = os.path.join(directory_path, file_name)
            
            # Read the file, skipping the header line that starts with #
            with open(file_path, 'r') as file:
                lines = file.readlines()
                data = [line.strip().split('\t') for line in lines if not line.startswith('#')]
                for row in data:
                    taxa, detected, number_of_reads = row[0], row[1], int(row[2])
                    all_data.append({"Taxa": taxa, "Threshold": threshold, "Number_of_reads": number_of_reads})

    # Convert the list of dictionaries to a DataFrame
    return pd.DataFrame(all_data)

def plot_data(df):
    # Calculate total reads per taxa
    df_total_reads = df.groupby('Taxa').sum().reset_index()

    # Sort taxa by total read count in descending order
    sorted_taxa = df_total_reads.sort_values(by='Number_of_reads', ascending=False)['Taxa']

    # Pivot the original dataframe for plotting
    df_pivot = df.pivot(index="Taxa", columns="Threshold", values="Number_of_reads").fillna(0)

    # Reorder the DataFrame based on the sorted taxa
    df_pivot = df_pivot.reindex(sorted_taxa)
    print(df_pivot)

    # Proceed with plotting as before
    df_pivot.plot(kind='bar', stacked=True, figsize=(14, 8), colormap=cm.viridis)
    plt.title('Number of Detected Reads per Taxon Across Thresholds')
    plt.xlabel('Taxa')
    plt.ylabel('Number of Reads')
    plt.xticks(rotation=45)
    plt.legend(title='Threshold', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("threshold_plot.png")


# Specify the directory path containing your files
directory_path = '.'

# Read and parse the data
df = read_and_parse_files(directory_path)

# Plot the data
plot_data(df)

