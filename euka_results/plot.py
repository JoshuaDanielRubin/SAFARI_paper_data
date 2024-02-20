import pandas as pd
import matplotlib.pyplot as plt
import os

# Directory containing your files (adjust as needed)
directory = '/home/projects/MAAG/Magpie/Magpie/euka_results'
# List all files matching the pattern
files = [f for f in os.listdir(directory) if f.endswith('_detected.tsv')]

# Manually specify column names to avoid issues with special characters
columns = ['Taxa', 'Detected', 'Number_of_reads', 'Proportion_estimate',
           '85_CI_lower', '85_CI_upper', '95_CI_lower', '95_CI_upper']

# Initialize a dictionary to accumulate data
data = {}

# Process each file
for file in files:
    if 'uncorrected' in file:
        continue
    # Extract the threshold value from the filename
    threshold = float(file.split('_j_')[1].split('_detected')[0])
    # Construct full file path
    path = os.path.join(directory, file)
    # Read the TSV file, specifying column names to avoid issues with special characters
    df = pd.read_csv(path, sep='\t', comment='#', names=columns, header=0)
    # Iterate through the DataFrame rows
    for index, row in df.iterrows():
        taxon = row['Taxa']
        reads = row['Number_of_reads']
        if taxon not in data:
            data[taxon] = {}
        data[taxon][threshold] = reads

# Convert the accumulated data into a DataFrame for plotting
thresholds = sorted(list(set(threshold for d in data.values() for threshold in d)))
taxa = sorted(data.keys())
df_plot = pd.DataFrame(index=taxa, columns=thresholds).fillna(0)

# Populate the DataFrame with the accumulated data
for taxon, thresholds_data in data.items():
    for threshold, reads in thresholds_data.items():
        df_plot.at[taxon, threshold] = reads

# Sort taxa by their total reads at the smallest threshold
df_plot['Total'] = df_plot.sum(axis=1)
df_plot = df_plot.sort_values('Total', ascending=False)
del df_plot['Total']

# Plot
fig, ax = plt.subplots(figsize=(10, 8))
df_plot.plot(kind='bar', stacked=True, ax=ax, colormap='viridis')
plt.title('Read Counts by Taxon across Thresholds')
plt.xlabel('Taxon')
plt.ylabel('Number of Reads')
plt.xticks(rotation=45, ha='right')
plt.legend(title='Threshold', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the figure
plot_path = 'threshold_plot.png'
plt.savefig(plot_path)
plt.close()

print(f"Plot saved as {plot_path}")

