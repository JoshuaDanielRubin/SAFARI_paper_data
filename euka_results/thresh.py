import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def parse_filename_for_threshold(filename):
    # This regex assumes the threshold is always formatted as a decimal in the filename.
    match = re.search(r'_j_(\d+(\.\d+)?)', filename)
    return float(match.group(1)) if match else None

def read_data(file_path):
    # Open the file and read the first line that starts with '#' for the header
    with open(file_path, 'r') as file:
        header_line = next((line for line in file if line.startswith('#')), None)
    if header_line:
        # Remove the '#' and split the header line into column names
        header_line = header_line[1:].strip()  # Remove '#' and trailing newline
        column_names = header_line.split('\t')  # Assuming the file is tab-delimited

        # Read the file into a DataFrame, skipping the commented header line, using the extracted column names
        return pd.read_csv(file_path, sep='\t', comment='\0', names=column_names, header=0)
    else:
        # Return an empty DataFrame or raise an error if no header line was found
        return pd.DataFrame()

def process_files(directory):
    data = {}
    for filename in os.listdir(directory):
        if filename.endswith(".tsv") and 'detected' in filename:
            threshold = parse_filename_for_threshold(filename)
            if threshold is not None:
                file_path = os.path.join(directory, filename)
                df = read_data(file_path)
                for index, row in df.iterrows():
                    taxon = row['Taxa']
                    if taxon not in data:
                        data[taxon] = []
                    data[taxon].append((threshold, row['Number_of_reads']))
    
    # Sort data by threshold for each taxon
    for taxon in data.keys():
        data[taxon].sort(key=lambda x: x[0])
    
    return data

def plot_data(data):
    for taxon, values in data.items():
        thresholds, reads = zip(*values)
        plt.figure(figsize=(10, 6))
        plt.plot(thresholds, reads, marker='o', linestyle='-')
        plt.title(f"Detected Reads vs. Threshold for {taxon}")
        plt.xlabel('Threshold')
        plt.ylabel('Number of Detected Reads')
        plt.grid(True)
        plt.savefig("threshold_plot.png")

# Assuming the script runs in the directory containing the TSV files
directory = '.'  # Use the current directory
data = process_files(directory)
plot_data(data)

