import re
import matplotlib.pyplot as plt
import numpy as np

def parse_threshold(header):
    match = re.search(r'_j_(\d+\.\d+)_detected', header)
    return float(match.group(1)) if match else None

def load_data(file_path):
    data_by_taxon = {}
    with open(file_path, 'r') as file:
        current_threshold = None
        for line in file:
            if line.startswith('euka_corrected_full_j_'):
                current_threshold = parse_threshold(line.strip())
            elif not line.startswith('#') and line.strip():
                parts = line.split('\t')
                print(parts)
                taxon = parts[0]
                number_of_reads = int(parts[2])
                if taxon not in data_by_taxon:
                    data_by_taxon[taxon] = {}
                data_by_taxon[taxon][current_threshold] = number_of_reads
    return data_by_taxon

def calculate_incremental_reads(data_by_taxon):
    incremental_data_by_taxon = {}
    for taxon, thresholds_reads in data_by_taxon.items():
        sorted_thresholds = sorted(thresholds_reads.keys(), reverse=True)
        incremental_reads = []
        previous_reads = 0
        for threshold in sorted_thresholds:
            current_reads = thresholds_reads[threshold]
            incremental_reads.append((threshold, current_reads - previous_reads))
            previous_reads = current_reads
        incremental_data_by_taxon[taxon] = incremental_reads
    return incremental_data_by_taxon

def plot_data(incremental_data_by_taxon, data_by_taxon):
    sorted_taxa = sorted(incremental_data_by_taxon.keys(), key=lambda x: data_by_taxon[x][0.0], reverse=True)
    cmap = plt.get_cmap("jet")
    thresholds = sorted(list(next(iter(incremental_data_by_taxon.values()))), key=lambda x: x[0])
    colors = {t[0]: cmap(i / len(thresholds)) for i, t in enumerate(thresholds)}

    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    bar_width = 1.0

    for i, taxon in enumerate(sorted_taxa):
        heights = np.array([max(0, val[1]) for val in incremental_data_by_taxon[taxon]])
        bottoms = np.cumsum(heights) - heights
        print(taxon, bottoms)
        for j, (threshold, _) in enumerate(incremental_data_by_taxon[taxon]):
            ax.bar(i, heights[j], bar_width, bottom=bottoms[j], color=colors[threshold])

      # Add a horizontal dashed line at y=100
    ax.axhline(y=100, color='green', linestyle='--')

    ax.set_xticks(range(len(sorted_taxa)))
    ax.set_xticklabels(sorted_taxa, rotation=90, fontstyle='italic')
    ax.set_ylabel('Number of Detected Reads')
    ax.set_title('Detected Reads per Taxon by Posterior Odds Threshold', pad=20)
    ax.legend([plt.Rectangle((0,0),1,1, color=colors[threshold[0]]) for threshold in thresholds], [f"Threshold {threshold[0]}" for threshold in thresholds], loc='upper left', bbox_to_anchor=(1,1))

    plt.tight_layout()
    plt.savefig("post_threshold_plot.png")

file_path = 'threshold_data.txt'  # Update with the path to your file
data_by_taxon = load_data(file_path)
incremental_data_by_taxon = calculate_incremental_reads(data_by_taxon)
plot_data(incremental_data_by_taxon, data_by_taxon)

