import gzip
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def read_gz_file(file_path):
    with gzip.open(file_path, 'rt') as f:  # 'rt' mode for text mode reading
        data = [int(line.strip()) for line in f if line.strip().isdigit()]
    return np.array(data)

def plot_density(distributions, labels, output_path):
    colors=["blue", "green"]
    fig, axs = plt.subplots(len(distributions), 1, figsize=(10, 8))
    x_range = np.linspace(0, max(max(d) for d in distributions), 1000)
    
    for i, (dist, label) in enumerate(zip(distributions, labels)):
        kde = gaussian_kde(dist)
        density = kde(x_range)
        axs[i].plot(x_range, density, color=colors[i])
        axs[i].fill_between(x_range, density, color=colors[i], alpha=0.3)
        axs[i].set_title(f'Density Plot of Fragment Lengths for {label}')
        axs[i].set_xlabel('Fragment Length')
        axs[i].set_ylabel('Density')
    
    plt.tight_layout()
    plt.savefig(output_path)

# Paths to the input files
file_paths = ["chagyrskaya8.gz", "Vi33.19.gz"]
labels = ["Chagyrskaya Cave", "Vindija Cave"]
distributions = [read_gz_file(fp) for fp in file_paths]

# Generate and save the plots
plot_density(distributions, labels, "dists.png")

