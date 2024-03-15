
import pandas as pd
import matplotlib.pyplot as plt

# Load data from CSV file
data = pd.read_csv('detected.csv')

# Sort the data by the total reads (SAFARI + vg giraffe)
data['total_reads'] = data['SAFARI_reads'] + data['vg_giraffe_reads']
data.sort_values(by='total_reads', ascending=False, inplace=True)

# Set up the bar chart
fig, ax = plt.subplots(figsize=(12, 8))
bar_width = 0.35
index = range(len(data))

# Plot the bars
bar1 = ax.bar(index, data['SAFARI_reads'], bar_width, label='SAFARI', color='g', edgecolor='grey')
bar2 = ax.bar([i + bar_width for i in index], data['vg_giraffe_reads'], bar_width, label='vg giraffe', color='orange', edgecolor='grey')

# Set the position and labels for the X-axis
ax.set_xlabel('Taxon', fontweight='bold')
ax.set_xticks([r + bar_width for r in range(len(data))])
ax.set_xticklabels(data['Taxon'], fontstyle='italic', rotation=45, ha='right')

# Add legend, title, and labels
ax.legend(fontsize='large')
ax.set_ylabel('Detected Reads')
ax.set_title('Comparison of Detected Reads in SAFARI vs. vg giraffe')

# Show the plot
plt.tight_layout()
plt.savefig("detected_taxa.png")

