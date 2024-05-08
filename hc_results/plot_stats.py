import matplotlib.pyplot as plt

# Data for Mean Levenshtein Distance
methods = ['SAFARI', 'vg giraffe']
leven_dist = [6.2375, 6.78]

# Data for box plot
box_data = [
    [82.75, 160.5, 310.0],  # SAFARI: lower quartile, median, upper quartile
    [73.75, 148.0, 289.75]  # vg giraffe: lower quartile, median, upper quartile
]

# Create figure and axes with a more suitable aspect ratio
fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=300)

# Boxplot for Aligned and Processed Reads
axs[0].boxplot(box_data, vert=False, whis=[0, 100], widths=0.6, patch_artist=True,
               medianprops=dict(color="orange"), boxprops=dict(facecolor='lightgrey'))
axs[0].set_yticklabels(['SAFARI', 'vg giraffe'], fontweight='bold', fontsize=12)
axs[0].set_title('Aligned and Processed Reads', fontweight='bold', fontsize=14)
axs[0].set_xlabel('Read Count', fontweight='bold', fontsize=12)

# Bar plot for Mean Levenshtein Distance
bars = axs[1].bar(methods, leven_dist, color=['green', 'orange'])
axs[1].set_title('Mean Levenshtein Distance', fontweight='bold', fontsize=14)
axs[1].set_ylim(0, 8)
axs[1].set_xticklabels(methods, fontweight='bold', fontsize=12)  # Bolded x-axis labels for clarity
for bar, value in zip(bars, leven_dist):
    axs[1].text(bar.get_x() + bar.get_width() / 2, value + 0.1, f'{value}', ha='center', va='bottom', fontweight='bold', fontsize=10)

# Layout and display
plt.tight_layout()
plt.show()

