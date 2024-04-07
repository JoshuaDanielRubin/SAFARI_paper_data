import matplotlib.pyplot as plt
import pandas as pd

# Data extracted from the file
data = {
    "Damage Coefficient": [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00],
    "hominin": [84809, 84809, 85472, 84455, 83729, 80843, 78597, 75429],
    "rCRS": [71539, 71539, 71605, 70610, 65091, 58719, 54407, 51097]
}

# Creating DataFrame
df = pd.DataFrame(data)

# Plotting adjustments for larger and bold text
plt.rcParams.update({'font.size': 14, 'font.weight': 'bold'})

# Adjusting subplots to be vertical and sharing x-axis, also adjusting for a wider figure
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 12), sharex=True)

df.plot(x='Damage Coefficient', y='hominin', ax=axes[0], marker='o', linestyle='-', color='blue', legend=False)
axes[0].set_title('Full Hominin Panmitogenome', fontsize=16, fontweight='bold')
axes[0].set_ylabel('Number Aligned', fontsize=14, fontweight='bold')

df.plot(x='Damage Coefficient', y='rCRS', ax=axes[1], marker='o', linestyle='-', color='red', legend=False)
axes[1].set_title('Single-haplotype Panmitogenome', fontsize=16, fontweight='bold')
axes[1].set_xlabel('Damage Coefficient', fontsize=14, fontweight='bold')
axes[1].set_ylabel('Number Aligned', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig('misspecified_vertical.png')

