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

# Plotting
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6), sharey=True)

df.plot(x='Damage Coefficient', y='hominin', ax=axes[0], marker='o', linestyle='-', color='blue', legend=False)
axes[0].set_title('Full Hominin Panmitogenome')
axes[0].set_xlabel('Damage Coefficient')
axes[0].set_ylabel('Number Aligned')

df.plot(x='Damage Coefficient', y='rCRS', ax=axes[1], marker='o', linestyle='-', color='red', legend=False)
axes[1].set_title('Single-haplotype Panmitogenome')
axes[1].set_xlabel('Damage Coefficient')

plt.tight_layout()
plt.savefig("misspecified_dm.png")

