import matplotlib.pyplot as plt
import pandas as pd

# Data for plotting, based on the table provided
data = {
    "Damage Coefficient": [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00],
    "vg giraffe Hominin": [101200, 88719, 77684, 68250, 60127, 53311, 47019, 41896],
    "SAFARI Hominin": [120660, 107364, 96179, 84455, 75319, 65623, 57320, 49846],
    "vg giraffe rCRS": [66992, 57406, 49214, 42445, 36638, 32076, 27681, 24244],
    "SAFARI rCRS": [104326, 92001, 80955, 70610, 58435, 48029, 39885, 33515]
}

# Convert data into a DataFrame
df = pd.DataFrame(data)

# Plotting
plt.figure(figsize=(10, 6))

# Plotting each line
plt.plot(df["Damage Coefficient"], df["vg giraffe Hominin"], color='orange', label='vg giraffe (Hominin graph)')
plt.plot(df["Damage Coefficient"], df["SAFARI Hominin"], color='green', label='SAFARI (Hominin graph)')
plt.plot(df["Damage Coefficient"], df["vg giraffe rCRS"], color='orange', linestyle='--', label='vg giraffe (rCRS graph)')
plt.plot(df["Damage Coefficient"], df["SAFARI rCRS"], color='green', linestyle='--', label='SAFARI (rCRS graph)')

# Adding titles and labels
plt.title('Aligned Reads: vg giraffe vs. SAFARI')
plt.xlabel('Damage Coefficient')
plt.ylabel('Number of Aligned Reads')
plt.legend()

# Show the plot
plt.tight_layout()
plt.savefig("hominin_results.png", dpi=300)

