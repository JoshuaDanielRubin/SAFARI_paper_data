import matplotlib.pyplot as plt
import re

def parse_f1_scores(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = {}
    current_damage = ''
    for line in lines:
        if 'Damage:' in line:
            current_damage = line.split(':')[1].strip()
            data[current_damage] = {}
        elif 'Tool:' in line:
            current_tool = line.split(':')[1].strip().lower()
            data[current_damage][current_tool] = {'all': 0, 'mitochondrial': 0}
        elif 'F1 Score:' in line:
            score = float(line.split(':')[1].strip())
            if '[Mitochondrial Reads]' in prev_line:
                data[current_damage][current_tool]['mitochondrial'] = score
            elif '[All Reads]' in prev_line:
                data[current_damage][current_tool]['all'] = score
        prev_line = line
    return data

def calculate_percent_changes(data):
    percent_changes = {'all': {}, 'mitochondrial': {}}
    for damage in data:
        for read_type in ['all', 'mitochondrial']:
            safari_score = data[damage]['safari'][read_type]
            giraffe_score = data[damage]['giraffe'][read_type]
            if giraffe_score != 0:  # Check to avoid division by zero
                percent_change = ((safari_score - giraffe_score) / giraffe_score) * 100
            else:
                percent_change = None  # Or set to a default value or skip
            percent_changes[read_type][damage] = percent_change
    return percent_changes

def plot_percent_changes(percent_changes):
    fig, axs = plt.subplots(2, 1, figsize=(10, 8))
    damage_levels = ['none', 'single', 'dmid', 'dhigh']
    damage_labels = ['None', 'Single-stranded', 'Mid', 'High']

    for i, read_type in enumerate(['all', 'mitochondrial']):
        changes = [percent_changes[read_type].get(damage, 0) for damage in damage_levels]  # Ensure valid number, using 0 as default
        # Filter out None values to avoid TypeError, can also decide to skip plotting if None values are critical
        changes = [change for change in changes if change is not None]
        if changes:  # Check if there are valid changes to plot
            axs[i].bar(damage_labels[:len(changes)], changes, color=['skyblue', 'lightgreen'][i])
            axs[i].set_title(f'F1 Score Change ({read_type.capitalize()} Reads)')
            axs[i].set_ylabel('Percent Change')
            axs[i].set_xlabel('Damage Level')

    plt.tight_layout()
    plt.savefig("linear.png")

# Adjust the file path as necessary
file_path = 'linear_stats.txt'
f1_scores = parse_f1_scores(file_path)
percent_changes = calculate_percent_changes(f1_scores)
plot_percent_changes(percent_changes)

