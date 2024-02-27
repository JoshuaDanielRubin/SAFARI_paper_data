import matplotlib.pyplot as plt
import numpy as np
import re
from collections import defaultdict

def parse_data(file_path):
    data_struct = defaultdict(lambda: defaultdict(dict))
    current_damage, current_tool, read_type = None, None, None
    re_damage = re.compile(r'^Damage: (\w+)$')
    re_tool = re.compile(r'^Tool: (\w+)$')
    re_f1_score = re.compile(r'^\s+F1 Score: ([0-9.]+)$')

    with open(file_path, 'r') as file:
        for line in file:
            match_damage = re_damage.match(line)
            if match_damage:
                current_damage = match_damage.group(1)
                continue
            match_tool = re_tool.match(line)
            if match_tool:
                current_tool = match_tool.group(1)
                continue
            if '[Mitochondrial Reads]' in line:
                read_type = 'mitochondrial'
                continue
            elif '[All Reads]' in line:
                read_type = 'all'
                continue
            match_f1_score = re_f1_score.match(line)
            if match_f1_score and current_damage and current_tool and read_type:
                f1_score = float(match_f1_score.group(1))
                data_struct[current_damage][current_tool][read_type] = f1_score
    return data_struct

def calculate_percent_increases(data_struct):
    percent_increases = defaultdict(lambda: defaultdict(float))
    for damage, tools in data_struct.items():
        if 'giraffe' in tools and 'safari' in tools:
            for read_type in ['mitochondrial', 'all']:
                giraffe_f1 = tools['giraffe'].get(read_type, 0)
                safari_f1 = tools['safari'].get(read_type, 0)
                if giraffe_f1 > 0:  # Avoid division by zero
                    percent_increase = ((safari_f1 - giraffe_f1) / giraffe_f1) * 100
                    percent_increases[damage][read_type] = percent_increase
    return percent_increases

def plot_data(percent_increases):
    damage_types_ordered = ['none', 'single', 'dmid', 'dhigh']
    damage_names = ['None', 'Single-stranded', 'Mid', 'High']
    mitochondrial_increases_ordered = [percent_increases[damage]['mitochondrial'] for damage in damage_types_ordered]
    all_increases_ordered = [percent_increases[damage]['all'] for damage in damage_types_ordered]
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    y_positions_ordered = np.arange(len(damage_names))
    
    axs[0].barh(y_positions_ordered, mitochondrial_increases_ordered, color='skyblue')
    axs[0].set_yticks(y_positions_ordered)
    axs[0].set_yticklabels(damage_names, rotation=45)
    axs[0].set_title('Increase of F1 Score for Mitochondrial Reads')
    
    axs[1].barh(y_positions_ordered, all_increases_ordered, color='lightgreen')
    axs[1].set_yticks(y_positions_ordered)
    axs[1].set_yticklabels(damage_names, rotation=45)
    axs[1].set_xlabel('% Increase from Giraffe to SAFARI')
    axs[1].set_ylabel('Damage Matrix')
    axs[1].set_title('Increase of F1 Score for All Reads')
    
    plt.tight_layout()
    plt.savefig("linear2.png")

if __name__ == "__main__":
    file_path = 'linear_stats.txt'  # Update this path to your data file
    data_struct = parse_data(file_path)
    percent_increases = calculate_percent_increases(data_struct)
    plot_data(percent_increases)

