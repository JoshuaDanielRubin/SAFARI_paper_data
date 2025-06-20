import os
import re
from collections import defaultdict
import sys

# Directory containing the stat files
directory = sys.argv[1]

# Data structure to hold parsed data
results = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

# Regular expression to extract info from filenames and file contents
file_re = re.compile(r'numtS_and_gen_0_n300000_(d\w+)_l\d+_s([\d.]+)_(\w+)\.stat')
content_re = re.compile(r'(\w+ reads \w+): (\d+)')

# Function to calculate precision, recall, and F1
def calculate_metrics(correct, incorrect, unmapped):
    precision = correct / (correct + incorrect) if correct + incorrect else 0
    recall = correct / (correct + unmapped) if correct + unmapped else 0
    f1 = 2 * (precision * recall) / (precision + recall) if precision + recall else 0
    return precision, recall, f1

# Parse each file
for filename in os.listdir(directory):
    match = file_re.search(filename)
    if match:
        damage, similarity, tool = match.groups()
        with open(os.path.join(directory, filename), 'r') as file:
            for line in file:
                if 'Correct mitochondrial mappings' in line or 'Incorrect mitochondrial mappings' in line:
                    key, value = line.strip().split(': ')
                    results[damage][tool][key] = int(value)
                else:
                    content_match = content_re.search(line)
                    if content_match:
                        key, value = content_match.groups()
                        results[damage][tool][key] += int(value)

# Process and print the results
for damage in sorted(results):
    print(f'Damage: {damage}')
    for tool in sorted(results[damage]):
        data = results[damage][tool]
        correct_mito = data['Correct mitochondrial mappings']
        incorrect_mito = data['Incorrect mitochondrial mappings']
        mapped_mito = correct_mito  # Assuming correct mappings are all that's considered "mapped"
        unmapped_mito = data.get('Mito reads unmapped', 0)

        # Calculate metrics for mitochondrial reads
        mito_precision, mito_recall, mito_f1 = calculate_metrics(correct_mito, incorrect_mito, unmapped_mito)

        # Calculate overall metrics (assuming all reads)
        total_mapped = sum(data[k] for k in data if 'mapped' in k)
        total_unmapped = sum(data[k] for k in data if 'unmapped' in k)
        all_precision, all_recall, all_f1 = calculate_metrics(total_mapped - incorrect_mito, incorrect_mito, total_unmapped)

        # Print formatted results for mitochondrial reads and all reads
        print(f'Tool: {tool}')
        print(f'  [Mitochondrial Reads]\n    Precision: {mito_precision:.5f}\n    Recall: {mito_recall:.5f}\n    F1 Score: {mito_f1:.5f}')
        print('  [Raw Counts]')
        for key, value in data.items():
            print(f'    {key}: {value}')
        print(f'  [All Reads]\n    Precision: {all_precision:.5f}\n    Recall: {all_recall:.5f}\n    F1 Score: {all_f1:.5f}\n')
