import os
import re
import numpy as np
import glob
import sys

def parse_estimated(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip header
    matrix = [list(map(float, line.strip().split('\t'))) for line in lines]
    return np.array(matrix)

def parse_ground_truth(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip header
        print(lines)
    matrix = []
    for line in lines:
        parts = line.strip().split('\t')[1:]  # Skip first column
        row = [float(part.split(' ')[0]) for part in parts]  # Take only the first (numeric) part of each cell
        matrix.append(row)
    return np.array(matrix)

def compare_matrices(estimated, ground_truth):
    return np.sqrt(np.mean((estimated - ground_truth) ** 2))  # RMSE

def list_prof_files(directory):
    """Recursively list .prof files in a given directory that contain '_3' or '_5'."""
    return [os.path.join(dirpath, f)
            for dirpath, _, files in os.walk(directory)
            for f in files if f.endswith('.prof') and ('_3' in f or '_5' in f)]

def process_files(estimated_files, ground_truth_dir, output_file):
    rmse_results = {}
    for file_path in estimated_files:
        filename = os.path.basename(file_path)
        damage_type = re.search(r'_d(.*?)_', filename).group(1)
        tool_name = re.search(r'_s0\.\d+_(.*?)_', filename).group(1)
        strand = '3' if '_3' in filename else '5'
        for gt_file in ground_truth_files:
            if damage_type in gt_file and strand in gt_file:
                ground_truth_matrix = parse_ground_truth(gt_file)
                estimated_matrix = parse_estimated(file_path)
                print(ground_truth_matrix)
                print("\n\n\n")
                print(estimated_matrix)
                rmse = compare_matrices(estimated_matrix, ground_truth_matrix)
                category = (damage_type, tool_name, strand)
                
                if category not in rmse_results:
                    rmse_results[category] = []
                rmse_results[category].append(rmse)
    
    # Calculate average RMSE for each category
    for category, rmses in rmse_results.items():
        rmse_results[category] = np.median(rmses)

    # Save results to the output file
    with open(output_file, 'w') as out_file:
        for category, rmse in sorted(rmse_results.items()):
            out_file.write(f"Category: {category}, RMSE: {rmse}\n")

# Example Usage
vin_dir = sys.argv[1]
chag_dir = sys.argv[2]
output_file = sys.argv[3]

estimated_files = list_prof_files(vin_dir)
process_files(estimated_files, chag_dir, output_file)

