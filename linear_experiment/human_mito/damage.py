import os
import re
import numpy as np
import sys
import glob

# Avoid scientific notation
np.set_printoptions(suppress=True)

def parse_estimated(file_path):
    #print(file_path)
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip header
    matrix = [list(map(float, line.strip().split('\t'))) for line in lines]
    return np.array(matrix)

def parse_ground_truth(file_path, strand):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Skip header
    matrix = []
    for line in lines:
        parts = line.strip().split('\t')[1:]  # Skip first column
        row = [float(part.split(' ')[0]) for part in parts]  # Take only the first (numeric) part of each cell
        matrix.append(row)

    matrix = matrix[::-1]
    
    return np.array(matrix)


def compare_matrices(estimated, ground_truth):
    return np.sqrt(np.mean((estimated - ground_truth) ** 2))

def list_files(base_directory):
    # Find all subdirectories starting with 'k'
    subdirs = [d for d in os.listdir(base_directory) if os.path.isdir(os.path.join(base_directory, d)) and d.startswith('k')]
    all_files = []
    # Iterate through each matching subdir to find .prof files
    for subdir in subdirs:
        subdir_path = os.path.join(base_directory, subdir)
        # Use glob to find all .prof files within the subdir
        files = glob.glob(os.path.join(subdir_path, '*.prof'))
        all_files.extend(files)
    return all_files

def find_ground_truth_files(damage_type, ground_truth_dir, strand):
    pattern = f"{damage_type}{strand}.dat"
    for file_name in os.listdir(ground_truth_dir):
        if file_name == pattern:
            return os.path.join(ground_truth_dir, file_name)
    return None

def parse_k_w_from_path(file_path):
    match = re.search(r'k(\d+)_w(\d+)', file_path)
    return match.groups() if match else ('NA', 'NA')

def process_files(file_paths, ground_truth_dir):
    rmse_results = []
    for file_path in file_paths:
        filename = os.path.basename(file_path)
        if not ('_3' in filename or '_5' in filename):
            continue
        damage_type = re.search(r'_d(.*?)_', filename).group(1)
        tool_name = re.search(r'_s0\.\d+_(.*?)_', filename).group(1)
        sub_rate = re.search(r'_s(0\.\d+)_', filename).group(1)
        strand = '3' if '_3' in filename else '5'
        k, w = parse_k_w_from_path(file_path)
        #print(k,w)

        gt_file_path = find_ground_truth_files(damage_type, ground_truth_dir, strand)
        if not gt_file_path:
            print(f"No ground truth file found for {filename}")
            continue

        ground_truth_matrix = parse_ground_truth(gt_file_path, strand)
        estimated_matrix = parse_estimated(file_path)

        if damage_type == 'dhigh' and strand == '3':
            print('\n\n\n')
            print(tool_name)
            print(file_path)
            print(ground_truth_matrix)
            print('\n')
            print(estimated_matrix)
            print('\n\n\n')

        rmse = compare_matrices(estimated_matrix, ground_truth_matrix)

        # Check for 'vin' in file_path, default to 'N/A' if neither 'vin' nor 'Chagyrskaya' are present
        fragment_len_dist = 'Vindija' if 'vin' in file_path else ('Chagyrskaya' if 'chag' in file_path else 'N/A')

        # Check if 'safari' is in tool_name, default to 'N/A' if neither 'safari' nor 'vg giraffe' are present
        tool_name = 'SAFARI' if 'safari' in tool_name else ('vg giraffe' if 'giraffe' in tool_name else 'N/A')


        # Update damage levels for output
        damage_level = {'none': 'None', 'dmid': 'Mid', 'dhigh': 'High', 'single': 'Single-stranded'}[damage_type]

        rmse_results.append((fragment_len_dist, k, w, damage_level, tool_name, sub_rate, rmse))
    
    return rmse_results

def write_results_to_file(rmse_results, output_file_path):
    with open(output_file_path, 'w') as f:
        # Write header
        f.write('fragment_len_dist,k,w,damage_level,tool,sub_rate,RMSE\n')
        for result in rmse_results:
            f.write(','.join(map(str, result)) + '\n')

if __name__ == "__main__":
    vin_dir = sys.argv[1]
    chag_dir = sys.argv[2]
    ground_truth_dir = '.'  # Assuming current directory as per your request
    output_file = sys.argv[3]

    vin_files = list_files(vin_dir)
    chag_files = list_files(chag_dir)

    rmse_results_vin = process_files(vin_files, ground_truth_dir)
    rmse_results_chag = process_files(chag_files, ground_truth_dir)

    # Combine results from both directories and write to output
    combined_results = rmse_results_vin + rmse_results_chag
    write_results_to_file(combined_results, output_file)

