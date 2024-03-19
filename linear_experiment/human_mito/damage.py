import os
import re
import numpy as np
import glob
import sys

def parse_estimated_matrix(data_str):
    """
    Parses the estimated matrix from a .prof file.
    """
    # Assuming data_str is a path to the .prof file for simplicity
    with open(data_str, 'r') as file:
        lines = file.readlines()[1:]  # Skip header line
    
    # Parse numerical values
    matrix = [list(map(float, line.strip().split('\t'))) for line in lines]
    
    return np.array(matrix)

def parse_ground_truth(file_path):
    """
    Parses the ground truth matrix from a .dat file.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]  # Assuming the first line is a header and should be ignored for consistency
    
    # Extract numerical values and ignore the ranges for now, focusing on the primary values
    matrix = []
    for line in lines:
        row_values = []
        for value in line.strip().split('\t')[1:]:  # Skip the first item which is not a numerical value
            primary_value = float(value.split(' ')[0])
            row_values.append(primary_value)
        matrix.append(row_values)
    
    return np.array(matrix)

def average_matrices(matrix_list):
    return np.mean(matrix_list, axis=0)

def find_ground_truth_files(damage_type):
    gt_mapping = {
        'dmid': ('dmid3.dat', 'dmid5.dat'),  # Corrected file extensions
        'dhigh': ('dhigh3.dat', 'dhigh5.dat'),
        'none': ('none3.dat', 'none5.dat'),
        'single': ('single3.dat', 'single5.dat'),
    }
    return gt_mapping.get(damage_type, (None, None))

def process_files(file_paths, ground_truth_dir):
    rmse_results = []
    ground_truth_matrices = {}
    
    for damage_type in ['dmid', 'dhigh', 'none', 'single']:
        gt_files = find_ground_truth_files(damage_type)
        matrices = [parse_ground_truth(os.path.join(ground_truth_dir, gt_file)) for gt_file in gt_files if gt_file]
        ground_truth_matrices[damage_type] = average_matrices(matrices)
    
    for file_path in file_paths:
        print(file_path)
        if '_3' not in file_path and '_5' not in file_path:
            continue
        path_parts = file_path.split(os.sep)
        fragment_dist = path_parts[-3]  
        if 'linear_results' in file_path:
            k, w = 'linear', 'linear'
        else:
            match = re.search(r'k(\d+)_w(\d+)', path_parts[-2])
            if match:
                k, w = match.groups()
            else:
                print(f"File does not match expected pattern: {file_path}")
                continue
        
        filename = os.path.basename(file_path)
        damage_type = re.search(r'_d(.*?)_', filename).group(1)
        tool_name = re.search(r'_s0\.\d+_(.*?)_', filename).group(1)
        
        estimated_matrix = parse_estimated_matrix(file_path)
        ground_truth_matrix = ground_truth_matrices[damage_type]
        rmse = compare_matrices(estimated_matrix, ground_truth_matrix)
        
        rmse_results.append((fragment_dist, k, w, damage_type, tool_name, rmse))
    
    return rmse_results

def compare_matrices(estimated, ground_truth):
    return np.sqrt(np.mean((estimated - ground_truth) ** 2))

def list_files(directory):
    return glob.glob(directory + '/**/*.prof', recursive=True)

def write_results_to_file(rmse_results, output_file_path):
    with open(output_file_path, 'w') as f:
        for result in rmse_results:
            f.write(','.join(map(str, result)) + '\n')

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python damage.py <alignments_dir_vin> <alignments_dir_chag> <output_file_path>")
        sys.exit(1)
    
    alignments_dir_vin = sys.argv[1]
    alignments_dir_chag = sys.argv[2]
    output_file_path = sys.argv[3]
    ground_truth_dir = './'  # Adjust as necessary

    vin_file_paths = list_files(alignments_dir_vin)
    chag_file_paths = list_files(alignments_dir_chag)
    all_file_paths = vin_file_paths + chag_file_paths

    rmse_results = process_files(all_file_paths, ground_truth_dir)
    write_results_to_file(rmse_results, output_file_path)

    print(f"Results written to {output_file_path}")

