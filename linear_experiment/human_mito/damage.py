import os
import re
import numpy as np
import glob
import sys

def parse_estimated_matrix(file_path):
    matrix = np.zeros((4, 4))  # A, C, G, T
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            row_idx = "ACGT".find(parts[0])
            for i, count in enumerate(parts[1:]):
                matrix[row_idx, i] = int(count)
    total_counts = np.sum(matrix)
    if total_counts > 0:
        return matrix / total_counts
    return matrix

def parse_ground_truth(file_path):
    conversion = {'A>C': (0, 1), 'A>G': (0, 2), 'A>T': (0, 3),
                  'C>A': (1, 0), 'C>G': (1, 2), 'C>T': (1, 3),
                  'G>A': (2, 0), 'G>C': (2, 1), 'G>T': (2, 3),
                  'T>A': (3, 0), 'T>C': (3, 1), 'T>G': (3, 2)}
    matrix = np.zeros((4, 4))
    with open(file_path, 'r') as file:
        headers = next(file).strip().split()
        for line in file:
            parts = re.split(r'\s+', line.strip())
            for i, header in enumerate(headers):
                if header in conversion:
                    row_idx, col_idx = conversion[header]
                    prob_part = parts[i + 1].split('[')[0]
                    if prob_part == "0.0":
                        prob = 0.0
                    else:
                        try:
                            prob = float(prob_part)
                        except ValueError:
                            continue
                    matrix[row_idx, col_idx] = prob
    return matrix

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
    return glob.glob(directory + '/**/*.txt', recursive=True)

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

