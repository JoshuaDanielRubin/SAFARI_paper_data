import os
import re
import numpy as np
import glob

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
        headers = next(file).strip().split()  # Read the first line to get headers
        for line in file:
            parts = re.split(r'\s+', line.strip())
            for i, header in enumerate(headers):
                if header in conversion:
                    row_idx, col_idx = conversion[header]
                    prob_part = parts[i + 1].split('[')[0]
                    if prob_part == "0.0":  # Directly handle the common case
                        prob = 0.0
                    else:
                        try:
                            prob = float(prob_part)
                        except ValueError:
                            # Skip or handle the specific range notation or log a more precise error
                            continue  # For now, just skip to avoid logging errors for [0..0]
                    matrix[row_idx, col_idx] = prob
    return matrix

def compare_matrices(estimated, ground_truth):
    return np.sqrt(np.mean((estimated - ground_truth) ** 2))  # RMSE

def list_files(directory):
    """List full file paths in a given directory."""
    return [os.path.join(directory, f) for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def find_ground_truth_file(damage_type, strand):
    """Map damage type and strand to the correct ground truth file."""
    gt_mapping = {
        ('dmid', '3p'): 'dmid3.dat',
        ('dmid', '5p'): 'dmid5.dat',
        ('dhigh', '3p'): 'dhigh3.dat',
        ('dhigh', '5p'): 'dhigh5.dat',
        ('none', '3p'): 'none3.dat',
        ('none', '5p'): 'none5.dat',
        ('single', '3p'): 'single3.dat',
        ('single', '5p'): 'single5.dat',
    }
    return gt_mapping.get((damage_type, strand), None)

def process_files(file_paths, ground_truth_dir):
    rmse_results = {}
    for file_path in file_paths:
        filename = os.path.basename(file_path)
        damage_type = re.search(r'_d(.*?)_', filename).group(1)
        tool_name = re.search(r'_s0\.\d+_(.*?)_', filename).group(1)
        strand = '3p' if '3p' in filename else '5p'
        
        ground_truth_file = find_ground_truth_file(damage_type, strand)
        if ground_truth_file:
            category = (damage_type, tool_name, strand)
            ground_truth_path = os.path.join(ground_truth_dir, ground_truth_file)
            ground_truth_matrix = parse_ground_truth(ground_truth_path)
            estimated_matrix = parse_estimated_matrix(file_path)
            rmse = compare_matrices(estimated_matrix, ground_truth_matrix)
            
            if category not in rmse_results:
                rmse_results[category] = []
            rmse_results[category].append(rmse)
    
    # Calculate average RMSE for each category
    for category, rmses in rmse_results.items():
        rmse_results[category] = np.median(rmses)
    
    return rmse_results

# Directory Paths
alignments_dir = 'alignments/profs/'
ground_truth_dir = './'  # Adjust as necessary

# List all substitution matrix files
file_paths = list_files(alignments_dir)

# Process files and calculate RMSE
rmse_results = process_files(file_paths, ground_truth_dir)

# Display Results
for category, rmse in sorted(rmse_results.items()):
    print(f"Category: {category}, RMSE: {rmse}")

