
import re

def compute_metrics(TP, FP, TN, FN):
    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    specificity = TN / (TN + FP) if (TN + FP) > 0 else 0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    f1_score = (2 * precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
    return sensitivity, specificity, f1_score

def process_file(file_path):
    data_dict = {'hominin': {}, 'rCRS': {}}
    with open(file_path, 'r') as file:
        for line in file:
            if '.bam' in line:
                parts = line.strip().split('_')
                graph_type = parts[1]
                damage_coeff = float(parts[2].split('.bam')[0])
                if damage_coeff not in data_dict[graph_type]:
                    data_dict[graph_type][damage_coeff] = {'TP': 0, 'FP': 0, 'TN': 0, 'FN': 0}
            else:
                if 'Correct mitochondrial mappings:' in line:
                    data_dict[graph_type][damage_coeff]['TP'] = int(line.split(': ')[1])
                elif 'Incorrect mitochondrial mappings:' in line:
                    data_dict[graph_type][damage_coeff]['FP'] = int(line.split(': ')[1])
                elif 'Bacteria reads mapped:' in line:
                    data_dict[graph_type][damage_coeff]['FN'] = int(line.split(': ')[1])
                elif 'Bacteria reads unmapped:' in line:
                    data_dict[graph_type][damage_coeff]['TN'] = int(line.split(': ')[1])

    metrics_dict = {'hominin': {}, 'rCRS': {}}
    for graph_type in ['hominin', 'rCRS']:
        for damage_coeff, stats in data_dict[graph_type].items():
            TP, FP, TN, FN = stats['TP'], stats['FP'], stats['TN'], stats['FN']
            sensitivity, specificity, f1_score = compute_metrics(TP, FP, TN, FN)
            metrics_dict[graph_type][damage_coeff] = {
                'Sensitivity': sensitivity,
                'Specificity': specificity,
                'F1 Score': f1_score
            }
    return metrics_dict

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script.py <file_path>")
    else:
        file_path = sys.argv[1]
        metrics = process_file(file_path)
        print(metrics)

