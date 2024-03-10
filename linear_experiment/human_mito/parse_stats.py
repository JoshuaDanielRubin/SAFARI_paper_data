import os
import csv
import numpy as np
import re

def parse_stat_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    mito_correct = 0
    mito_incorrect = 0
    mito_mapped = 0
    mito_unmapped = 0
    bacteria_mapped = 0
    bacteria_unmapped = 0
    numt_mapped = 0
    numt_unmapped = 0

    for line in lines:
        if line.startswith('Correct mitochondrial mappings:'):
            mito_correct = int(line.split(':')[1].strip())
        elif line.startswith('Incorrect mitochondrial mappings:'):
            mito_incorrect = int(line.split(':')[1].strip())
        elif line.startswith('Mito reads mapped:'):
            mito_mapped = int(line.split(':')[1].strip())
        elif line.startswith('Mito reads unmapped:'):
            mito_unmapped = int(line.split(':')[1].strip())
        elif line.startswith('Bacteria reads mapped:'):
            bacteria_mapped = int(line.split(':')[1].strip())
        elif line.startswith('Bacteria reads unmapped:'):
            bacteria_unmapped = int(line.split(':')[1].strip())
        elif line.startswith('NuMT reads mapped:'):
            numt_mapped = int(line.split(':')[1].strip())
        elif line.startswith('NuMT reads unmapped:'):
            numt_unmapped = int(line.split(':')[1].strip())

    return {
        'mito_correct': mito_correct,
        'mito_incorrect': mito_incorrect,
        'mito_mapped': mito_mapped,
        'mito_unmapped': mito_unmapped,
        'bacteria_mapped': bacteria_mapped,
        'bacteria_unmapped': bacteria_unmapped,
        'numt_mapped': numt_mapped,
        'numt_unmapped': numt_unmapped
    }

def compute_confusion_matrix_elements(y_true, y_pred):
    TP = np.sum((y_true == 1) & (y_pred == 1))
    TN = np.sum((y_true == 0) & (y_pred == 0))
    FP = np.sum((y_true == 0) & (y_pred == 1))
    FN = np.sum((y_true == 1) & (y_pred == 0))
    return TP, FP, TN, FN

def precision(TP, FP):
    if TP + FP == 0:
        return 0
    return TP / (TP + FP)

def sensitivity(TP, FN):  # Also known as recall
    if TP + FN == 0:
        return 0
    return TP / (TP + FN)

def specificity(TN, FP):
    if TN + FP == 0:
        return 0
    return TN / (TN + FP)

def f1_score(precision, sensitivity):
    if precision == 0 and sensitivity == 0:
        return 0
    return 2 * (precision * sensitivity) / (precision + sensitivity)

def accuracy(TP, TN, FP, FN):
    return (TP + TN) / (TP + TN + FP + FN)

def compute_metrics(stats):
    mito_correct = stats['mito_correct']
    mito_incorrect = stats['mito_incorrect']
    mito_unmapped = stats['mito_unmapped']
    bacteria_mapped = stats['bacteria_mapped']
    bacteria_unmapped = stats['bacteria_unmapped']
    numt_mapped = stats['numt_mapped']
    numt_unmapped = stats['numt_unmapped']

    y_true = np.array([1] * (mito_correct + mito_unmapped) + [0] * (mito_incorrect + bacteria_mapped + bacteria_unmapped + numt_mapped + numt_unmapped))
    y_pred = np.array([1] * (mito_correct + mito_incorrect + bacteria_mapped + numt_mapped) + [0] * (mito_unmapped + bacteria_unmapped + numt_unmapped))

    TP, FP, TN, FN = compute_confusion_matrix_elements(y_true, y_pred)
    precision_score = precision(TP, FP)
    sensitivity_score = sensitivity(TP, FN)
    specificity_score = specificity(TN, FP)
    f1 = f1_score(precision_score, sensitivity_score)
    accuracy_score = accuracy(TP, TN, FP, FN)

    return {
        'precision': precision_score,
        'sensitivity': sensitivity_score,
        'specificity': specificity_score,
        'f1': f1,
        'accuracy': accuracy_score,
        'TP': TP,
        'FP': FP,
        'TN': TN,
        'FN': FN
    }


def process_subfolder(subfolder_path, csv_writer):
    k, w = subfolder_path.split('/')[-1].split('_')
    k = int(k[1:]) if k != 'linear' else None
    w = int(w[1:]) if w != 'results' else None

    for i, file_name in enumerate(os.listdir(subfolder_path)):
        if file_name.endswith('.stat'):
            print(i)
            file_path = os.path.join(subfolder_path, file_name)
            tool_name = file_name.split('_')[-1].split('.')[0]
            damage_level = file_name.split('_')[5][1:]

            # Use regex to find read length (lXX) and subsampling rate (sX.X or sX)
            read_length_match = re.search(r'_l(\d+)_', file_name)
            subsampling_rate_match = re.search(r'_s([0-9]+\.?[0-9]*)_', file_name)
            read_length = int(read_length_match.group(1)) if read_length_match else None
            subsampling_rate = float(subsampling_rate_match.group(1)) if subsampling_rate_match else None

            stats = parse_stat_file(file_path)
            metrics = compute_metrics(stats)

            row = [k, w, tool_name, damage_level, read_length, subsampling_rate]
            row.extend(stats.values())
            row.extend(metrics.values())
            csv_writer.writerow(row)


def main():
    output_file = 'newresults.csv'
    fieldnames = ['k', 'w', 'tool', 'read_length', 'subsampling_rate', 'damage_level', 'mito_correct', 'mito_incorrect', 'mito_mapped', 'mito_unmapped', 'bacteria_mapped', 'bacteria_unmapped', \
                  'numt_mapped', 'numt_unmapped', 'precision', 'sensitivity', 'specificity', 'f1', 'accuracy', 'TP', 'FP', 'TN', 'FN']

    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(fieldnames)

        alignments_dir = 'alignments'
        for subfolder in os.listdir(alignments_dir):
            subfolder_path = os.path.join(alignments_dir, subfolder)
            if os.path.isdir(subfolder_path) and (subfolder.startswith('k') or subfolder == 'linear_results'):
                process_subfolder(subfolder_path, csv_writer)

if __name__ == '__main__':
    main()

