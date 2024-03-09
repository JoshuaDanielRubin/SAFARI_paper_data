import os
import csv
from sklearn.metrics import precision_score, f1_score, accuracy_score
import numpy as np

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

def precision(y_true, y_pred):
    true_positives = np.sum((y_true == 1) & (y_pred == 1))
    predicted_positives = np.sum(y_pred == 1)
    if predicted_positives == 0:
        return 0
    return true_positives / predicted_positives

def sensitivity(y_true, y_pred): # Renamed from recall to sensitivity
    true_positives = np.sum((y_true == 1) & (y_pred == 1))
    actual_positives = np.sum(y_true == 1)
    if actual_positives == 0:
        return 0
    return true_positives / actual_positives

def specificity(y_true, y_pred):
    true_negatives = np.sum((y_true == 0) & (y_pred == 0))
    actual_negatives = np.sum(y_true == 0)
    if actual_negatives == 0:
        return 0
    return true_negatives / actual_negatives

def f1_score(y_true, y_pred):
    p = precision(y_true, y_pred)
    r = sensitivity(y_true, y_pred) # Updated to use sensitivity
    if p == 0 and r == 0:
        return 0
    return 2 * (p * r) / (p + r)

def accuracy(y_true, y_pred):
    return np.sum(y_true == y_pred) / len(y_true)

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

    precision_score = precision(y_true, y_pred)
    sensitivity_score = sensitivity(y_true, y_pred)
    specificity_score = specificity(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    accuracy_score = accuracy(y_true, y_pred)

    return {
        'precision': precision_score,
        'sensitivity': sensitivity_score,
        'specificity': specificity_score,
        'f1': f1,
        'accuracy': accuracy_score
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

            stats = parse_stat_file(file_path)
            metrics = compute_metrics(stats)

            row = [k, w, tool_name, damage_level]
            row.extend(stats.values())
            row.extend(metrics.values())
            csv_writer.writerow(row)

def main():
    output_file = 'results.csv'
    fieldnames = ['k', 'w', 'tool', 'damage_level', 'mito_correct', 'mito_incorrect', 'mito_mapped', 'mito_unmapped', 'bacteria_mapped', 'bacteria_unmapped', 'numt_mapped', 'numt_unmapped', 'precision', 'sensitivity', 'specificity', 'f1', 'accuracy']

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

