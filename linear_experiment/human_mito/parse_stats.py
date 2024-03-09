import os
import csv
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score
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

def recall(y_true, y_pred):
    true_positives = np.sum((y_true == 1) & (y_pred == 1))
    actual_positives = np.sum(y_true == 1)
    if actual_positives == 0:
        return 0
    return true_positives / actual_positives

def f1_score(y_true, y_pred):
    p = precision(y_true, y_pred)
    r = recall(y_true, y_pred)
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

    total_samples = mito_correct + mito_incorrect + mito_unmapped + bacteria_mapped + bacteria_unmapped + numt_mapped + numt_unmapped

    y_true = [1] * (mito_correct + mito_incorrect + mito_unmapped) + [0] * (bacteria_mapped + bacteria_unmapped + numt_mapped + numt_unmapped)
    y_pred = [1] * (mito_correct + mito_incorrect + bacteria_mapped + numt_mapped) + [0] * (mito_unmapped + bacteria_unmapped + numt_unmapped)

    precision_score = precision(y_true, y_pred)
    recall_score = recall(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    accuracy_score = accuracy(y_true, y_pred)

    return {
        'precision': precision_score,
        'recall': recall_score,
        'f1': f1,
        'accuracy': accuracy_score
    }

    return {}

def process_subfolder(subfolder_path, csv_writer):
    print(subfolder_path)
    k, w = subfolder_path.split('/')[-1].split('_')
    k = int(k[1:]) if k != 'linear' else None
    w = int(w[1:]) if w != 'results' else None

    print(len(os.listdir(subfolder_path)))
    for i, file_name in enumerate(os.listdir(subfolder_path)):
        print(i)
        if file_name.endswith('.stat'):
            file_path = os.path.join(subfolder_path, file_name)
            tool_name = file_name.split('_')[-1].split('.')[0]
            damage_level = file_name.split('_')[3][1:]

            stats = parse_stat_file(file_path)
            metrics = compute_metrics(stats)

            row = [k, w, tool_name, damage_level]
            row.extend(stats.values())
            row.extend(metrics.values())
            csv_writer.writerow(row)

def main():
    output_file = 'results.csv'
    fieldnames = ['k', 'w', 'tool', 'damage_level', 'mito_correct', 'mito_incorrect', 'mito_mapped', 'mito_unmapped', 'bacteria_mapped', 'bacteria_unmapped', 'numt_mapped', 'numt_unmapped', 'precision', 'recall', 'f1', 'accuracy']

    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(fieldnames)

        alignments_dir = 'alignments'
        for subfolder in os.listdir(alignments_dir):
            subfolder_path = os.path.join(alignments_dir, subfolder)
            if os.path.isdir(subfolder_path) and (subfolder.startswith('k') or subfolder == 'linear_stats'):
                process_subfolder(subfolder_path, csv_writer)



if __name__ == '__main__':
    main()
