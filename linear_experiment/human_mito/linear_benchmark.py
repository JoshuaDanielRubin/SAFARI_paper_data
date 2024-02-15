import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_data(file_path):
    """Load CSV data from the specified file path."""
    return pd.read_csv(file_path)

def calculate_metrics(df):
    """Calculate sensitivity, specificity, and their MQ30 equivalents."""
    sensitivity = df['TP'] / (df['TP'] + df['FN'])
    specificity = df['TN'] / (df['TN'] + df['FP'])
    sensitivity_mq30 = df['TP_MQ30'] / (df['TP_MQ30'] + df['FN_MQ30'])
    specificity_mq30 = df['TN_MQ30'] / (df['TN_MQ30'] + df['FP_MQ30'])
    return sensitivity.mean(), specificity.mean(), sensitivity_mq30.mean(), specificity_mq30.mean()

def plot_metrics(safari_metrics, giraffe_metrics, safari_fpr, giraffe_fpr):
    """Plot sensitivity and false positive rate for SAFARI and vg giraffe."""
    labels = ['Overall', 'MQ30']
    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Sensitivity
    ax[0].bar(x - width/2, safari_metrics[:2], width, color='green', label='SAFARI')
    ax[0].bar(x + width/2, giraffe_metrics[:2], width, color='orange', label='vg giraffe')
    ax[0].set_ylabel('Sensitivity')
    ax[0].set_title('Sensitivity Comparison')
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(labels)
    ax[0].legend()

    # False Positive Rate
    ax[1].bar(x - width/2, safari_fpr, width, color='green', label='SAFARI')
    ax[1].bar(x + width/2, giraffe_fpr, width, color='orange', label='vg giraffe')
    ax[1].set_ylabel('1 - Specificity (FPR)')
    ax[1].set_title('False Positive Rate Comparison')
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(labels)
    ax[1].legend()

    plt.tight_layout()
    plt.show()

def main(file_path):
    data = load_data(file_path)

    # Filter data for SAFARI and vg giraffe
    safari_data = data[data['Aligner_Name'] == 'safari']
    giraffe_data = data[data['Aligner_Name'] == 'giraffe']

    # Calculate metrics
    safari_sensitivity, safari_specificity, safari_sensitivity_mq30, safari_specificity_mq30 = calculate_metrics(safari_data)
    giraffe_sensitivity, giraffe_specificity, giraffe_sensitivity_mq30, giraffe_specificity_mq30 = calculate_metrics(giraffe_data)

    safari_metrics = [safari_sensitivity, safari_sensitivity_mq30, safari_specificity, safari_specificity_mq30]
    giraffe_metrics = [giraffe_sensitivity, giraffe_sensitivity_mq30, giraffe_specificity, giraffe_specificity_mq30]

    # Calculate the complement of specificity (False Positive Rate)
    safari_fpr = [1 - safari_specificity, 1 - safari_specificity_mq30]
    giraffe_fpr = [1 - giraffe_specificity, 1 - giraffe_specificity_mq30]

    # Plot
    plot_metrics(safari_metrics, giraffe_metrics, safari_fpr, giraffe_fpr)

if __name__ == "__main__":
    file_path = 'alignment_stats.csv'  # Update this with the path to your CSV file
    main(file_path)

