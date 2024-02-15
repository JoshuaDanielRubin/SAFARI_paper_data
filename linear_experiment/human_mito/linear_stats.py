import csv
import numpy as np
import matplotlib.pyplot as plt

def calculate_metrics(TP, FP, TN, FN):
    sensitivity = TP / (TP + FN) if TP + FN != 0 else 0
    specificity = TN / (TN + FP) if TN + FP != 0 else 0
    precision = TP / (TP + FP) if TP + FP != 0 else 0
    accuracy = (TP + TN) / (TP + FP + TN + FN) if TP + FP + TN + FN != 0 else 0
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if precision + sensitivity != 0 else 0
    return f1_score, sensitivity, specificity, precision, accuracy

def rename_damage_level(level):
    mapping = {
        "none": "None",
        "dmid": "Mid",
        "dhigh": "High",
        "single": "Single-stranded"
    }
    return mapping.get(level, level)

def rename_aligner(aligner):
    mapping = {
        "safari": "SAFARI",
        "bb": "BBMap",
        "aln_anc": "BWA-aln (anc)",
        "aln": "BWA-aln",
        "mem": "BWA-MEM",
        "Bowtie2": "Bowtie2",
        "shrimp": "SHRiMP"
    }
    return mapping.get(aligner, aligner)

with open('alignment_stats.csv', 'r') as file:
    csv_reader = csv.DictReader(file)
    
    summary = {}
    
    for row in csv_reader:
        damage_type = rename_damage_level(row["Damage_Type"])
        aligner = rename_aligner(row["Aligner_Name"])
        
        if damage_type not in summary:
            summary[damage_type] = {}
        
        if aligner not in summary[damage_type]:
            summary[damage_type][aligner] = {
                "TP": 0, "FP": 0, "TN": 0, "FN": 0,
                "TP_MQ30": 0, "FP_MQ30": 0, "TN_MQ30": 0, "FN_MQ30": 0, "Support": 0
            }

        for key in ["TP", "FP", "TN", "FN", "TP_MQ30", "FP_MQ30", "TN_MQ30", "FN_MQ30"]:
            summary[damage_type][aligner][key] += int(row[key])

        # Add this line to calculate the support for each aligner and damage type:
        summary[damage_type][aligner]["Support"] = summary[damage_type][aligner]["TP"] + summary[damage_type][aligner]["FP"] + summary[damage_type][aligner]["TN"] + summary[damage_type][aligner]["FN"]


    aligner_order = ['SAFARI', 'giraffe'] + [aligner for aligner in summary[next(iter(summary))].keys() if aligner not in ['SAFARI', 'giraffe']]
    metrics = ["F1 Score", "Sensitivity", "Specificity", "Precision", "Accuracy"]
    
    for damage_type in summary:
        for metric_type, suffix in [("Overall", ""), ("MQ > 30", "_MQ30")]:
            data_to_plot = {metric: [] for metric in metrics}
            
            # Check which aligners have data for this damage type
            available_aligners = [aligner for aligner in aligner_order if aligner in summary[damage_type]]
            
            for aligner in available_aligners:
                results = calculate_metrics(summary[damage_type][aligner]["TP" + suffix], summary[damage_type][aligner]["FP" + suffix], \
                                            summary[damage_type][aligner]["TN" + suffix], summary[damage_type][aligner]["FN" + suffix])
                
                # Print metrics to terminal
                print(f"Metrics for Damage Type: {damage_type}, Aligner: {aligner}, Metric Type: {metric_type}")
                for metric, value in zip(metrics, results):
                    print(f"{metric}: {value:.4f}")
                print(f"Support: {summary[damage_type][aligner]['Support']}")
                print("-" * 50)

                
