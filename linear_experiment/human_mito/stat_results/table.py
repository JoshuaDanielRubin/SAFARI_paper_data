import pandas as pd
import sys

def analyze_and_save(input_file_path, output_file_path):
    # Load the data
    data = pd.read_csv(input_file_path)
    
    # Filter for subsampling_rate=0.9
    data = data[data['subsampling_rate'] == 0.9]
    
    # Calculate totals for bacteria, numt, and mito for each row
    data['bacteria_total'] = data['bacteria_mapped'] + data['bacteria_unmapped']
    data['numt_total'] = data['numt_mapped'] + data['numt_unmapped']
    data['mito_total_sum'] = data['mito_correct'] + data['mito_incorrect'] + data['mito_unmapped']

    # Filter for SAFARI and vg giraffe tools
    filtered_data = data[data['tool'].isin(['SAFARI', 'vg giraffe'])]

    # Since there's only one sample per tool per (k, w) pair, aggregate operations are not needed for medians
    grouped = filtered_data.groupby(['k', 'w', 'tool']).first().reset_index()

    # Pivot the table for side by side comparison, including total counts
    pivot_table = grouped.pivot_table(index=['k', 'w'], columns='tool', 
                                      values=['bacteria_mapped', 'bacteria_total', 
                                              'numt_mapped', 'numt_total', 
                                              'mito_correct', 'mito_total_sum']).reset_index()

    # Convert to integers and format for "mapped/total" or "correct/total"
    for metric in ['bacteria_mapped', 'numt_mapped', 'mito_correct']:
        for tool in ['SAFARI', 'vg giraffe']:
            total_col = metric.replace('mapped', 'total').replace('correct', 'total_sum')
            pivot_table[(metric, tool)] = pivot_table[(metric, tool)].astype(int).astype(str) + '/' + pivot_table[(total_col, tool)].astype(int).astype(str)

    # Save the table to a LaTeX file
    with open(output_file_path, 'w') as latex_file:
        latex_file.write(pivot_table.to_latex(index=False, escape=False))

if __name__ == "__main__":
    input_file_path = sys.argv[1]  # The input CSV file name
    output_file_path = sys.argv[2]  # The output LaTeX file name
    analyze_and_save(input_file_path, output_file_path)

