import pandas as pd
import sys

def analyze_and_save(input_file_path, output_file_path):
    # Load the data
    data = pd.read_csv(input_file_path)
    
    # Calculate totals for bacteria, numt, and mito for each row
    data['bacteria_total'] = data['bacteria_mapped'] + data['bacteria_unmapped']
    data['numt_total'] = data['numt_mapped'] + data['numt_unmapped']
    data['mito_total_sum'] = data['mito_correct'] + data['mito_incorrect'] + data['mito_unmapped']

    # Filter for SAFARI and vg giraffe tools
    filtered_data = data[data['tool'].isin(['SAFARI', 'vg giraffe'])]

    # Ensure aggregation includes the median for totals for each (k, w) pair
    agg_operations = {
        'bacteria_mapped': 'median',
        'bacteria_unmapped': 'sum',
        'numt_mapped': 'median',
        'numt_unmapped': 'sum',
        'mito_correct': 'median',
        'mito_incorrect': 'sum',
        'mito_unmapped': 'sum',
        'mito_mapped': 'sum',
        'bacteria_total': 'median',  # Taking median of the pre-calculated total
        'numt_total': 'median',      # Ditto for numt
        'mito_total_sum': 'median'   # Ditto for mito
    }
    grouped = filtered_data.groupby(['k', 'w', 'tool']).agg(agg_operations).reset_index()

    # Pivot the table for side by side comparison, including total counts
    pivot_table = grouped.pivot_table(index=['k', 'w'], columns='tool', 
                                      values=['bacteria_mapped', 'bacteria_total', 
                                              'numt_mapped', 'numt_total', 
                                              'mito_correct', 'mito_total_sum']).reset_index()

    # Formatting for "mapped/total" or "correct/total"
    for metric in ['bacteria_mapped', 'numt_mapped', 'mito_correct']:
        for tool in ['SAFARI', 'vg giraffe']:
            total_col = metric.replace('mapped', 'total').replace('correct', 'total_sum')
            pivot_table[(metric, tool)] = pivot_table[(metric, tool)].astype(str) + '/' + pivot_table[(total_col, tool)].astype(str)

    # Define the columns order with updated formatting
    columns_order = [
        ('k', ''), ('w', ''),
        ('bacteria_mapped', 'SAFARI'), ('bacteria_mapped', 'vg giraffe'),
        ('numt_mapped', 'SAFARI'), ('numt_mapped', 'vg giraffe'),
        ('mito_correct', 'SAFARI'), ('mito_correct', 'vg giraffe')
    ]
    pivot_table = pivot_table.reindex(columns=pd.MultiIndex.from_tuples(columns_order))

    # Save the table to a CSV file
    pivot_table.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    input_file_path = sys.argv[1]  # The input CSV file name
    output_file_path = sys.argv[2]  # The output CSV file name
    analyze_and_save(input_file_path, output_file_path)

