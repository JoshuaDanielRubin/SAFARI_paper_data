import pandas as pd
import sys

def analyze_and_save(input_file_path, output_file_path):
    # Load the data
    data = pd.read_csv(input_file_path)
    
    # Add a new column for the sum of mito-related counts
    data['mito_total_sum'] = data['mito_correct'] + data['mito_incorrect'] + data['mito_unmapped']

    # Filter for SAFARI and vg giraffe tools
    filtered_data = data[data['tool'].isin(['SAFARI', 'vg giraffe'])]

    # Ensure aggregation includes the median for mito_total_sum for each (k, w) pair
    agg_operations = {
        'bacteria_mapped': 'median',
        'bacteria_unmapped': 'sum',
        'numt_mapped': 'median',
        'numt_unmapped': 'sum',
        'mito_correct': 'median',
        'mito_incorrect': 'sum',
        'mito_unmapped': 'sum',
        'mito_mapped': 'sum',
        'mito_total_sum': 'median'  # This already ensures the median for each (k, w) pair
    }
    grouped = filtered_data.groupby(['k', 'w', 'tool']).agg(agg_operations).reset_index()

    # Calculate total counts
    grouped['bacteria_total'] = grouped['bacteria_mapped'] + grouped['bacteria_unmapped']
    grouped['numt_total'] = grouped['numt_mapped'] + grouped['numt_unmapped']
    # Update mito_total to use mito_total_sum directly as it is the median value we need
    grouped['mito_total'] = grouped['mito_total_sum']

    # Pivot the table for side by side comparison, including total counts
    pivot_table = grouped.pivot_table(index=['k', 'w'], columns='tool', 
                                      values=['bacteria_mapped', 'bacteria_total', 
                                              'numt_mapped', 'numt_total', 
                                              'mito_correct', 'mito_total']).reset_index()

    # Ensure all numeric columns are integers
    for col in pivot_table.columns:
        if col[0] in ['bacteria_mapped', 'bacteria_total', 'numt_mapped', 'numt_total', 'mito_correct', 'mito_total']:
            pivot_table[col] = pivot_table[col].fillna(0).astype(int)

    # Format the columns to "mapped/total" or "correct/total"
    for metric in ['bacteria_mapped', 'numt_mapped', 'mito_correct']:
        for tool in ['SAFARI', 'vg giraffe']:
            total_col = metric.replace('mapped', 'total').replace('correct', 'total')
            pivot_table[(metric, tool)] = pivot_table[(metric, tool)].astype(str) + '/' + pivot_table[(total_col, tool)].astype(str)

    # Define the columns order with updated formatting
    columns_order = [
        ('k', ''), ('w', ''),
        ('bacteria_mapped', 'SAFARI'), ('bacteria_mapped', 'vg giraffe'),
        ('numt_mapped', 'SAFARI'), ('numt_mapped', 'vg giraffe'),
        ('mito_correct', 'SAFARI'), ('mito_correct', 'vg giraffe')
    ]
    pivot_table = pivot_table.reindex(columns=pd.MultiIndex.from_tuples(columns_order))

    # Convert to LaTeX with custom formatting to ensure no .0 for integers
    latex_table = pivot_table.to_latex(index=False, formatters={
        col: lambda x: x.split('.')[0] if '.' in x else x for col in pivot_table.columns
    })

    # Save the LaTeX table to file
    with open(output_file_path, 'w') as f:
        f.write(latex_table)

if __name__ == "__main__":
    input_file_path = sys.argv[1]  # Update this path as needed
    output_file_path = sys.argv[2]  # The output LaTeX file name
    analyze_and_save(input_file_path, output_file_path)

