import pandas as pd
import sys

def analyze_and_save(input_file_path, output_file_path):
    # Load the data
    data = pd.read_csv(input_file_path)

    # Filter for SAFARI and vg giraffe tools
    filtered_data = data[data['tool'].isin(['SAFARI', 'vg giraffe'])]

    # Extend the aggregation to sum up unmapped counts and calculate total for mito
    agg_operations = {
        'bacteria_mapped': 'median',
        'bacteria_unmapped': 'sum',
        'numt_mapped': 'median',
        'numt_unmapped': 'sum',
        'mito_correct': 'median',
        'mito_incorrect': 'sum',
        'mito_unmapped': 'sum'
    }
    grouped = filtered_data.groupby(['k', 'w', 'tool']).agg(agg_operations).reset_index()

    # Calculate total counts
    grouped['bacteria_total'] = grouped['bacteria_mapped'] + grouped['bacteria_unmapped']
    grouped['numt_total'] = grouped['numt_mapped'] + grouped['numt_unmapped']
    grouped['mito_total'] = grouped['mito_correct'] + grouped['mito_unmapped']

    # Pivot the table for side by side comparison, including total counts
    pivot_table = grouped.pivot_table(index=['k', 'w'], columns='tool', 
                                      values=['bacteria_mapped', 'bacteria_total', 
                                              'numt_mapped', 'numt_total', 
                                              'mito_correct', 'mito_total']).reset_index()

    # Convert numeric columns that should be integers back to integers and prepare for "out of" formatting
    for col in ['bacteria_mapped', 'bacteria_total', 'numt_mapped', 'numt_total', 'mito_correct', 'mito_total']:
        for tool in ['SAFARI', 'vg giraffe']:
            pivot_table[(col, tool)] = pivot_table[(col, tool)].fillna(0).astype(int)

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

    # Convert to LaTeX with custom formatting
    latex_table = pivot_table.to_latex(index=False)

    # Save the LaTeX table to file
    with open(output_file_path, 'w') as f:
        f.write(latex_table)

if __name__ == "__main__":
    input_file_path = sys.argv[1]  # Update this path as needed
    output_file_path = sys.argv[2]  # The output LaTeX file name
    analyze_and_save(input_file_path, output_file_path)

