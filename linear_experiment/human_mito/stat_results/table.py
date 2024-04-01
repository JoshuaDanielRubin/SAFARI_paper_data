import pandas as pd
import sys

def analyze_and_save(input_file_path, output_file_path):
    # Load the data
    data = pd.read_csv(input_file_path)

    # Filter for SAFARI and vg giraffe tools
    filtered_data = data[data['tool'].isin(['SAFARI', 'vg giraffe'])]

    # Group by k, w values and tool, then calculate medians for specified metrics
    grouped = filtered_data.groupby(['k', 'w', 'tool']).median().reset_index()

    # Pivot the table for side by side comparison
    pivot_table = grouped.pivot_table(index=['k', 'w'], columns='tool', values=['bacteria_mapped', 'numt_mapped', 'mito_correct']).reset_index()

    # Convert numeric columns that should be integers back to integers
    int_cols = ['k', 'w']  # Add metric columns here if they should also be integers
    for col in int_cols:
        pivot_table[col] = pivot_table[col].astype(int)
    metrics_int_cols = ['bacteria_mapped', 'numt_mapped', 'mito_correct']
    for col in metrics_int_cols:
        for tool in ['SAFARI', 'vg giraffe']:
            pivot_table[(col, tool)] = pivot_table[(col, tool)].fillna(0).astype(int)

    # Reorder columns as specified
    columns_order = [('k', ''), ('w', ''), ('bacteria_mapped', 'SAFARI'), ('bacteria_mapped', 'vg giraffe'),
                     ('numt_mapped', 'SAFARI'), ('numt_mapped', 'vg giraffe'),
                     ('mito_correct', 'SAFARI'), ('mito_correct', 'vg giraffe')]
    pivot_table = pivot_table.reindex(columns=pd.MultiIndex.from_tuples(columns_order))

    # Convert to LaTeX with integer formatting
    latex_table = pivot_table.to_latex(index=False, formatters={
        ('k', ''): lambda x: f"{x:d}",
        ('w', ''): lambda x: f"{x:d}",
        # Add formatters for metrics columns if necessary
    })

    # Save the LaTeX table to file
    with open(output_file_path, 'w') as f:
        f.write(latex_table)

if __name__ == "__main__":
    input_file_path = sys.argv[1]  # Update this path as needed
    output_file_path = sys.argv[2]  # The output LaTeX file name
    analyze_and_save(input_file_path, output_file_path)

