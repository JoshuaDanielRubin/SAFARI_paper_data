import pandas as pd
import sys

def analyze_and_save(input_file_path, output_file_path):
    # Load the data
    data = pd.read_csv(input_file_path)
    
    # Filter for subsampling_rate=0.9
    data = data[data['subsampling_rate'] == 0.9]
    data = data[data['damage_level'] == 'High']
    
    # Assuming 'bacteria_mapped' as 'bacteria_correct' for simplification
    data['bacteria_total'] = data['bacteria_mapped'] + data['bacteria_unmapped']
    data['numt_total'] = data['numt_mapped'] + data['numt_unmapped']
    data['mito_total'] = data['mito_correct'] + data['mito_incorrect'] + data['mito_unmapped']

    # Filter for SAFARI and vg giraffe tools
    filtered_data = data[data['tool'].isin(['SAFARI', 'vg giraffe'])]

    # Select necessary columns and convert to integer for aggregation
    columns = ['k', 'w', 'tool', 'bacteria_mapped', 'bacteria_total', 
               'numt_mapped', 'numt_total', 'mito_correct', 'mito_total']
    filtered_data = filtered_data[columns].astype({'bacteria_mapped': 'int', 'bacteria_total': 'int',
                                                   'numt_mapped': 'int', 'numt_total': 'int',
                                                   'mito_correct': 'int', 'mito_total': 'int'})
    
    # Pivot the table for side by side comparison
    pivot_table = pd.pivot_table(filtered_data, values=['bacteria_mapped', 'bacteria_total', 
                                                        'numt_mapped', 'numt_total', 
                                                        'mito_correct', 'mito_total'], 
                                 index=['k', 'w'], columns=['tool'], aggfunc='first').reset_index()

    # Formatting for "correct/total"
    for metric in ['bacteria_mapped', 'numt_mapped', 'mito_correct']:
        total_col = metric.replace('_mapped', '_total').replace('_correct', '_total')
        for tool in ['SAFARI', 'vg giraffe']:
            pivot_table[(metric, tool)] = pivot_table[(metric, tool)].astype(str) + '/' + pivot_table[(total_col, tool)].astype(str)

    # Select and rename columns appropriately
    pivot_table.columns = [' '.join(col).strip() for col in pivot_table.columns.values]
    pivot_table.rename(columns={
        'k ': 'k', 'w ': 'w',
        'bacteria_mapped SAFARI': 'Bacteria SAFARI', 'bacteria_mapped vg giraffe': 'Bacteria vg giraffe',
        'numt_mapped SAFARI': 'NuMT SAFARI', 'numt_mapped vg giraffe': 'NuMT vg giraffe',
        'mito_correct SAFARI': 'Mito SAFARI', 'mito_correct vg giraffe': 'Mito vg giraffe'}, inplace=True)

    # Simplify the table to only include necessary comparisons
    pivot_table = pivot_table[['k', 'w', 'Bacteria SAFARI', 'Bacteria vg giraffe',
                               'NuMT SAFARI', 'NuMT vg giraffe', 'Mito SAFARI', 'Mito vg giraffe']]

    # Save the table to a LaTeX file
    with open(output_file_path, 'w') as latex_file:
        latex_file.write(pivot_table.to_latex(index=False, escape=False))

if __name__ == "__main__":
    input_file_path = sys.argv[1]  # The input CSV file name
    output_file_path = sys.argv[2]  # The output LaTeX file name
    analyze_and_save(input_file_path, output_file_path)

