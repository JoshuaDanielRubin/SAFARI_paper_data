import pandas as pd
from io import StringIO

# Function to process LaTeX data into a CSV format acceptable for Pandas
def process_latex_data_for_pandas(latex_table):
    lines = latex_table.split('\n')
    data_lines = [line for line in lines if line.strip() and not line.strip().startswith('\\')]
    data_lines = [line for line in data_lines if not line.strip().startswith('end')]
    processed_lines = []
    for line in data_lines:
        cells = [cell.replace('\\', '').replace('BETTER', '').replace('WORSE', '').strip() for cell in line.split('&') if cell.strip()]
        processed_line = ','.join(cells)
        processed_lines.append(processed_line)
    csv_table = '\n'.join(processed_lines)
    return csv_table

# Example usage with a hypothetical path to your LaTeX table
file_path = 'table.txt'
with open(file_path, 'r') as file:
    full_latex_table = file.read()

csv_table_improved = process_latex_data_for_pandas(full_latex_table)
df_full_improved = pd.read_csv(StringIO(csv_table_improved))

# Convert to numeric, coerce non-numeric to NaN
df_full_improved['Reads_corrected'] = pd.to_numeric(df_full_improved['Reads_corrected'], errors='coerce')
df_full_improved['Reads_uncorrected'] = pd.to_numeric(df_full_improved['Reads_uncorrected'], errors='coerce')

# Statistical calculations for 'Reads_corrected' and 'Reads_uncorrected'
def stats(data):
    desc = data.describe()
    iqr = desc['75%'] - desc['25%']
    stats_data = {
        'mean': desc['mean'],
        'std': desc['std'],
        'median': desc['50%'],
        'lower_quartile': desc['25%'],
        'upper_quartile': desc['75%'],
        'iqr': iqr
    }
    return stats_data

stats_corrected = stats(df_full_improved['Reads_corrected'])
stats_uncorrected = stats(df_full_improved['Reads_uncorrected'])

print("Stats for Reads Corrected:", stats_corrected)
print("Stats for Reads Uncorrected:", stats_uncorrected)

