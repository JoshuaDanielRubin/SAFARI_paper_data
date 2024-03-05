import pandas as pd
import os

# Specify the folder path containing the CSV files
folder_path = '.'

# List all CSV files in the folder
csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

# Initialize an empty DataFrame to store combined data
combined_df = pd.DataFrame()

# Iterate over the list of CSV files and append them to the combined DataFrame
# An additional column 'Source File' is added to identify the source of each row
for file in csv_files:
    file_path = os.path.join(folder_path, file)
    temp_df = pd.read_csv(file_path)
    temp_df['Source File'] = file  # Add source file identifier
    combined_df = pd.concat([combined_df, temp_df], ignore_index=True)

# Save the combined DataFrame to a new CSV file
output_file_path = os.path.join(folder_path, 'combined_data.csv')
combined_df.to_csv(output_file_path, index=False)

# Print the path to the combined CSV file
print(output_file_path)

