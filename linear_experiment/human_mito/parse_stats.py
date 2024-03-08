import re
import csv
import sys

def parse_data(file_path):
    # Regular expressions for extracting data
    damage_pattern = re.compile(r'Damage: (\w+)')
    tool_pattern = re.compile(r'Tool: (\w+)')
    metric_pattern = re.compile(r'([a-zA-Z ]+): ([0-9.]+)')
    
    data = []

    with open(file_path, 'r') as file:
        current_damage = ''
        current_tool = ''
        for line in file:
            # Check for damage
            damage_match = damage_pattern.match(line)
            if damage_match:
                current_damage = damage_match.group(1)
                continue

            # Check for tool
            tool_match = tool_pattern.match(line)
            if tool_match:
                current_tool = tool_match.group(1)
                continue

            # Check for metrics
            metric_match = metric_pattern.match(line.strip())
            if metric_match:
                metric_name = metric_match.group(1).strip()
                metric_value = metric_match.group(2)
                data.append([current_damage, current_tool, metric_name, metric_value])
                
    return data

def write_to_csv(data, output_file_path):
    # Define CSV file headers
    headers = ['Damage', 'Tool', 'Metric', 'Value']
    
    # Write data to CSV
    with open(output_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for row in data:
            writer.writerow(row)

if __name__ == "__main__":
    input_file_path = sys.argv[1]  # Adjust this to your actual input file path
    output_file_path = sys.argv[2]    # Adjust this if you want a different output file name
    
    # Parse the input file
    data = parse_data(input_file_path)
    
    # Write the data to a CSV file
    write_to_csv(data, output_file_path)
    print("Conversion complete. Data saved to:", output_file_path)

