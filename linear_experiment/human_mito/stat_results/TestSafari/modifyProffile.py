import sys

def multiply_and_output(file_path, factor):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line_index, line in enumerate(lines):
        values = line.split()
        if not values:  # Skip empty lines
            continue
        if line_index == 0:
            # Print header as is
            print("\t".join(values))
        else:
            # Multiply numeric values by factor
            modified_values = [str(float(value) * factor) for value in values]
            print("\t".join(modified_values))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3.8 modify.py [input prof] [factor]")
        sys.exit(1)

    input_file = sys.argv[1]
    try:
        multiplication_factor = float(sys.argv[2])
    except ValueError:
        print("Error: Factor must be a number.")
        sys.exit(1)

    multiply_and_output(input_file, multiplication_factor)
