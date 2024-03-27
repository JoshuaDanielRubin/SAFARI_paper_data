#!/bin/bash

# Function to process each .prof file
process_file() {
    local file="$1"
    echo "$file"

    # Skip processing for files that are already split
    if [[ $file =~ _[35]\.prof$ ]]; then
        echo "Skip processing for already split file: ${file}"
        return
    fi

    # Extract the base filename without the .prof extension
    local base_filename=$(basename "$file" .prof)
    local dir=$(dirname "$file")
    local base_path="${dir}/${base_filename}"

    # Check if the split files already exist to skip processing
    if [[ -f "${base_path}_3.prof" && -f "${base_path}_5.prof" ]]; then
        echo "Skip processing for ${file} as split files already exist."
        return
    fi

    # Split the file into two parts based on the repeated header line
    awk -v base="${base_path}" 'NR==1{print > base"_3.prof"; print > base"_5.prof"} NR>1 && NR<=6{print > base"_3.prof"; next} NR>7{print > base"_5.prof"}' "$file"
}

export -f process_file

# Use GNU Parallel to find and process all .prof files in parallel
find . -type f -name '*.prof' | parallel -j 40 process_file

