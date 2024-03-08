#!/usr/bin/python
import pysam
import sys

def intersects(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

def is_numt(query_name):
    return ':' in query_name and any(chromosome in query_name.split(':')[0] for chromosome in ['1', 'X', 'Y'])

def main(bam_file_path):
    bamInputFile = pysam.Samfile(bam_file_path, "rb")

    # Initialize counters for all categories
    correct_mappings, incorrect_mappings = 0, 0
    mito_mapped, mito_unmapped = 0, 0  # Added counters for mitochondrial reads
    bacteria_mapped, bacteria_unmapped, numts_mapped, numts_unmapped = 0, 0, 0, 0

    for read in bamInputFile:
        is_from_mt = read.query_name.startswith("generation")
        is_bacteria = read.query_name.startswith("NC_")
        is_numt_read = is_numt(read.query_name)

        if is_from_mt:
            if read.is_unmapped:
                mito_unmapped += 1  # Count unmapped mitochondrial reads
            else:
                fields = read.query_name.split(":")
                if len(fields) >= 4:
                    expected_start, expected_end = int(fields[2]), int(fields[3])
                    actual_start, actual_end = read.reference_start - 50, read.reference_start + read.reference_length + 50
                    if intersects(expected_start, expected_end, actual_start, actual_end):
                        correct_mappings += 1
                        mito_mapped += 1  # Count correctly mapped mitochondrial reads
                    else:
                        incorrect_mappings += 1
                        # Optionally, consider incorrectly mapped as a separate category
                else:
                    incorrect_mappings += 1
        elif is_bacteria:
            if read.is_unmapped:
                bacteria_unmapped += 1
            else:
                #print(str(read.mapping_quality))
                bacteria_mapped += 1
        elif is_numt_read:
            if read.is_unmapped:
                numts_unmapped += 1
            else:
                numts_mapped += 1

    bamInputFile.close()

    # Output comprehensive stats, including for mitochondrial reads
    print(f"Correct mitochondrial mappings: {correct_mappings}")
    print(f"Incorrect mitochondrial mappings: {incorrect_mappings}")
    print(f"Mito reads mapped: {mito_mapped}")
    print(f"Mito reads unmapped: {mito_unmapped}\n")
    print("Bacteria reads mapped:", bacteria_mapped)
    print("Bacteria reads unmapped:", bacteria_unmapped)
    print("NuMT reads mapped:", numts_mapped)
    print("NuMT reads unmapped:", numts_unmapped)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <bam_file_path>")
    else:
        main(sys.argv[1])

