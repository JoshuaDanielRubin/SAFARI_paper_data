import pysam
import os
from Bio import SeqIO

# Function to compute substitution matrices for 5' and 3' ends
def compute_end_substitution_matrices(bam_file, reference_genome, reference_name, threshold=5):
    substitution_matrix_5p = {base: {sub_base: 0 for sub_base in "ACGT"} for base in "ACGT"}
    substitution_matrix_3p = {base: {sub_base: 0 for sub_base in "ACGT"} for base in "ACGT"}
    
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(reference_genome), 'fasta'))
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            seq = read.query_sequence
            ref_seq = fasta_sequences[reference_name].seq
            read_length = len(seq)

            for qpos, refpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or refpos is None:
                    continue

                ref_base = ref_seq[refpos].upper()
                read_base = seq[qpos].upper()

                if ref_base in "ACGT" and read_base in "ACGT":
                    if qpos < threshold:
                        substitution_matrix_5p[ref_base][read_base] += 1
                    elif read_length - qpos <= threshold:
                        substitution_matrix_3p[ref_base][read_base] += 1

    return substitution_matrix_5p, substitution_matrix_3p

# Function to sort and index a BAM file
def sort_and_index_bam(original_bam_path):
    sorted_bam_path = original_bam_path.replace(".bam", "_sorted.bam")
    pysam.sort("-o", sorted_bam_path, original_bam_path)
    pysam.index(sorted_bam_path)
    return sorted_bam_path

# Main processing function with modifications
def process_bam_files(alignments_dir, output_dir, reference_genome_base):
    os.makedirs(output_dir, exist_ok=True)

    for bam_file in os.listdir(alignments_dir):
        if 'sorted' in bam_file:
            continue
        if bam_file.endswith(".bam"):
            output_file_path_5p = os.path.join(output_dir, f"{bam_file.replace('.bam', '')}_5p_substitution_matrix.txt")
            output_file_path_3p = os.path.join(output_dir, f"{bam_file.replace('.bam', '')}_3p_substitution_matrix.txt")

            if os.path.exists(output_file_path_5p) and os.path.exists(output_file_path_3p):
                print(f"Skipping {bam_file} as it has already been processed.")
                print(output_file_path_5p)
                continue

            print(f"Processing {bam_file}...")

            reference_genome = os.path.join(reference_genome_base, "H2a2a1.fa")
            reference_name = "generation_0"

            full_bam_path = os.path.join(alignments_dir, bam_file)
            sorted_bam_path = sort_and_index_bam(full_bam_path)

            matrix_5p, matrix_3p = compute_end_substitution_matrices(sorted_bam_path, reference_genome, reference_name)

            # Writing after each file processed to manage memory usage
            with open(output_file_path_5p, 'w') as f5p, open(output_file_path_3p, 'w') as f3p:
                for base in "ACGT":
                    f5p.write("\t".join([base] + [str(matrix_5p[base][sub_base]) for sub_base in "ACGT"]) + "\n")
                for base in "ACGT":
                    f3p.write("\t".join([base] + [str(matrix_3p[base][sub_base]) for sub_base in "ACGT"]) + "\n")

            # Optional: Remove sorted and indexed BAM files if not needed
            # os.remove(sorted_bam_path)
            # os.remove(sorted_bam_path + ".bai")

    print("Processing complete.")

# Example usage
alignments_dir = "alignments/"
output_dir = os.path.join(alignments_dir, "subs/")
reference_genome_base = "."

process_bam_files(alignments_dir, output_dir, reference_genome_base)

