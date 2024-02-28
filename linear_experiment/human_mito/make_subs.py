import pysam
import os
from Bio import SeqIO

# Function to compute substitution matrices for 5' and 3' ends
def compute_end_substitution_matrices(bam_file, reference_genome, reference_name, threshold=5):
    # Initialize matrices for 5' and 3' ends
    substitution_matrix_5p = {base: {sub_base: 0 for sub_base in "ACGT"} for base in "ACGT"}
    substitution_matrix_3p = {base: {sub_base: 0 for sub_base in "ACGT"} for base in "ACGT"}
    
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(reference_genome), 'fasta'))
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped, secondary, or supplementary alignments

            seq = read.query_sequence
            ref_seq = fasta_sequences[reference_name].seq  # Use the correct reference name
            read_length = len(seq)

            for qpos, refpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or refpos is None:
                    continue

                ref_base = ref_seq[refpos].upper()
                read_base = seq[qpos].upper()

                if ref_base in "ACGT" and read_base in "ACGT":
                    # Check if position is within threshold bases of the 5' end
                    if qpos < threshold:
                        substitution_matrix_5p[ref_base][read_base] += 1
                    # Check if position is within threshold bases of the 3' end
                    elif read_length - qpos <= threshold:
                        substitution_matrix_3p[ref_base][read_base] += 1

    return substitution_matrix_5p, substitution_matrix_3p

# Function to sort and index a BAM file
def sort_and_index_bam(original_bam_path):
    sorted_bam_path = original_bam_path.replace(".bam", "_sorted.bam")
    pysam.sort("-o", sorted_bam_path, original_bam_path)
    pysam.index(sorted_bam_path)
    return sorted_bam_path

# Main processing function
def process_bam_files(alignments_dir, output_dir, reference_genome_base):
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    for bam_file in os.listdir(alignments_dir):
        if bam_file.endswith(".bam"):
            print(f"Processing {bam_file}...")

            # Choose reference genome and name based on file name
            if 'safari' in bam_file or 'giraffe' in bam_file:
                reference_genome = os.path.join(reference_genome_base, "H2a2a1.fa")
                reference_name = "generation_0"
            else:
                reference_genome = os.path.join(reference_genome_base, "H2a2a1.fa")
                reference_name = "generation_0"

            full_bam_path = os.path.join(alignments_dir, bam_file)
            sorted_bam_path = sort_and_index_bam(full_bam_path)

            # Compute substitution matrices for 5' and 3' ends
            matrix_5p, matrix_3p = compute_end_substitution_matrices(sorted_bam_path, reference_genome, reference_name)

            # Save the substitution matrices to files
            for end, matrix in [("5p", matrix_5p), ("3p", matrix_3p)]:
                output_file_path = os.path.join(output_dir, f"{bam_file.replace('.bam', '')}_{end}_substitution_matrix.txt")
                with open(output_file_path, 'w') as f:
                    for base in "ACGT":
                        f.write("\t".join([base] + [str(matrix[base][sub_base]) for sub_base in "ACGT"]) + "\n")

            # Optional: Remove sorted and indexed BAM files if not needed
            # os.remove(sorted_bam_path)
            # os.remove(sorted_bam_path + ".bai")

    print("Processing complete.")

# Example usage
alignments_dir = "alignments/"
output_dir = os.path.join(alignments_dir, "profs/")
reference_genome_base = "."

process_bam_files(alignments_dir, output_dir, reference_genome_base)

