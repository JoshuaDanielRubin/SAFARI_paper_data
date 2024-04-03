#/usr/bin/python

import sys

def fasta_to_fastq(fasta_file, quality_score, output_file):
    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        seq_id = ''
        seq = ''
        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    write_fastq_record(f_out, seq_id, seq, quality_score)
                seq_id = line[1:]
                seq = ''
            else:
                seq += line
        if seq:
            write_fastq_record(f_out, seq_id, seq, quality_score)

def write_fastq_record(f_out, seq_id, seq, quality_score):
    f_out.write('@{}\n{}\n+\n{}\n'.format(seq_id, seq, quality_score * len(seq)))

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python script.py input.fasta quality_score output.fastq')
        sys.exit(1)

    fasta_file = sys.argv[1]
    quality_score = sys.argv[2]
    output_file = sys.argv[3]

    fasta_to_fastq(fasta_file, quality_score, output_file)

    
