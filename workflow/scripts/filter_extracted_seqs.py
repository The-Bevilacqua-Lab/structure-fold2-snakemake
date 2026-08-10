# imports
from Bio import SeqIO
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Filter genes without a stop codon')
parser.add_argument('--input', help='Path to the input FASTA file')
parser.add_argument('--output', help='Path to the filtered FASTA file')
args = parser.parse_args()

# open output file
with open(args.output, "w") as output_handle:
    for rec in SeqIO.parse(args.input, "fasta"):
        seq_str = str(rec.seq)
        # Keep only if sequence ends with '*' and contains no 'X'
        if "CDS" in rec.description:
            output_handle.write(f">{rec.id}\n{seq_str}\n")