#########################################################################
# Converts a reactivity file to a csv file
#########################################################################

import argparse
from itertools import islice
from Bio import SeqIO
from tqdm import tqdm

def read_derived_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile, 'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript, reactivities = [n.strip() for n in next_n_lines]
            reactivities = [float(x) if x != 'NA' else 'NA' for x in reactivities.split()]
            information[transcript] = reactivities
    return information

def read_in_fasta(fasta_file):
    '''Reads in a fasta file to a dictionary'''
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, 'fasta')}

# Parse Arguments
parser = argparse.ArgumentParser(description='Converts a reactivity file to a csv file')
parser.add_argument('--react', type=str, help='<.react> file to convert')
parser.add_argument('--fasta', type=str, help='<.fasta> file used to generate the <.react>')
parser.add_argument('--output', type=str, help='Output file name')
args = parser.parse_args()

sequences = read_in_fasta(args.fasta)
reactivities = read_derived_reactivities(args.react)

with open(args.output, 'w') as fn:
    fn.write("transcript,position,base,reactivity\n")
    for transcript_id in tqdm(reactivities, desc="Processing transcripts"):
        seq = sequences.get(transcript_id)
        reacts = reactivities[transcript_id]
        if seq is None:
            print(f'Warning: {transcript_id} not found in FASTA, skipping.')
            continue
        if len(seq) != len(reacts):
            print(f'Warning: length mismatch for {transcript_id}, skipping.')
            continue
        for position, (base, value) in enumerate(zip(seq, reacts), start=1):
            fn.write(f"{transcript_id},{position},{base},{value}\n")