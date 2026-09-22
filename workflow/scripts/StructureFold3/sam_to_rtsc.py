#!/usr/bin/env python3

'''
Converts filtered <.sam> files into <.rtsc> files, requires the original <.fasta> containing the transcripts which 
were mapped against to generate the <.sam>. This script will not work on a <.sam> with a sam header; this is 
intentional as to ensure users have done the filtering step via sam_batch_filter.py script before generating <.rtsc>.
Be aware that batch mode will attempt to use all <.sam> in the directory, so make sure to have either moved the intermediate
non-filtered <.sam> to another directory before using this mode, or use the -suffix option to only operate on filtered files.

Python 3 port. Needs no Biopython (sf3libs reads the fasta). Transcripts are written in fasta order.
'''

#Imports
import glob
import argparse
from collections import Counter
from sf3libs import read_fasta

#Functions
def sam_to_rtsc_dict(samfyle):
    '''Reads in the RT stops from a <.sam> file, returns a dictionary with nested Counter objects logging the stops.'''
    stops = {}
    with open(samfyle,'r') as f:
        for line in f:
            #Only the reference name and position (fields 3 and 4) are needed
            bits = line.strip().split('\t',4)
            gene,position = bits[2],bits[3]
            stops.setdefault(gene, Counter())[int(position)-1]+=1
    return stops

def write_rtsc(stop_counts,transcript_limits,outfyle='stops.rtsc'):
    '''Writes the stops up to the limit of the transcript length, as an <.rtsc>'''
    with open(outfyle,'w') as g:
        for transcript, limit in transcript_limits.items():
            out_values = ['0']*limit
            if transcript in stop_counts:
                for index,number in stop_counts[transcript].items():
                    if 0 <= index < limit:
                        out_values[index] = str(number)
            g.write(transcript+'\n')
            g.write('\t'.join(out_values)+'\n'+'\n')

def main():
    parser = argparse.ArgumentParser(description='Creates <.rtsc> file(s) from filtered <.sam> file(s) and a <.fasta> file')
    parser.add_argument('index',type=str,help='<.fasta> file containing the transcripts mapped against')
    parser.add_argument('-single',default=None,help='Operate on this single file, rather than the directory')
    parser.add_argument('-suffix',default=None,help='Operate only on <.sam> with this suffix before the extension')
    parser.add_argument('-trim',default=None,help='Remove this suffix from output file name before writing')
    args = parser.parse_args()
    
    #Generate transcript lengths
    fasta_limits = dict([(name,len(seq)) for name,seq in read_fasta(args.index).items()])

    #Batch mode
    if args.single == None:
        sam_files = glob.glob('*.sam') if args.suffix == None else glob.glob(''.join(['*',args.suffix,'.sam']))
        for sam_file in sorted(sam_files):
            stops = sam_to_rtsc_dict(sam_file)
            outfyle = sam_file.replace('.sam','.rtsc') if args.trim == None else sam_file.replace(args.trim+'.sam','')+'.rtsc'
            write_rtsc(stops,fasta_limits,outfyle)

    #Single Mode
    else:
        stops = sam_to_rtsc_dict(args.single)
        outfyle = args.single.replace('.sam','.rtsc') if args.trim == None else args.single.replace(args.trim+'.sam','')+'.rtsc'
        write_rtsc(stops,fasta_limits,outfyle)


if __name__ == '__main__':
    main()
