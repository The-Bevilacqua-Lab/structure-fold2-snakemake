#!/usr/bin/env python3

'''
Reformats <.rtsc> files such that correlation may be easily caluclated, either genome-wide or on a per transcript basis.
Output is a <.csv>.
You may filter output by using a coverage overlap. where you only want to get the correlation between transcripts mutually above a certain threshold.
You may filter with a nucleotide specificity filter, but this requires use of the <.fasta> mapped against via Bowtie.

Python 3 port. Rows are assembled with C-level iterators (zip/compress/join) straight from the tab-separated stop text,
and no Biopython is needed. Transcripts are written in the order they first appear in the input files, where Python 2
used arbitrary dict order.
'''

#Imports
import argparse
from itertools import compress, count, repeat
from sf3libs import read_fasta, read_rtsc_records, read_restrict, specificity_mask

#Functions
def generate_repeatabilty_dictionary(rtsc_fyles,restrict_dict=None):
    '''Read in as many <.rtsc> as you need to. Returns transcript:{file label:stops text}, keeping only restricted transcripts'''
    all_data = {}
    for rtsc_fyle in rtsc_fyles:
        label = rtsc_fyle.replace('.rtsc','')
        for transcript, stops in read_rtsc_records(rtsc_fyle,restrict_dict):
            all_data.setdefault(transcript, {})[label] = stops
    return all_data

def write_out_repeatability(stops,in_files,outfile='OUT.csv'):
    '''Writes out a nested RTSC dictionary'''
    with open(outfile,'w',buffering=1 << 20) as g:
        header = ','.join(['transcript','position']+sorted([x.replace('.rtsc','') for x in in_files]))+'\n'
        g.write(header)
        for transcript, subdict in stops.items():
            columns = [subdict[key].split('\t')[1:] for key in sorted(subdict.keys())]
            rows = '\n'.join(map(','.join,zip(repeat(transcript),map(str,count(1)),*columns)))
            if rows:
                g.write(rows+'\n')

def write_out_repeatability_spec(stops,in_files,seqs,specificity,outfile='OUT.csv'):
    '''Writes out a nested RTSC dictionary with a given specificity'''
    with open(outfile,'w',buffering=1 << 20) as g:
        header = ','.join(['transcript','position','base']+sorted([x.replace('.rtsc','') for x in in_files]))+'\n'
        g.write(header)
        for transcript, subdict in stops.items():
            base_seq = seqs[transcript][:-1]
            mask = specificity_mask(base_seq,specificity)
            columns = [compress(subdict[key].split('\t')[1:],mask) for key in sorted(subdict.keys())]
            rows = '\n'.join(map(','.join,zip(repeat(transcript),map(str,compress(count(1),mask)),compress(base_seq,mask),*columns)))
            if rows:
                g.write(rows+'\n')

def main():
    parser = argparse.ArgumentParser(description='Reformats <.rtsc> for easy correlation analysis')
    parser.add_argument('rtsc',help='Input <.rtsc> files', nargs='+')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    parser.add_argument('-name',default=None, help='Specify output file name')
    parser.add_argument('-fasta',default=None, help='<.fasta> to apply specificity')
    parser.add_argument('-spec',default=None, help='[ACGT] Nucleotide Specifictiy')
    parser.add_argument('-restrict',default=None, help='Filter to these transcripts via coverage file')
    args = parser.parse_args()
    
    #Outfile nomenclature
    default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'_correlation.csv'
    if args.spec != None:
        new_suffix = '_'+''.join(sorted(list(args.spec)))+'spec_correlation.csv'
        default_name = default_name.replace('_correlation.csv',new_suffix)
    out_name = default_name if args.name == None else args.name
    
    #Restrictions, if they exist
    restrict_dict = None if args.restrict == None else read_restrict(args.restrict)
    
    #Read in data, leaving out transcripts that are filtered out
    data = generate_repeatabilty_dictionary(args.rtsc,restrict_dict)

    #No specificity path
    if args.fasta == None and args.spec == None:
        write_out_repeatability(data,args.rtsc,out_name)
    
    #Specificity path
    elif args.fasta != None and args.spec != None:
        sequences = read_fasta(args.fasta)
        write_out_repeatability_spec(data,args.rtsc,sequences,args.spec,out_name)
    
    else:
        print('Invalid command combination.')
        print('-spec and -fasta must be invoked together')
        print('')


if __name__ == '__main__': 
    main()
