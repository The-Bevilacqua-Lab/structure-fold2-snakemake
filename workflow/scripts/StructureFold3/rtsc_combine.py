#!/usr/bin/env python3

'''
This script combines any numbers of <.rtsc> files.
It is imperative that all the files used to generate the <.rtsc> used an identical reference in generation.

Python 3 port. Counts are overwhelmingly zero, so each transcript is merged as its short list of non-zero positions.
Transcripts are written in first-seen order (files taken in sorted order), where Python 2 used arbitrary dict order.
'''

#Imports
import argparse
from sf3libs import read_rtsc_records, nonzero_counts, check_extension

#Functions
def merge_rtsc(rtsc_lyst):
    '''Takes a list of <.rtsc> files, returns a dictionary transcript:(length,{position:stops}) that is the sum of the RT stops.
    As in the original, merging two vectors of different lengths keeps only the shorter length.'''
    master_dictionary = {}
    for rtsc in sorted(rtsc_lyst):
        data = {}
        for transcript,stops in read_rtsc_records(rtsc):
            length,positions,values = nonzero_counts(stops)
            data[transcript] = (length,dict(zip(positions,values)))
        for transcript,(length,counts) in data.items():
            if transcript in master_dictionary:
                current_length,current = master_dictionary[transcript]
                shortest = min(current_length,length)
                merged = {}
                for position in current.keys() | counts.keys():
                    if position < shortest:
                        merged[position] = current.get(position,0) + counts.get(position,0)
                master_dictionary[transcript] = (shortest,merged)
            else:
                master_dictionary[transcript] = (length,counts)
    return master_dictionary

def write_rtsc(rtsc_dictionary,sort_flag=False,outfile='data.rtsc'):
    '''Takes a dictionary, writes to a file'''
    items = rtsc_dictionary.items() if sort_flag == False else sorted(rtsc_dictionary.items())
    with open(outfile,'w',buffering=1 << 20) as g:
        for transcript,(length,counts) in items:
            values = ['0']*length
            for position,number in counts.items():
                values[position] = str(number)
            g.write(transcript+'\n'+'\t'.join(values)+'\n\n')

def main():
    parser = argparse.ArgumentParser(description='Combines <.rtsc> files, typically replicates of the same sample')
    parser.add_argument('rtsc',help='Input <.rtsc> files', nargs='+')
    parser.add_argument('-sort',action='store_true',default=False,help = 'Sort output by transcript name')
    parser.add_argument('-name',default=None, help='Specify output file name')
    args = parser.parse_args()
    #Generate name or assign the user provided name
    default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'.rtsc'
    out_name = default_name if args.name == None else check_extension(args.name,'.rtsc')
    #Pool all <.rtsc> into a dictionary
    all_stops = merge_rtsc(args.rtsc)
    #Write out the dictionary
    write_rtsc(all_stops,args.sort,out_name)


if __name__ == '__main__': 
    main()
