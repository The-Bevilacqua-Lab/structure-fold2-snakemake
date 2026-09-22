#!/usr/bin/env python3

'''
By default, this will compare every <.rtsc> in the directory to the -index fasta file, and generate a small report.
A custom set of files can instead be given to the program.

Python 3 port. Only the non-zero stop counts are visited, and no Biopython is needed. Fractions are written like
Python 2's str() (12 significant digits), so the <.csv> matches the Python 2 version.
'''

#Imports
import glob
import argparse
from collections import Counter
from sf3libs import read_fasta, read_rtsc_records, nonzero_counts, check_extension, py2_str

#Functions
def rtsc_to_specificity(rtsc_fyle,fasta_index):
    '''Returns a Counter object of the frequencies of the nucleotide specificities'''
    Xcounter = Counter()
    for transcript, stops in read_rtsc_records(rtsc_fyle):
        seq = fasta_index[transcript].upper()[:-1]
        length,positions,values = nonzero_counts(stops)
        #Position i of the stops (i >= 1) pairs with seq[i-1]; zip stopped at the shorter of the two
        for position,value in zip(positions,values):
            if 1 <= position <= len(seq):
                Xcounter[seq[position-1]] += value
    return Xcounter
            
def write_specificity_data(info,outfyle='specificity.csv'):
    '''Dumps to a .csv'''
    header_keys = sorted([x.replace('.rtsc','') for x in info.keys()])
    types = ['_count','_specificity']
    header = ','.join(['base']+[fyle+mod for fyle in header_keys for mod in types])
    f_keys = sorted(info.keys())
    x_data = [[base]+[(info[q][base],info[q][base]/float(sum(info[q].values()))) for q in f_keys] for base in 'ATGC']
    with open(outfyle,'w') as g:
        g.write(header+'\n')
        for data in x_data:
            g.write(','.join(data[0:1]+[str(count)+','+py2_str(fraction) for count,fraction in data[1:]])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Analyzes native/reagent nucleotide stop specificity')
    parser.add_argument('-index',type=str,help='<.fasta> file used to generate the <.rtsc>')
    parser.add_argument('-rtsc',default = None, help='Operate on specific <.rtsc>', nargs='+')
    parser.add_argument('-name',default = None, help='Specify output file name')
    args = parser.parse_args()
    
    #Read in fasta sequences
    seqs = read_fasta(args.index)

    #Batch Mode
    if args.rtsc == None:
        #Generate Out Name
        fyles,info = sorted(glob.glob('*.rtsc')),{}
        default_name = '_'.join([x.replace('.rtsc','') for x in fyles])+'_specificity.csv'
        out_name = default_name if args.name == None else check_extension(args.name,'.rtsc')
        for fyle in fyles:
            info[fyle] = rtsc_to_specificity(fyle,seqs)
        write_specificity_data(info,out_name)
    
    #Specific Mode
    if args.rtsc != None:
        #Generate Out Name
        info = {}
        default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'_specificity.csv'
        out_name = default_name if args.name == None else check_extension(args.name,'.csv')
        for fyle in args.rtsc:
            info[fyle] = rtsc_to_specificity(fyle,seqs)
        write_specificity_data(info,out_name)

if __name__ == '__main__': 
    main()
