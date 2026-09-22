#!/usr/bin/env python3

'''
Creates a <.csv> of stop coverages (stops per specific base) from <.rtsc> files.

Python 3 port of StructureFold2's rtsc_coverage.py. Does not need Biopython or sf2libs, and streams each
<.rtsc> one transcript at a time instead of loading every stop count into memory. Coverage values are written
rounded to -digits decimal places (default 5). Rounding applies to the written <.csv> only; the -ol overlap
file is still computed from the unrounded coverages.
'''

#Imports
import argparse
from itertools import compress, islice
from multiprocessing import Pool

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    return astring if astring.endswith(extension) else astring + extension

def read_fasta(afasta):
    '''Fasta to Python dictionary, transcript_name:sequence (name is the first word of the header)'''
    sequences,name,parts = {},None,[]
    with open(afasta,'r') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = ''.join(parts)
                header = line[1:].split(None,1)
                name,parts = (header[0] if header else ''),[]
            elif name is not None:
                parts.append(line.strip().replace(' ',''))
    if name is not None:
        sequences[name] = ''.join(parts)
    return sequences

def read_rtsc(rtsc_fyle):
    '''Streams a <.rtsc> file, yielding (transcript_name,[stop tokens]) -- each record is a transcript line,
    a tab-separated counts line and a blank line.'''
    with open(rtsc_fyle,'r') as f:
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line)
        for transcript,stops in zip(lines,lines):
            yield transcript,stops.split('\t')

def rtsc_coverage(rtsc_file,fasta_index,specificity='AC'):
    '''Generates a coverage dictionary for an <.rtsc> file with a given specificity'''
    #Maps specific bases to 1 and everything else to 0
    table = bytes(1 if chr(i) in specificity else 0 for i in range(256))
    coverage = {}
    for transcript,stops in read_rtsc(rtsc_file):
        effective_sequence = fasta_index[transcript].upper()[:-1]
        mask = effective_sequence.encode('ascii','replace').translate(table)
        specific_bases_seq = mask.count(1)
        #compress stops at the shorter of the two, like zip; zero counts need no conversion
        selected = compress(islice(stops,1,None),mask)
        specific_bases_stops = sum(map(int,[x for x in selected if x != '0']))
        try:
            coverage[transcript] = float(specific_bases_stops)/specific_bases_seq
        except ZeroDivisionError:
            coverage[transcript] = 0
    return coverage

_FASTA_INDEX = None

def _coverage_job(args):
    fyle,specificity = args
    return rtsc_coverage(fyle,_FASTA_INDEX,specificity)

def collect_coverages(fyle_list,fasta_fyle,specificity='AC',threads=1):
    '''Applies rtsc_coverage to files'''
    global _FASTA_INDEX
    _FASTA_INDEX = read_fasta(fasta_fyle)
    jobs = [(fyle,specificity) for fyle in fyle_list]
    if threads > 1 and len(jobs) > 1:
        with Pool(min(threads,len(jobs))) as pool: #forked workers share the fasta index
            results = pool.map(_coverage_job,jobs)
    else:
        results = [_coverage_job(job) for job in jobs]
    return {fyle.replace('.rtsc',''):result for fyle,result in zip(fyle_list,results)}

def write_coverage(data,out_fyle='derp.csv',digits=5):
    '''Writes out coverages, rounded to the given number of decimal places'''
    h_keys = sorted(data.keys())
    v_keys = sorted(list(set.union(*map(set, data.values()))))
    header = ','.join(['transcript']+[derp+'_coverage' for derp in h_keys])
    with open(out_fyle,'w') as g:
        g.write(header+'\n')
        for transcript in v_keys:
            entry = [format(data[h_key][transcript],'.%df' % digits) if transcript in data[h_key] else 'NA' for h_key in h_keys]
            out_line = ','.join([transcript]+entry)
            g.write(out_line+'\n')

def write_ol(data,out_fyle='derp.txt',threshold=1.0):
    '''Writes out an overlap file'''
    h_keys = sorted(data.keys())
    shared_transcripts = sorted(list(set.intersection(*map(set, data.values()))))
    with open(out_fyle,'w') as g:
        for transcript in shared_transcripts:
            passing = all([data[h_key][transcript] >= threshold for h_key in h_keys])
            if passing:
                g.write(transcript+'\n')

def main():
    parser = argparse.ArgumentParser(description='Creates a <.csv> of stop coverages from <.rtsc> files')
    parser.add_argument('-f',type=str,help='<.rtsc> files to operate on', nargs='+')
    parser.add_argument('index',type=str,help='<.fasta> file used to generate <.rtsc>')
    parser.add_argument('-bases',type=str,default='AC', help='[default = AC] Coverage Specificity')
    parser.add_argument('-name',type=str,default = None, help='Output file name')
    parser.add_argument('-ol',action='store_true', help='Create an overlap file')
    parser.add_argument('-ot',type=float,default=1.0, help='[default = 1.0] Overlap file threshold')
    parser.add_argument('-on',type=str,default=None, help='Overlap file name')
    parser.add_argument('-digits',type=int,default=5, help='[default = 5] Decimal places written for coverage values')
    parser.add_argument('-threads',type=int,default=1, help='[default = 1] Number of <.rtsc> files processed in parallel')
    args = parser.parse_args()
    
    #Generate or assign name
    default_name = '_'.join(sorted([fyle.replace('.rtsc','') for fyle in args.f])+['coverage'])+'.csv'
    out_name = default_name if args.name == None else check_extension(args.name,'.csv')
    
    #Collect Data
    coverage_data = collect_coverages(args.f,args.index,args.bases,args.threads)
    
    #Write Data
    write_coverage(coverage_data,out_name,args.digits)
    
    #Create overlap file
    if args.ol:
        default_ol = '_'.join(sorted([fyle.replace('.rtsc','') for fyle in args.f])+['overlap',str(args.ot)])+'.txt'
        out_ol = default_ol if args.on == None else check_extension(args.on,'.txt')
        write_ol(coverage_data,out_ol,args.ot)

if __name__ == '__main__':
    main()
