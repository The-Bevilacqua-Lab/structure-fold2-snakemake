#!/usr/bin/env python3

'''
Finds the 5' or 3' coverage of all transcripts in a .rtsc file.
It is highly recommended that only DMS(-) .rtsc files are used for input.

5' coverage is calculated by dividing the number of stops in the 5' most n number of nucleotides
(n is set by the user with -length) by the average number of stops for an equal sized region of the transcript.
To calculate 5' coverage, set FP as mode and set the value of length. -tp_l and -trim are not required.

3' coverage divides the number of stops in the 3' most n number of transcripts (n is again set by the user with -length)
by the number of stops in larger region of the 3' end of the transcript (set by the user with -tp_l).
It is recommended that 3' coverage is calculated following the trimming of increasing numbers of nucleotides from
the 3' end (set by the user with -trim) to assess the appropriate number of nucleotides to trim for all downstream analysis.
Any transcript that is shorter than the sum of the trim and tp_l values is returned as NA.
To calculate 3' coverage, set TP as mode and set the length, tp_l and trim values.

Idea and some code by Joseph Waldron

Python 3 port. Only the tokens a calculation needs are converted to integers (3' coverage reads just the tail of each
transcript), and -sweep N computes the mean 3' coverage for every trim from 0 to N in a single pass over the <.rtsc>.
Values are written with Python 3's shortest round-trip float repr (Python 2's str() kept only 12 significant digits).
As in the original, -trim 0 slices stops[:-0], which is empty, so every eligible transcript gets 0.0.
'''

#Imports
import argparse

#Functions
def read_in_rtsc(rtsc_fyle):
    '''Streams a <.rtsc> file, yielding (transcript_name,tab-separated stops string) -- each record is a
    transcript line, a counts line and a blank line.'''
    with open(rtsc_fyle,'r') as f:
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line)
        for transcript,stops in zip(lines,lines):
            yield transcript,stops

def average(alist):
    '''Empty Docstring'''
    return sum(alist)/float(len(alist))

def check_extension(astr,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astr if astr.endswith(extension) else astr+extension
    return out_string

def tail_tokens(stops,count):
    '''The last count tab-separated tokens of stops (all of them if there are fewer), without splitting the rest'''
    parts = stops.rsplit('\t',count)
    return parts[-count:] if len(parts) > count else parts

def three_prime_value(stops,parameters):
    '''3' coverage of one transcript: a float, or 'NA' if the transcript is too short'''
    trim,tp_l,length = parameters['trim'],parameters['tp_l'],parameters['length']
    n = stops.count('\t')+1
    if n + -trim + -tp_l > 0:
        if trim >= 1 and tp_l >= 1 and length >= 1:
            #Only the tail of the transcript takes part, so only the tail is converted
            end = n - trim
            k_part = min(length,end)
            m = trim + max(tp_l,k_part)
            tail = [int(x) for x in tail_tokens(stops,m)]
            e_t = len(tail) - trim
            try:
                return (sum(tail[e_t-k_part:e_t])/float(k_part))/(sum(tail[e_t-tp_l:e_t])/float(tp_l))
            except ZeroDivisionError:
                return 0.0
        #Unusual parameters (including -trim 0, where stops[:-0] is empty): use the slicing exactly as written
        try:
            trimmed_stops = [int(x) for x in stops.split('\t')][:-trim]
            full_end = trimmed_stops[-tp_l:]
            partial_end = trimmed_stops[-length:]
            return average(partial_end)/average(full_end)
        except ZeroDivisionError:
            return 0.0
    return 'NA'

def five_prime_coverage(rtsc_data,parameters):
    '''Calculates FP coverage by Joseph's reckoning'''
    new_data,length = {},parameters['length']
    for name, stops in rtsc_data:
        tokens = stops.split('\t')
        #Counts are mostly 0, so only the non-zero tokens are converted
        total = sum(map(int,[x for x in tokens if x != '0']))
        try:
            new_data[name] = sum(map(int,tokens[:length]))/((total/float(len(tokens)))*length)
        except ZeroDivisionError:
            new_data[name] = 0.0
    extra_params = [str(parameters[q]) for q in ['length']]
    out_name = '_'.join([parameters['rtsc'].replace('.rtsc',''),'FP']+extra_params)
    out_name = check_extension(out_name,'.csv')
    return out_name,new_data

def three_prime_coverage(rtsc_data,parameters):
    '''Calculates TP coverage by Joseph's reckoning'''
    new_data = {}
    for name, stops in rtsc_data:
        new_data[name] = three_prime_value(stops,parameters)
    extra_params = [str(parameters[q]) for q in ['length','tp_l','trim']]
    out_name = '_'.join([parameters['rtsc'].replace('.rtsc',''),'TP']+extra_params)
    out_name = check_extension(out_name,'.csv')
    return out_name, new_data

def sweep_three_prime(rtsc_fyle,max_trim,parameters):
    '''Mean 3' coverage over the non-NA transcripts for every trim from 0 to max_trim, in one pass over the file.
    Returns a list of (trim,mean or None if every transcript is NA).'''
    import numpy as np
    tp_l,length = parameters['tp_l'],parameters['length']
    trims = list(range(max_trim+1))
    if tp_l < 1 or length < 1:
        #Unusual parameters: evaluate each trim with the generic code
        sums,counts = [0.0]*len(trims),[0]*len(trims)
        for name,stops in read_in_rtsc(rtsc_fyle):
            for i,trim in enumerate(trims):
                value = three_prime_value(stops,dict(parameters,trim=trim))
                if value != 'NA':
                    sums[i]+=value
                    counts[i]+=1
        return [(t,(s/c if c else None)) for t,s,c in zip(trims,sums,counts)]

    t_arr = np.arange(1,max_trim+1)
    sums,counts = np.zeros(len(trims)),np.zeros(len(trims),dtype=np.int64)
    n_tail = max_trim + max(tp_l,length)
    for name,stops in read_in_rtsc(rtsc_fyle):
        n = stops.count('\t')+1
        if n <= tp_l:
            continue #NA for every trim
        counts[0]+=1 #trim 0 always gives 0.0 for an eligible transcript (empty slice)
        valid = t_arr[n - t_arr - tp_l > 0]
        if valid.size == 0:
            continue
        tail = np.array(tail_tokens(stops,min(n_tail,n)),dtype=np.int64)
        cs = np.concatenate(([0],np.cumsum(tail)))
        end_in_tail = len(tail) - valid
        k_part = np.minimum(length,n - valid)
        part = (cs[end_in_tail] - cs[end_in_tail - k_part])/k_part.astype(float)
        full = (cs[end_in_tail] - cs[end_in_tail - tp_l])/float(tp_l)
        with np.errstate(divide='ignore',invalid='ignore'):
            values = np.where(full == 0,0.0,part/full)
        sums[valid] += values
        counts[valid] += 1
    return [(t,(sums[t]/counts[t] if counts[t] else None)) for t in trims]

def write_sweep(results,out_file):
    '''Writes the sweep as a tab-separated table, the mean formatted as awk's printf "%.6f"'''
    with open(out_file,'w') as g:
        g.write('trim\tmean_tp_coverage\n')
        for trim,mean in results:
            g.write('%d\t%s\n' % (trim,'NA' if mean is None else '%.6f' % mean))

def write_coverage(coverage_dict,out_file,mode):
    '''Writes out in csv format'''
    header = ','.join(['transcript',mode+'_coverage'])
    #Python 2 ordered strings ('NA') above every number; keep that (highest first)
    order = lambda x: (1,0.0) if x[1] == 'NA' else (0,x[1])
    with open(out_file,'w') as g:
        g.write(header+'\n')
        g.write(''.join([transcript+','+str(value)+'\n' for transcript, value in sorted(coverage_dict.items(), key=order,reverse=True)]))

def main():
    parser = argparse.ArgumentParser(description='Creates a <.csv> of 5\' or 3\' end coverage from an <.rtsc> file')
    parser.add_argument('rtsc',help='file to operate on')
    parser.add_argument('mode',type=str.upper,help='5\' or 3\' Prime',choices = ['FP','TP'])
    parser.add_argument('-length',type=int,default=50,help='[5\'& 3\', default=50] Number of bases from the end')
    parser.add_argument('-tp_l',type=int,default=300,help='[3\', default=300] Number of comparative bases')
    parser.add_argument('-trim',type=int,default=30,help='[3\', default=30] Number of bases to ignore')
    parser.add_argument('-name',type=str,default=None,help='outfile name, overrides auto-gen name')
    parser.add_argument('-sweep',type=int,default=None,help='[TP only] Instead of one <.csv>, write a table of the mean 3\' coverage for every trim from 0 to this value (-name is the table path)')
    args = parser.parse_args()
    functs = {'FP':five_prime_coverage,'TP':three_prime_coverage}

    #Sweep every trim value in a single pass
    if args.sweep is not None:
        if args.mode != 'TP':
            parser.error('-sweep is only available in TP mode')
        out_file = args.name if args.name != None else args.rtsc.replace('.rtsc','')+'_TP_sweep.tsv'
        write_sweep(sweep_three_prime(args.rtsc,args.sweep,vars(args)),out_file)
        return
    
    #Read in an rtsc and calculate coverage on the given end.
    out_file, coverage = functs[args.mode](read_in_rtsc(args.rtsc),vars(args))
    
    #Choose final out name
    out_file = out_file if args.name == None else check_extension(args.name,'.csv')
    
    #Write Out
    write_coverage(coverage,out_file,args.mode)

if __name__ == '__main__':
    main()
