#!/usr/bin/env python3

'''
This script takes two <.rtsc> along with the base <.fasta> to calculate the reactivity scores of each nucelotide.
If not given a <.scale> file, the script will generate a <.scale> to be used with this script when calculating any reactvities you wish to do a comparison to.
Transcripts which produce a 0 average on the normalization scale will not have their final reactivity calculated, and these transcripts may be logged to a file (optional).
If given a <.scale> file, it will apply it as the normalization scale. Transcripts not included in the <.scale> will not be calculated and can be logged to a file (optional).
Normalization (2-8%) may be turned off for the reactivity calculation, as can taking the natural log of the reactivity values, via options.

Python 3 port of StructureFold2's rtsc_to_react.py, written to give the same <.react> and <.scale> content. Stop counts
are overwhelmingly zero, so each transcript is handled as the short list of its non-zero positions instead of full-length
vectors. Sums are added sequentially (Python 3.12+ sum() compensates float sums, which would change the last bits) and
<.scale> values are written like Python 2's str() (12 significant digits), so results match the Python 2 version.
Transcripts are written in file order rather than Python 2's arbitrary dict order, and the Biopython dependency is gone.
'''

#Imports
import math
import argparse
import os
import re
from functools import lru_cache, reduce
from operator import add

#Functions
def py2_str(value):
    '''str() of a float as Python 2 wrote it: 12 significant digits, and a trailing .0 on whole numbers'''
    text = '%.12g' % value
    return text + '.0' if text.lstrip('-').isdigit() else text

def sequential_sum(values):
    '''sum() that always adds strictly left to right, starting from the int 0 (like Python 2)'''
    return reduce(add,values,0)

def get_covered_transcripts(coverage_fyle):
    '''Reads in a standard overlapped coverage file'''
    info = {}
    with open(coverage_fyle,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

def read_in_fasta(fasta_fyle):
    '''Reads in a fasta file to a dictionary, transcript_name:transcript_sequence'''
    fasta_dict,name,parts = {},None,[]
    with open(fasta_fyle,'r') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    fasta_dict[name] = ''.join(parts)
                header = line[1:].split(None,1)
                name,parts = (header[0] if header else ''),[]
            elif name is not None:
                parts.append(line.strip().replace(' ',''))
    if name is not None:
        fasta_dict[name] = ''.join(parts)
    return fasta_dict

def read_in_rtsc(rtsc_fyle,keep=None):
    '''Reads a <.rtsc> file into a dictionary, transcript_name:tab-separated stops (kept as text and only
    parsed for the transcripts that are actually used). keep optionally limits which transcripts are read.'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line)
        for transcript,stops in zip(lines,lines):
            if keep is None or transcript in keep:
                information[transcript] = stops
    return information

_NONZERO_DIGIT = re.compile('[1-9]').finditer

def nonzero_counts(stops):
    '''Returns (number of positions, [positions with a non-zero count], [those counts as floats]).
    Counts are overwhelmingly 0, so rather than splitting every token this finds the non-zero digits and
    converts just the tokens containing them.'''
    positions,values = [],[]
    token_end,tabs,scanned = -1,0,0
    for match in _NONZERO_DIGIT(stops):
        digit = match.start()
        if digit < token_end:
            continue #another digit of a token already taken
        token_start = stops.rfind('\t',0,digit)+1
        token_end = stops.find('\t',digit)
        if token_end == -1:
            token_end = len(stops)
        tabs += stops.count('\t',scanned,token_start)
        scanned = token_start
        positions.append(tabs)
        values.append(float(stops[token_start:token_end]))
    return stops.count('\t')+1,positions,values

def calculate_raw_reactivity(reagent_minus,reagent_plus,nlog_off=False):
    '''Calculates raw reactivity, with or without the natural log. Each transcript is returned as
    (length,positions,values): the vector has length entries, all 0 except the listed positive values.'''
    data_out,log_cache = {},{}
    def log_plus_one(value):
        try:
            return log_cache[value]
        except KeyError:
            result = log_cache[value] = math.log(value+1,math.e)
            return result
    for key,plus_stops in reagent_plus.items():
        if key not in reagent_minus:
            continue
        length,plus_pos,plus_vals = nonzero_counts(plus_stops)
        minus_length,minus_pos,minus_vals = nonzero_counts(reagent_minus[key])
        if nlog_off == False:
            plus_vals = [log_plus_one(v) for v in plus_vals]
            minus_vals = [log_plus_one(v) for v in minus_vals]
        sum_plus,sum_minus = sequential_sum(plus_vals),sequential_sum(minus_vals)
        if sum_plus != 0 and sum_minus != 0:
            minus_nrm = {p:float(z)/float(sum_minus)*length for p,z in zip(minus_pos,minus_vals)}
            shortest = min(length,minus_length)
            positions,values = [],[]
            for p,y in zip(plus_pos,plus_vals):
                #Only positions with plus signal can be positive after subtracting the minus
                if p < shortest:
                    difference = float(y)/float(sum_plus)*length - minus_nrm.get(p,0.0)
                    if difference > 0:
                        positions.append(p)
                        values.append(difference)
            data_out[key] = (shortest,positions,values)
    return data_out

@lru_cache(maxsize=None)
def specificity_table(specificity):
    '''Translation table mapping the specificity bases to 1 and every other byte to 0'''
    return bytes(1 if chr(i) in specificity else 0 for i in range(256))

def specificity_mask(sequence,specificity,needed):
    '''bytes of 1 where sequence[i] is one of the specificity bases (else 0), for the first needed positions.
    Like indexing the sequence, running off its end is an IndexError.'''
    if needed > len(sequence):
        raise IndexError('string index out of range')
    return sequence[:needed].encode('ascii','replace').translate(specificity_table(specificity))

def generate_normalization_scale(derived_reactivities,transcript_seqs,specificity,trim3=0):
    '''Generates the 2-8% scale to normalize against'''
    data = {}
    for transcript, (length,positions,values) in derived_reactivities.items():
        sequence = transcript_seqs[transcript]
        stop = max(1, length - trim3)
        mask = specificity_mask(sequence,specificity,stop-1)
        #The accepted list is every specific position (mostly 0s); only its positive values need sorting
        accepted_count = mask.count(1)
        positive = sorted([v for p,v in zip(positions,values) if 1 <= p < stop and mask[p-1]],reverse=True)
        first,last = int(accepted_count*0.02),int(accepted_count*0.1)
        top_length = last - first
        top_average = sequential_sum(positive[first:last])/top_length if top_length > 0 else 0
        if top_average > 0:
            data[transcript] = top_average
    return data

def read_normalization_scale(normalization_file):
    '''Reads in a normalization scale file'''
    info = {}
    with open(normalization_file, 'r') as f:
        for line in f:
            if line.startswith('transcript,'):
                continue
            else:
                transcript,value = line.strip().split(',')
                info[transcript] = float(value)
    return info

def write_normalization_scale(scale_dictionary,outfile):
    '''Writes out a normalization scale file'''
    with open(outfile,'w') as g:
        g.write(','.join(['transcript','value'])+'\n')
        g.write(''.join([transcript+','+py2_str(value)+'\n' for transcript, value in scale_dictionary.items()]))

def format_reactivity(value,threshold):
    '''One reactivity as written to a <.react> file'''
    return str(float('%.3f'%min(value, threshold)))

def calculate_final_reactivity(derived_reactivities,transcript_sequences,specificity,threshold,nrm_scale,norm_off=False):
    '''Calculates the final reactivity'''
    data_out,missing_transcripts = {},{}
    for transcript, (length,positions,values) in derived_reactivities.items():
        if transcript in nrm_scale:
            normalizer = nrm_scale[transcript] if norm_off == False else 1
            sequence = transcript_sequences[transcript]
            mask = specificity_mask(sequence,specificity,length-1)
            #Positions 1..length-1 are reported (position i is the base at sequence[i-1]), then a final NA
            if mask.count(1):
                zero_text = format_reactivity(0/normalizer,threshold)
                normalized_values = [zero_text if m else 'NA' for m in mask]
                for p,v in zip(positions,values):
                    if 1 <= p < length and mask[p-1]:
                        normalized_values[p-1] = format_reactivity(v/normalizer,threshold)
            else:
                normalized_values = ['NA']*len(mask)
            normalized_values.append('NA')
            data_out[transcript] = normalized_values
        else:
            missing_transcripts[transcript] = None
    return data_out,missing_transcripts

def write_out_reactivity_file(info, outfyle):
    '''Writes out a <.react> file'''
    with open(outfyle,'w') as g:
        for name, values in info.items():
            g.write(name+'\n')
            g.write('\t'.join(values)+'\n')

def write_out_missing_file(info,outfyle):
    '''Writes out the transcripts which were missing from the normalizaition scale.'''
    with open(outfyle,'w') as g:
        for transcript in info.keys():
            g.write(transcript+'\n')

def check_extension(astring,extension):
    '''Checks and fixes things to have the proper extension'''
    out_string = astring if astring.endswith(extension) else astring + extension
    return out_string

def main():
    parser = argparse.ArgumentParser(description='Generates <.react> files from two <.rtsc> files')
    parser.add_argument('control',type=str,help='Control <.rtsc> file')
    parser.add_argument('treatment',type=str,help='Reagent <.rtsc> file')
    parser.add_argument('fasta',type=str,help='Transcript <.fasta> file')
    parser.add_argument('-threshold',type=float,default=7.0,help='[default = 7.0] Reactivity Cap')
    parser.add_argument('-ln_off',action='store_true', help='Do not take the natural log of the stop counts')
    parser.add_argument('-nrm_off',action='store_true',help='Turn off 2-8'+u"％"+' normalization of the derived reactivity')
    parser.add_argument('-save_fails',action='store_true',help='Log transcripts with zero or missing scales')
    parser.add_argument('-scale',type=str,default=None, help='Provide a normalizaiton <.scale> for calculation')
    parser.add_argument('-bases',type=str,default='AC', help='[default = AC] Reaction Specificity, (AGCT) for SHAPE')
    parser.add_argument('-name',type=str,default=None, help='Change the name of the outfile, overrides default')
    parser.add_argument('-restrict',default = None, help = 'Limit analysis to these specific transcripts <.txt> ')
    parser.add_argument('-trim3',type=int,default=0, help='[default = 0] Number of bases to exclude from the 3\' end before normalization')
    args = parser.parse_args()

    #Create output name
    base_name = [x.split(os.sep)[-1].replace('.rtsc','') for x in [args.control,args.treatment]]
    log_tag = ['ln'] if args.ln_off==False else []
    nrm_tag = ['nrm'] if args.ln_off==False else []
    out_name = '_'.join(base_name+log_tag+nrm_tag)+'.react' if args.name == None else check_extension(args.name,'.react')

    #Read in data, apply restrictions if applicable
    covered = get_covered_transcripts(args.restrict) if args.restrict != None else None
    control_data,treatment_data = read_in_rtsc(args.control,covered),read_in_rtsc(args.treatment,covered)

    #Calculate Derived Reactivity
    data = calculate_raw_reactivity(control_data,treatment_data,args.ln_off)

    #Read in transcript sequences
    seqs = read_in_fasta(args.fasta)

    #Generate and write scale, or read a <.scale> file in
    normalizaiton_scale = generate_normalization_scale(data,seqs,args.bases,args.trim3) if args.scale == None else read_normalization_scale(args.scale)
    if args.scale == None:
        write_normalization_scale(normalizaiton_scale,out_name.replace('.react','.scale'))

    #Calculate Final Reactivity
    out_reactivity,out_missing = calculate_final_reactivity(data,seqs,args.bases,args.threshold,normalizaiton_scale,args.nrm_off)

    #Write Out
    write_out_reactivity_file(out_reactivity,out_name)
    if args.save_fails:
        write_out_missing_file(out_missing,out_name.replace('.react','_unresolvable_transcripts.txt'))

if __name__ == '__main__':
    main()
