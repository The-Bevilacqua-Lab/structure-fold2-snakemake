#!/usr/bin/env python3

'''
Must be used on <.react> files that have already been filtered to only 
include the transcripts above the accepted coverage threshold, and must be the exact same transcripts!
Supports multiple replicates per temperature condition.

Python 3 port. The <.react> files are kept as text; reactivities are summed strictly left to right (as Python 2's sum did),
and each distinct value is corrected and formatted once and looked up after that. Numbers are written like Python 2's
str() (12 significant digits), so the corrected files and printed factors match the Python 2 version.
'''

#Imports
import argparse
import sys
from sf3libs import read_react_records, Py2FloatText, py2_str, sequential_sum, py2_dict_order

class FloatCache(dict):
    '''token -> float, parsed once per distinct token'''
    def __missing__(self,token):
        value = self[token] = float(token)
        return value

#Functions
def sum_react(react_dict,floats):
    '''Sum of all the reactivities (NA is skipped): a sequential sum per transcript, then over the transcripts. The
    transcripts are taken in the order Python 2 iterated the dict (float addition is order dependent), which reproduces its sums.'''
    per_transcript = (sequential_sum(map(floats.__getitem__,filter('NA'.__ne__,react_dict[key].split()))) for key in py2_dict_order(react_dict))
    return sequential_sum(per_transcript)

def apply_correction(react_dict, correction):
    '''Applies a correction'''
    text = Py2FloatText(correction)
    return {k: '\t'.join(map(text.__getitem__,v.split())) for k, v in react_dict.items()}

def write_react(react_dictionary,outfile):
    '''Writes out corrected <.react> text'''
    with open(outfile,'w',buffering=1 << 20) as g:
        for transcript, entry in react_dictionary.items():
            g.write(transcript+'\n'+entry+'\n')

def check_transcript_sets(react_dicts, labels):
    '''Verify all react dicts share the same transcript set'''
    key_sets = [set(d.keys()) for d in react_dicts]
    reference = key_sets[0]
    for i, ks in enumerate(key_sets[1:], 1):
        if ks != reference:
            print('Warning! Non-parallel transcript sets for sample {}! Quitting...'.format(labels[i]))
            sys.exit(1)

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Corrects <.react>s for differential temperature, supporting multiple replicates per condition')
    parser.add_argument('-lower', type=str, nargs='+', required=True, help='Lower temp <.react> file(s)')
    parser.add_argument('-higher', type=str, nargs='+', required=True, help='Higher temp <.react> file(s)')
    parser.add_argument('-suffix', type=str, default='corrected', help='[default = corrected] Suffix for out files')
    args = parser.parse_args()

    #Read in all replicates
    cold_reacts = [dict(read_react_records(f)) for f in args.lower]
    hot_reacts  = [dict(read_react_records(f)) for f in args.higher]

    #Check transcript sets are parallel across all samples
    all_reacts  = cold_reacts + hot_reacts
    all_labels  = args.lower + args.higher
    check_transcript_sets(all_reacts, all_labels)

    #Compute per-replicate sums, then average within each condition
    floats = FloatCache()
    cold_sums = [sum_react(r,floats) for r in cold_reacts]
    hot_sums  = [sum_react(r,floats) for r in hot_reacts]
    cold_mean = sequential_sum(cold_sums) / float(len(cold_sums))
    hot_mean  = sequential_sum(hot_sums)  / float(len(hot_sums))

    #Calculate corrections from condition means
    grand_mean      = (cold_mean + hot_mean) / 2.0
    cold_correction = grand_mean / cold_mean
    heat_correction = grand_mean / hot_mean

    print('Higher temp values to be scaled by factor: {}'.format(py2_str(heat_correction)))
    print('Lower temp values to be scaled by factor: {}'.format(py2_str(cold_correction)))

    #Apply per-replicate corrections and write out
    for fyle, react in zip(args.lower, cold_reacts):
        corrected = apply_correction(react, cold_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

    for fyle, react in zip(args.higher, hot_reacts):
        corrected = apply_correction(react, heat_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

if __name__ == '__main__':
    main()
