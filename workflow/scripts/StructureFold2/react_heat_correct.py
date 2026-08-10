#!/usr/bin/env python2

'''
Must be used on <.react> files that have already been filtered to only 
include the transcripts above the accepted coverage threshold, and must be the exact same transcripts!
Supports multiple replicates per temperature condition.
'''

#Imports
from sf2libs.structure_io import read_react, write_react
import argparse
import sys

#Functions
def sum_react(react_dict):
    '''Sum of all the reactivities'''
    return sum([sum(filter(lambda x: isinstance(x, float), v)) for v in react_dict.values()])

def apply_correction(react_dict, correction):
    '''Applies a correction'''
    atarashii = {}
    for k, v in react_dict.items():
        atarashii[k] = [x*correction if x != 'NA' else 'NA' for x in v]
    return atarashii

def check_transcript_sets(react_dicts, labels):
    '''Verify all react dicts share the same transcript set'''
    key_sets = [set(d.keys()) for d in react_dicts]
    reference = key_sets[0]
    for i, ks in enumerate(key_sets[1:], 1):
        if ks != reference:
            print 'Warning! Non-parallel transcript sets for sample {}! Quitting...'.format(labels[i])
            sys.exit()

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Corrects <.react>s for differential temperature, supporting multiple replicates per condition')
    parser.add_argument('-lower', type=str, nargs='+', required=True, help='Lower temp <.react> file(s)')
    parser.add_argument('-higher', type=str, nargs='+', required=True, help='Higher temp <.react> file(s)')
    parser.add_argument('-suffix', type=str, default='corrected', help='[default = corrected] Suffix for out files')
    args = parser.parse_args()

    #Read in all replicates
    cold_reacts = [read_react(f) for f in args.lower]
    hot_reacts  = [read_react(f) for f in args.higher]

    #Check transcript sets are parallel across all samples
    all_reacts  = cold_reacts + hot_reacts
    all_labels  = args.lower + args.higher
    check_transcript_sets(all_reacts, all_labels)

    #Compute per-replicate sums, then average within each condition
    cold_sums = map(sum_react, cold_reacts)
    hot_sums  = map(sum_react, hot_reacts)
    cold_mean = sum(cold_sums) / float(len(cold_sums))
    hot_mean  = sum(hot_sums)  / float(len(hot_sums))

    #Calculate corrections from condition means
    grand_mean      = (cold_mean + hot_mean) / 2.0
    cold_correction = grand_mean / cold_mean
    heat_correction = grand_mean / hot_mean

    print 'Higher temp values to be scaled by factor: {}'.format(heat_correction)
    print 'Lower temp values to be scaled by factor: {}'.format(cold_correction)

    #Apply per-replicate corrections and write out
    for fyle, react in zip(args.lower, cold_reacts):
        corrected = apply_correction(react, cold_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

    for fyle, react in zip(args.higher, hot_reacts):
        corrected = apply_correction(react, heat_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

if __name__ == '__main__':
    main()