#!/usr/bin/env python2

'''
Same as StructureFold2's react_heat_correct.py (must be used on <.react>
files already filtered to the same, coverage-qualified transcript set;
supports multiple replicates per temperature condition), except the
correction FACTORS themselves are computed from a stricter position set:
only positions where every provided lower-temperature replicate AND every
provided higher-temperature replicate has a reactivity > 0 (not just a
valid, non-NA value -- exact 0.0s are excluded too). The correction is
still applied to every position of every replicate as usual; only which
positions feed the factor calculation changes.
'''

#Imports
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'StructureFold2'))
from sf2libs.structure_io import read_react, write_react

#Functions
def check_transcript_sets(react_dicts, labels):
    '''Verify all react dicts share the same transcript set'''
    key_sets = [set(d.keys()) for d in react_dicts]
    reference = key_sets[0]
    for i, ks in enumerate(key_sets[1:], 1):
        if ks != reference:
            print 'Warning! Non-parallel transcript sets for sample {}! Quitting...'.format(labels[i])
            sys.exit()

def build_positive_mask(cold_reacts, hot_reacts):
    '''
    For each transcript, returns a list of booleans -- True at positions
    where every lower AND every higher replicate has a float reactivity > 0.
    '''
    mask = {}
    reference = cold_reacts[0]
    all_reacts = cold_reacts + hot_reacts
    for transcript, values in reference.items():
        flags = []
        for i in range(len(values)):
            keep = True
            for react_dict in all_reacts:
                v = react_dict[transcript][i]
                if not (isinstance(v, float) and v > 0):
                    keep = False
                    break
            flags.append(keep)
        mask[transcript] = flags
    return mask

def masked_sum(react_dict, mask):
    '''Sum of only the reactivities at masked (both-temperatures-positive) positions'''
    total = 0.0
    for transcript, flags in mask.items():
        values = react_dict[transcript]
        for i, keep in enumerate(flags):
            if keep:
                total += values[i]
    return total

def apply_correction(react_dict, correction):
    '''Applies a correction'''
    atarashii = {}
    for k, v in react_dict.items():
        atarashii[k] = [x*correction if x != 'NA' else 'NA' for x in v]
    return atarashii

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Corrects <.react>s for differential temperature, using only both-temperatures-positive bases to compute the correction factors')
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

    #Positions where every replicate on both sides has reactivity > 0
    mask = build_positive_mask(cold_reacts, hot_reacts)

    #Compute per-replicate sums over the masked positions only, then average within each condition
    cold_sums = [masked_sum(d, mask) for d in cold_reacts]
    hot_sums  = [masked_sum(d, mask) for d in hot_reacts]
    cold_mean = sum(cold_sums) / float(len(cold_sums))
    hot_mean  = sum(hot_sums)  / float(len(hot_sums))

    #Calculate corrections from condition means
    grand_mean      = (cold_mean + hot_mean) / 2.0
    cold_correction = grand_mean / cold_mean
    heat_correction = grand_mean / hot_mean

    n_positive = sum(1 for flags in mask.values() for keep in flags if keep)
    print '{} positions with reactivity > 0 at both temperatures used for factor calculation'.format(n_positive)
    print 'Higher temp values to be scaled by factor: {}'.format(heat_correction)
    print 'Lower temp values to be scaled by factor: {}'.format(cold_correction)

    #Apply per-replicate corrections (to EVERY position, not just the masked ones) and write out
    for fyle, react in zip(args.lower, cold_reacts):
        corrected = apply_correction(react, cold_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

    for fyle, react in zip(args.higher, hot_reacts):
        corrected = apply_correction(react, heat_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

if __name__ == '__main__':
    main()
