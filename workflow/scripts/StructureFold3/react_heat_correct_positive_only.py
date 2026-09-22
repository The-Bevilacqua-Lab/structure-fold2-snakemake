#!/usr/bin/env python3

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

Python 3 port (see react_heat_correct.py): text-based, sequential sums, Python 2 style number formatting.
'''

#Imports
import argparse
import sys
import numpy as np
from sf3libs import read_react_records, py2_str, sequential_sum, py2_dict_order
from react_heat_correct import FloatCache, apply_correction, write_react, check_transcript_sets

#Functions
class PositiveCache(dict):
    '''token -> float, with NA as 0.0 so that it never counts as > 0'''
    def __missing__(self,token):
        value = self[token] = 0.0 if token == 'NA' else float(token)
        return value

def build_positive_mask(cold_reacts, hot_reacts):
    '''
    For each transcript, returns (float arrays per replicate, boolean array) -- True at positions
    where every lower AND every higher replicate has a float reactivity > 0. The mask dict is filled in the order Python 2
    iterated the reference dict, as the original did, so that iterating it later matches the Python 2 summation order.
    '''
    mask,values = {},{}
    parse = PositiveCache().__getitem__
    all_reacts = cold_reacts + hot_reacts
    for transcript in py2_dict_order(cold_reacts[0]):
        reference = cold_reacts[0][transcript]
        length = len(reference.split())
        arrays = []
        for react_dict in all_reacts:
            tokens = react_dict[transcript].split()
            if len(tokens) < length:
                raise IndexError('list index out of range')
            arrays.append(np.array(list(map(parse,tokens[:length])),dtype=float))
        values[transcript] = arrays
        flags = arrays[0] > 0
        for array in arrays[1:]:
            flags &= array > 0
        mask[transcript] = flags
    return mask,values

def masked_sum(index, mask, values):
    '''Sum of only the reactivities at masked (both-temperatures-positive) positions, added strictly in order'''
    total = 0.0
    for transcript in py2_dict_order(mask):
        total = sequential_sum(values[transcript][index][mask[transcript]].tolist(),total)
    return total

#Main Function
def main():
    parser = argparse.ArgumentParser(description='Corrects <.react>s for differential temperature, using only both-temperatures-positive bases to compute the correction factors')
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

    #Positions where every replicate on both sides has reactivity > 0
    mask,values = build_positive_mask(cold_reacts, hot_reacts)

    #Compute per-replicate sums over the masked positions only, then average within each condition
    n_cold = len(cold_reacts)
    cold_sums = [masked_sum(i, mask, values) for i in range(n_cold)]
    hot_sums  = [masked_sum(n_cold+i, mask, values) for i in range(len(hot_reacts))]
    cold_mean = sequential_sum(cold_sums) / float(len(cold_sums))
    hot_mean  = sequential_sum(hot_sums)  / float(len(hot_sums))

    #Calculate corrections from condition means
    grand_mean      = (cold_mean + hot_mean) / 2.0
    cold_correction = grand_mean / cold_mean
    heat_correction = grand_mean / hot_mean

    n_positive = sum(int(flags.sum()) for flags in mask.values())
    print('{} positions with reactivity > 0 at both temperatures used for factor calculation'.format(n_positive))
    print('Higher temp values to be scaled by factor: {}'.format(py2_str(heat_correction)))
    print('Lower temp values to be scaled by factor: {}'.format(py2_str(cold_correction)))

    #Apply per-replicate corrections (to EVERY position, not just the masked ones) and write out
    for fyle, react in zip(args.lower, cold_reacts):
        corrected = apply_correction(react, cold_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

    for fyle, react in zip(args.higher, hot_reacts):
        corrected = apply_correction(react, heat_correction)
        write_react(corrected, fyle.replace('.react', '_'+args.suffix+'.react'))

if __name__ == '__main__':
    main()
