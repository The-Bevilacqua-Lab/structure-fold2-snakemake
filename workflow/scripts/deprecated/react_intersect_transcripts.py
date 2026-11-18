#!/usr/bin/env python2

'''
Restricts a set of <.react> files to the transcripts present in ALL of
them, and rewrites each one (to a possibly different path) containing only
that shared transcript set. Used to give react_heat_correct.py inputs that
satisfy its "must be the exact same transcripts" requirement, since
rtsc_to_react.py's own per-condition filtering (zero +DMS signal, missing
normalization scale) can otherwise leave two coverage-restricted .react
files with slightly different transcript sets.
'''

import argparse
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'StructureFold2'))
from sf2libs.structure_io import read_react, write_react


def main():
    parser = argparse.ArgumentParser(description='Restrict several <.react> files to their shared transcript set')
    parser.add_argument('-in', dest='infiles', type=str, nargs='+', required=True, help='Input <.react> files')
    parser.add_argument('-out', dest='outfiles', type=str, nargs='+', required=True, help='Output <.react> files, same order/count as -in')
    args = parser.parse_args()

    if len(args.infiles) != len(args.outfiles):
        sys.exit('-in and -out must have the same number of files')

    reacts = [read_react(f) for f in args.infiles]
    shared = set(reacts[0].keys())
    for r in reacts[1:]:
        shared &= set(r.keys())

    print 'Shared transcripts across {} file(s): {}'.format(len(reacts), len(shared))

    for infile, outfile, react in zip(args.infiles, args.outfiles, reacts):
        restricted = {k: v for k, v in react.items() if k in shared}
        write_react(restricted, outfile, sort_flag=True)


if __name__ == '__main__':
    main()
