#!/usr/bin/env python3

'''
Restricts a set of <.react> files to the transcripts present in ALL of
them, and rewrites each one (to a possibly different path) containing only
that shared transcript set. Used to give react_heat_correct.py inputs that
satisfy its "must be the exact same transcripts" requirement, since
rtsc_to_react.py's own per-condition filtering (zero +DMS signal, missing
normalization scale) can otherwise leave two coverage-restricted .react
files with slightly different transcript sets.

Python 3 port. The files are kept as text and only re-formatted on the way out (each distinct value is formatted once and
looked up after that), with numbers written like Python 2's str() (12 significant digits) as before.
'''

import argparse
import sys
from sf3libs import read_react_records, Py2FloatText


def main():
    parser = argparse.ArgumentParser(description='Restrict several <.react> files to their shared transcript set')
    parser.add_argument('-in', dest='infiles', type=str, nargs='+', required=True, help='Input <.react> files')
    parser.add_argument('-out', dest='outfiles', type=str, nargs='+', required=True, help='Output <.react> files, same order/count as -in')
    args = parser.parse_args()

    if len(args.infiles) != len(args.outfiles):
        sys.exit('-in and -out must have the same number of files')

    reacts = [dict(read_react_records(f)) for f in args.infiles]
    shared = set(reacts[0].keys())
    for r in reacts[1:]:
        shared &= set(r.keys())

    print('Shared transcripts across {} file(s): {}'.format(len(reacts), len(shared)))

    text = Py2FloatText()
    for outfile, react in zip(args.outfiles, reacts):
        with open(outfile,'w',buffering=1 << 20) as g:
            for transcript in sorted(shared):
                g.write(transcript+'\n'+'\t'.join(map(text.__getitem__,react[transcript].split()))+'\n')


if __name__ == '__main__':
    main()
