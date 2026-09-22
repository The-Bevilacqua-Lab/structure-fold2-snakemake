'''
Small helpers shared by the StructureFold3 (Python 3) scripts.

These replace StructureFold2's sf2libs for the parts the pipeline uses. Two things are deliberate so that results match the
Python 2 versions: floats are written like Python 2's str() (py2_str), and float sums are always added strictly left to right
(sequential_sum), because Python 3.12+'s sum() compensates float sums and would change the last bits.
'''

import re
from functools import lru_cache, reduce
from operator import add

def py2_str(value):
    '''str() of a float as Python 2 wrote it: 12 significant digits, and a trailing .0 on whole numbers'''
    text = '%.12g' % value
    return text + '.0' if text.lstrip('-').isdigit() else text

class Py2FloatText(dict):
    '''Cache mapping a numeric token to py2_str(float(token) * factor); 'NA' stays 'NA'. Reactivity files have few distinct
    values, so formatting each distinct token once and looking the rest up is much faster than formatting every one.'''
    def __init__(self,factor=None):
        super().__init__()
        self.factor = factor
        self['NA'] = 'NA'
    def __missing__(self,token):
        value = float(token)
        text = self[token] = py2_str(value if self.factor is None else value*self.factor)
        return text

def sequential_sum(values,start=0):
    '''sum() that always adds strictly left to right, starting from the int 0 (like Python 2)'''
    return reduce(add,values,start)

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

def read_rtsc_records(rtsc_fyle,keep=None):
    '''Streams a <.rtsc> file as (transcript_name, tab-separated stops text); each record is a transcript line, a counts
    line and a blank line. keep optionally limits which transcripts are yielded.'''
    with open(rtsc_fyle,'r') as f:
        lines = (line.strip() for line in f)
        lines = (line for line in lines if line)
        for transcript,stops in zip(lines,lines):
            if keep is None or transcript in keep:
                yield transcript,stops

def read_react_records(react_fyle):
    '''Streams a <.react> file as (transcript_name, whitespace-separated reactivities text); two lines per record'''
    with open(react_fyle,'r') as f:
        lines = (line.strip() for line in f)
        yield from zip(lines,lines)

def read_restrict(txt_file):
    '''Reads in a standard overlapped coverage file (one transcript per line)'''
    info = {}
    with open(txt_file,'r') as f:
        for line in f:
            info[line.strip()] = None
    return info

_NONZERO_DIGIT = re.compile('[1-9]').finditer

def nonzero_counts(stops):
    '''Returns (number of positions, [positions with a non-zero count], [those counts as ints]) for a tab-separated line.
    Counts are overwhelmingly 0, so rather than splitting every token this finds the non-zero digits and converts just
    the tokens containing them.'''
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
        values.append(int(stops[token_start:token_end]))
    return stops.count('\t')+1,positions,values

@lru_cache(maxsize=None)
def specificity_table(specificity):
    '''Translation table mapping the specificity bases to 1 and every other byte to 0'''
    return bytes(1 if chr(i) in specificity else 0 for i in range(256))

def specificity_mask(sequence,specificity):
    '''bytes of 1 where sequence[i] is one of the specificity bases (else 0)'''
    return sequence.encode('ascii','replace').translate(specificity_table(specificity))

def py2_dict_order(keys):
    '''The order in which Python 2.7 iterates a dict of str keys that were inserted in the given order (with hash
    randomization off, its default). Python 2's sums of per-transcript values followed this order, and float addition is not
    associative, so summing in the same order reproduces its results to the last bit. Emulates CPython 2.7's string hash,
    open-addressing probe and resize rules.'''
    mask64 = (1 << 64) - 1
    def py2_hash(text):
        data = text.encode('utf-8')
        if not data:
            return 0
        x = data[0] << 7
        for byte in data:
            x = ((1000003 * x) ^ byte) & mask64
        x ^= len(data)
        return mask64 - 1 if x == mask64 else x #a hash of -1 is stored as -2
    def place(table,hash_value,key):
        mask = len(table) - 1
        i,perturb = hash_value & mask,hash_value
        while table[i] is not None:
            if table[i][0] == key:
                return False
            i = (5 * i + perturb + 1) & mask64
            perturb >>= 5
            i &= mask
        table[i] = (key,hash_value)
        return True
    table,used = [None]*8,0
    for key in keys:
        if place(table,py2_hash(key),key):
            used += 1
            if used * 3 >= len(table) * 2:
                #Resize: rebuild in old slot order, into a table larger than 4x (2x above 50000) the entries
                size = 8
                while size <= (2 if used > 50000 else 4) * used:
                    size <<= 1
                new_table = [None]*size
                for entry in table:
                    if entry is not None:
                        place(new_table,entry[1],entry[0])
                table = new_table
    return [entry[0] for entry in table if entry is not None]
