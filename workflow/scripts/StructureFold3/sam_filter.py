#!/usr/bin/env python3

'''
StructureFold2 batch running script for post-processing <.sam> alignment files before generating <.rtsc> files
This is considered an important part of the StructureFold2 pipleline and should be used even if you chose not to use bowtie2 as your
aligner or used the other batch scripts in previous steps
This script serves two primary purposes
1)Filter out reads in <.sam> files which are either unmapped or improperly mapped
2)Filter out reads in <.sam> files with more than x mismatches, or a mismatch on base position 1 (improper RT stop)

Python 3 port. Does not need samtools: the bitflag filter (samtools view -h -F <flag>) is applied
directly, in the same pass as the mismatch filter, and the file is processed in parallel
byte-range chunks (-threads). No intermediate file is written unless -keep_all is given.
'''

#Imports
import glob
import os
import re
import argparse
from multiprocessing import Pool

#Functions
CHUNK_BYTES = 64 << 20
_NM = re.compile(r'\tNM:i:(-?\d+)').search
_MD = re.compile(r'\tMD:Z:(\d*)').search
_LEADING_INT = re.compile(r'^(\d+)').match

def _scan_all_fields(line,max_mismatch):
    '''Header-line scan: every field is checked for NM:i:/MD:Z: substrings.
    Returns (flag_nmmd,flag1,flag2).'''
    flag1,flag2,flag_nmmd=0,0,0
    for field in line.strip().split('\t'):
        if "NM:i:" in field:
            if int(field.rpartition(':')[2].strip()) > max_mismatch:
                flag1 = 1
            flag_nmmd = 1
        if "MD:Z:" in field:
            match = _LEADING_INT(field.rpartition(':')[2].strip())
            if match and int(match.group(1)) == 0:
                flag2 = 1
            flag_nmmd = 1
    return flag_nmmd,flag1,flag2

def process_chunk(args):
    '''Filters one byte range of a SAM file (which starts and ends on line boundaries).
    Returns (kept_text,tee_text,exceedandfirst,exceedmismatch,firstmismatch,header,lines_in,lines_out,newlines)'''
    path,start,end,bitflag,max_mismatch,allow_bp_first,want_tee = args
    with open(path,'rb') as f:
        f.seek(start)
        raw = f.read(end-start)
    newlines = raw.count(b'\n')
    lines = raw.decode().split('\n')
    if lines and lines[-1] == '':
        lines.pop()

    exceedmismatch,firstmismatch,exceedandfirst,header=0,0,0,0
    lines_in=0
    out,tee=[],[]
    for line in lines:
        if line[:1] == '@':
            header+=1
            lines_in+=1
            if want_tee:
                tee.append(line)
            flag_nmmd,flag1,flag2 = _scan_all_fields(line,max_mismatch)
        else:
            #Same as samtools view -F: drop the read if any bit of bitflag is set
            t1 = line.index('\t')
            t2 = line.index('\t',t1+1)
            if int(line[t1+1:t2]) & bitflag:
                continue
            lines_in+=1
            if want_tee:
                tee.append(line)
            #Determine whether number of mismatches is over limitation
            nm = _NM(line)
            #Determine whether the read has first nt mismatch
            md = _MD(line)
            if nm is None and md is None:
                continue
            flag_nmmd = 1
            flag1 = 1 if nm is not None and int(nm.group(1)) > max_mismatch else 0
            flag2 = 1 if md is not None and md.group(1) != '' and int(md.group(1)) == 0 else 0

        if flag_nmmd == 1:
            # If don't allow 1st nt mismatch
            if not allow_bp_first:
                if flag1 == 1 and flag2 == 1:
                    exceedandfirst+=1
                elif flag1 == 1 and flag2 == 0:
                    exceedmismatch+=1 #number of reads with mismatches is over limitation excluding those with 1st mismatch
                elif flag1 == 0 and flag2 == 1:
                    firstmismatch+=1  #number of reads with 1st nt mismatch excluding those with mismatches over limitation
            else:
                #If allows 1st nt mismatch
                if flag1 == 1:
                    exceedmismatch+=1

            if flag1 == 0 and (allow_bp_first or flag2 == 0):
                out.append(line)

    kept = '\n'.join(out)+'\n' if out else ''
    tee_text = '\n'.join(tee)+'\n' if tee else ''
    return kept,tee_text,exceedandfirst,exceedmismatch,firstmismatch,header,lines_in,len(out),newlines

def chunk_ranges(path,chunk_bytes=CHUNK_BYTES):
    '''Splits a file into (start,end) byte ranges that begin on a line boundary'''
    size = os.path.getsize(path)
    bounds = [0]
    with open(path,'rb') as f:
        pos = chunk_bytes
        while pos < size:
            f.seek(pos)
            f.readline() #advance to the start of the next line
            nxt = f.tell()
            if nxt >= size:
                break
            if nxt > bounds[-1]:
                bounds.append(nxt)
            pos = nxt + chunk_bytes
    bounds.append(size)
    return list(zip(bounds[:-1],bounds[1:]))

def filter_sam_file(samfyle,bitflag,out_suffix='filtered',max_mismatch=3,allow_bp_first=False,keep_all=False,threads=1):
    '''Applies the bitflag filter (samtools view -h -F <bitflag>) and the mismatch filter to a SAM file.
    Returns (intermediate_sam,filtered_sam,exceedandfirst,exceedmismatch,firstmismatch,header,input_lines,
    lines_after_bitflag,filtered_lines)'''
    intermediate = samfyle.replace('.sam','_samtools.sam')
    outfile = intermediate.replace('samtools.sam',out_suffix+'.sam')
    ranges = chunk_ranges(samfyle)
    jobs = [(samfyle,s,e,int(bitflag),max_mismatch,allow_bp_first,keep_all) for s,e in ranges]

    totals = [0]*7 #exceedandfirst,exceedmismatch,firstmismatch,header,lines_in,lines_out,newlines
    tee = open(intermediate,'w',buffering=1 << 20) if keep_all else None
    pool = Pool(threads) if threads > 1 and len(jobs) > 1 else None
    try:
        with open(outfile,'w',buffering=1 << 20) as g:
            results = pool.imap(process_chunk,jobs) if pool else map(process_chunk,jobs)
            for kept,tee_text,*counts in results:
                g.write(kept)
                if tee is not None:
                    tee.write(tee_text)
                totals = [a+b for a,b in zip(totals,counts)]
    finally:
        if tee is not None:
            tee.close()
        if pool:
            pool.close()
            pool.join()
    exceedandfirst,exceedmismatch,firstmismatch,header,lines_in,lines_out,newlines = totals
    return intermediate,outfile,exceedandfirst,exceedmismatch,firstmismatch,header,newlines,lines_in,lines_out

def write_log_file(info_bucket,outfile):
    '''Takes the info collected from stout and makes a nice consistent <.csv>'''
    sub_buckets = [info_bucket[x:x+12] for x in range(0,len(info_bucket),12)]
    x_header = ','.join(['in_sam','sam_lines','filtered_lines','sam_filter_flag','flag_options',
                         'mismatches_and_first_mismatch','mismatches','first_mismatch',
                         'header_lines','max_mismatch','out_sam','out_sam_lines'])
    with open(outfile,'w') as g:
        g.write(x_header+'\n')
        for record in sub_buckets:
            out_fields = ','.join(str(x) for x in record)
            g.write(out_fields+'\n')

def main():
    parser = argparse.ArgumentParser(description='\033[1;4;94mBatch filter all <.sam> files within a given directory. This will make them ready to be converted into <.rtsc>\033[0m')
    parser.add_argument('-turbo',action='store_true',default=False,help='All filter options ignored, default settings, no log')
    parser.add_argument('-sam',default=None, nargs='+', help='Specific files to operate on')
    parser.add_argument('-keep_all',action='store_true',default=False,help='Keep all intermediate files')
    parser.add_argument('-keep_reverse',action='store_true',default=False,help='Keep mappings from the reverse strand')
    parser.add_argument('-remove_secondary',action='store_true',default=False,help='Remove secondary alignments')
    parser.add_argument('-allow_bp1_mismatch',action='store_true',default=False,help='Accept mappings with first base mismatches',dest='firstmm')
    parser.add_argument('-logname',type=str,default='filter_log.csv', help='[default = filter_log.csv] Name of the log file')
    parser.add_argument('-suffix',type=str,default='filtered', help='[default = filtered] filtered <.sam> file suffix')
    parser.add_argument('-max_mismatch',type=int, default=3, help='[default = 3] Maximum allowed mismatches/indels',dest='max_mm')
    parser.add_argument('-threads',type=int, default=1, help='[default = 1] Worker processes used per file')
    args = parser.parse_args()

    #Make list of files to process
    fyle_lyst = sorted(glob.glob('*.sam')) if args.sam == None else sorted(args.sam)

    #Turbo Mode, no logs, no options.
    if args.turbo == True:
        for fyle in fyle_lyst:
            filter_sam_file(fyle,20,args.suffix,keep_all=args.keep_all,threads=args.threads)

    #Detailed Mode, explicit log.
    else:
        #Assemble a bitflag to use for all processing.
        unmapped = 4 #This currently must be set to 4.
        reverse = 16 if args.keep_reverse == False else 0
        secondary = 256 if args.remove_secondary == True else 0
        descriptions = {4:'unmapped',16:'reverse_strand',256:'secondary_alignments'}
        bitflag = str(sum([unmapped,reverse,secondary]))
        options_column = ' '.join([descriptions[flag] for flag in [x for x in [unmapped,reverse,secondary] if x!=0]])

        #Make a list to fill with information
        bucket = []
        for fyle in fyle_lyst:
            fyle_1,fyle_2,firstrsnp,snp,first,headlines,count_1,count_2,lines_out = filter_sam_file(
                fyle,bitflag,args.suffix,args.max_mm,args.firstmm,args.keep_all,args.threads)

            alignments_after = lines_out - headlines

            bucket.extend([fyle,count_1,count_2,bitflag,options_column])
            bucket.extend([firstrsnp,snp,first,headlines,str(args.max_mm),fyle_2,str(alignments_after)])

        write_log_file(bucket,args.logname)

if __name__ == '__main__':
    main()
