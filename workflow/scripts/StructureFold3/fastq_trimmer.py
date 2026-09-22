#!/usr/bin/env python3
#StructureFold2 batch running script for using cutadapt on <.fastq> before aligning to a reference
#This is an accessory script and is not a core part of StructureFold2, but you may find it useful.
#
#Python 3 port. Without -log the three cutadapt steps (5' adapter, 3' adapter, length/quality) are chained with pipes instead
#of writing two intermediate <.fastq> files, so they run at the same time; every step still sees exactly what it did before,
#so the trimmed <.fastq> is identical. With -log the original step-by-step commands are kept so the logs read the same.

#Imports
import glob
import subprocess
import time
import argparse
import sys

def trim_fastq(afile,fiveprime,threeprime,minlen,minqual,maxlen,log_flag,suffix,nextseqflag):
    '''Runs the cutadapt steps in sequence'''
    out1,out2,out3 = [afile.split('.')[0]+'_'+str(x)+'.fastq' for x in ['1','2',suffix]]
    Q = ['-q',minqual] if nextseqflag == False else ['--nextseq-trim='+str(minqual)]
    M = ['-M',maxlen] if maxlen else []
    if log_flag == False:
        #Pipe the steps together (a step without -o writes to stdout, and '-' reads stdin); the last one writes the trimmed file
        steps = [['cutadapt','-n','2','-g',fiveprime,afile],
                 ['cutadapt','-a',threeprime,'-'],
                 ['cutadapt','-m',minlen]+Q+M+['-o',out3,'-']]
        processes = []
        for i,step in enumerate(steps):
            processes.append(subprocess.Popen(step,stdin=processes[-1].stdout if processes else None,
                                              stdout=subprocess.PIPE if i < len(steps)-1 else None))
            if i > 0:
                processes[-2].stdout.close() #so an earlier step gets SIGPIPE if a later one dies
        for process in processes:
            process.wait()
        for name,process in zip(['5\' adapter','3\' adapter','length/quality'],processes):
            if process.returncode != 0:
                raise RuntimeError('cutadapt ({} step) failed with exit status {} for {}'.format(name,process.returncode,afile))
    if log_flag == True:
        commands = [['cutadapt','-n','2','-g',fiveprime,'-o',out1,afile],
                    ['cutadapt','-a',threeprime,'-o',out2,out1],
                    ['cutadapt','-m',minlen]+Q+M+['-o',out3,out2]]
        logs = []
        for command in commands:
            logs.append(subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True).stdout)
        subprocess.call(['rm',out1])
        subprocess.call(['rm',out2])
        return ''.join(logs)

def dump_bucket(alist,outfyle):
    '''Dumps the logs of stdout to a file for concise record keeping.'''
    with open(outfyle,'w') as g:
        for line in alist:
            g.write(line)

def main():
    parser = argparse.ArgumentParser(description='Batch run cutadapt with given settings on all <.fastq> in a directory')
    parser.add_argument('-log',action="store_true",default=False,help = 'Create an explicit log of the trimming')
    parser.add_argument('-nextseq',action="store_true",default=False,help = 'Use NextSeq/NovaSeq quality scores')
    parser.add_argument('-fp',type=str, default='TGAACAGCGACTAGGCTCTTCA', help='[default = TGAACAGCGACTAGGCTCTTCA] 5\' adapter',dest='fp_adapt')
    parser.add_argument('-tp',type=str, default='GATCGGAAGAGCACACGTCTG', help='[default = GATCGGAAGAGCACACGTCTG] 3\' adapter',dest='tp_adapt')
    parser.add_argument('-minlen',type=str, default='20', help='[default = 20] minimum accepted sequence length', dest='min_len')
    parser.add_argument('-minqual',type=str, default='30', help='[default = 30] minumum accepted base quality',dest='min_qual')
    parser.add_argument('-maxlen',type=str, default=None, help='[default = None] maximum seq length',dest='max_len')#Make this optional
    parser.add_argument('-suffix',type=str,default='trimmed', help='[default = trimmed] trimmed <.fastq> file suffix',dest='suffix')
    parser.add_argument('-logname',type=str,default='trim_log', help='[default = trim_log] Name of the log file')
    args = parser.parse_args()
    #
    if args.log == False:
        print('')
        print('\033[1;4;94mStructure Fold2:\033[0;0;92m trim_fastqs.py\033[0m')
        print('')
        print('\033[94mMass trimming ALL <.fastq> files\033[0m')
        print('')
        start_time = time.asctime()
        fastqs,count = sorted(glob.glob('*.fastq')),0
        for fyle in fastqs:
            trim_fastq(fyle,args.fp_adapt,args.tp_adapt,args.min_len,args.min_qual,args.max_len,args.log,args.suffix,args.nextseq)
            count+=1
        print('\033[1;4;94mStart:\033[0m',start_time)
        print('\033[1;4;94mFinished:\033[0m',time.asctime())
        print('\033[1;4;92mFASTQ Processed:\033[0m',count)
        print('')
    else:
        bucket = []
        fastqs,count = sorted(glob.glob('*.fastq')),0
        for fyle in fastqs:
            start_time = time.asctime()+'\n'
            all_lines = trim_fastq(fyle,args.fp_adapt,args.tp_adapt,args.min_len,args.min_qual,args.max_len,args.log,args.suffix,args.nextseq)
            end_time = time.asctime()+'\n'
            bucket.append(start_time+all_lines+end_time)
        dump_bucket(bucket,args.logname+'.txt')

if __name__ == '__main__':
    main()
